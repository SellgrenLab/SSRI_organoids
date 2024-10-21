########### Pseudobulk analysis ###########

DefaultAssay(project1) <- "RNA"
results <- list()

for(i in 1:length(celltype)){ 
    sub <- subset(project1, subset = types == celltype[i])
    DefaultAssay(sub) <- "RNA"
    sub_df <- as.data.frame(sub@meta.data)
    sub_sampled <- sub_df %>%
        group_by(Condition, assignment) %>%
        sample_n(size = min(100, n()), replace = FALSE) %>%
        ungroup()
    
    # Ensure equal number of cells in each condition
    min_cells <- min(table(sub_sampled$Condition))
    sub_sampled <- sub_sampled %>%
        group_by(Condition) %>%
        sample_n(size = min_cells, replace = FALSE) %>%
        ungroup()
    
    sub <- subset(sub, subset= barcode%in%sub_sampled$barcode)
    results[[i]] <- DElegate::findDE(sub, group_column = "Condition", method = "deseq", replicate_column = "assignment",)
    results[[i]]$celltype <- celltype[i]
}

names(results) <- celltype
results_df <- do.call(rbind, results)

write.xlsx(results_df, "./sertraline/201024_SSRI_pseudobulk_results_all.xlsx")


### Plot an upset plot to show the overlap of DE genes across cell types
results_df_sig <- dplyr::filter(results_df, padj < 0.05)

library(UpSetR)

# Create a list of DE genes for each cell type
de_genes_list <- split(results_df_sig$feature, results_df_sig$celltype)

# Create an UpSet plot
upset(fromList(de_genes_list), order.by = "freq", main.bar.color = "dodgerblue", sets.bar.color = "darkorange", text.scale = 1.5, nsets = 11)



########### Pathway analysis ###########

results_df <- read.xlsx("./sertraline/201024_SSRI_pseudobulk_results_all.xlsx")

results_pathway <- list()

for(i in 1:length(celltype)){
  df = dplyr::filter(results_df, celltype == celltype[i])
   #subset to include only GO:BP
   geneSets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
   ### filter background to only include genes that we captured.
   geneSets$gene_symbol <- toupper(geneSets$gene_symbol)
   geneSets <- geneSets[geneSets$gene_symbol %in% df$feature,]
   m_list <- geneSets %>% split(x = .$gene_symbol, f = .$gs_name)
   stats <- df$stat
   names(stats) <- df$feature
   stats <- na.omit(stats)
  results_pathway[[i]] <- fgsea(pathways = m_list, stats = stats, nperm = 5e4, minSize = 10)
  results_pathway[[i]]<- results_pathway[[i]] %>% filter(padj<0.01) %>% arrange(desc(NES)) 
}


names(results_pathway) <- celltype
#for(i in 1:length(results_pathway)){ results_pathway[[i]]$celltype <- celltype[i]  }

results_pathway_df <- do.call(rbind, results_pathway)
write.xlsx(results_pathway_df, "./sertraline/FGSEA_pathway_results.xlsx")



# Create a new workbook
wb <- createWorkbook()
# Loop through each element in the results list
for (i in 1:length(results_pathway)) {
  pathway <- results_pathway[[i]]  
  # Create a new sheet in the workbook
  addWorksheet(wb, sheetName = celltype[i])
  # Write the pathway data to the sheet
  writeData(wb, sheet = i, x = pathway)
}

# Save the workbook as an Excel file
saveWorkbook(wb, file = "SSRI_pathway_results.xlsx")

## Overrepresentation analysis
library(enrichR)    
# Define the databases to use for enrichment analysis
dbs <- c("GO_Biological_Process_2021", "KEGG_2021_Human")

results_ORA <- list()
for (i in 1:length(results)) {

    # Separate upregulated and downregulated genes
    upregulated_genes <- (results[[i]] %>% filter(padj<0.05 & log_fc > 0))$feature
    downregulated_genes <- (results[[i]] %>% filter(padj<0.05 & log_fc < 0))$feature
    
    # Perform enrichment analysis for upregulated genes
    if (length(upregulated_genes) > 0) {
        enriched_up <- enrichr(upregulated_genes, dbs)
        enriched_up <- do.call(rbind, enriched_up)
        enriched_up$direction <- "Upregulated"
        enriched_up$celltype <- celltype[i]
        results_ORA[[paste0(celltype[i], "_up")]] <- enriched_up
    }
    
    # Perform enrichment analysis for downregulated genes
    if (length(downregulated_genes) > 0) {
        enriched_down <- enrichr(downregulated_genes, dbs)
        enriched_down <- do.call(rbind, enriched_down)
        enriched_down$direction <- "Downregulated"
        enriched_down$celltype <- celltype[i]
        results_ORA[[paste0(celltype[i], "_down")]] <- enriched_down
    }
}

# Combine all results into a single data frame
results_ORA_df <- do.call(rbind, results_ORA)

# Save the results to an Excel file
write.xlsx(results_ORA_df, "./sertraline/EnrichR_pathway_results.xlsx")

# Plot top pathways
top_pathways <- results_ORA_df %>%
    group_by(celltype, direction) %>%
    top_n(10, wt = Combined.Score) %>%
    ungroup()

ggplot(top_pathways, aes(x = reorder(Term, Combined.Score), y = Combined.Score, fill = direction)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    facet_wrap(~ celltype, scales = "free_y") +
    theme_minimal() +
    labs(title = "Top Enriched Pathways", x = "Pathway", y = "Combined Score") +
    scale_fill_manual(values = c("Upregulated" = "firebrick", "Downregulated" = "dodgerblue"))
