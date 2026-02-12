library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(DESeq2)
library(RColorBrewer)
set.seed(123)  # For reproducibility

ssri <- readRDS("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/Internal data/CMV & SSRI/scRNAseq/SSRI/SSRI/sertraline/201024_SSRI_final_clusters.rds")
ssri$sex <- ifelse(ssri$assignment %in% c("0","2","3"), "M", "F")
ssri$barcode <- rownames(ssri@meta.data)
Idents(ssri) <- "types"
celltypes <- levels(ssri)
celltypes <- celltypes[-c(7,11,13)]

results_list_sex <- list()

perform_DE_sex_analysis <- function(cluster) {
  # Subset the Seurat object for the current cluster
  cluster_cells <- WhichCells(ssri, idents = cluster)
  cluster_data <- subset(ssri, cells = cluster_cells)

  sub_df <- as.data.frame(cluster_data@meta.data)
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
    
    cluster_data <- subset(cluster_data, subset= barcode%in%sub_sampled$barcode)
  
  print(paste0("Running cluster ID ", cluster))
  # Aggregate counts by replicate
  cluster_counts <- AggregateExpression(cluster_data, group.by = c("Condition", "assignment"), assays = "RNA", slot = "counts")$RNA
  
  # Create metadata for DESeq2
  cluster_metadata <- cluster_data@meta.data %>%
    group_by(assignment) %>%
    summarize(group_id = unique(Condition), sex = unique(sex))
  
  # Ensure row names of metadata match column names of counts
  cluster_metadata<- as.data.frame(cluster_metadata)
  rownames(cluster_metadata) <- paste0(cluster_metadata$group_id, "_", cluster_metadata$assignment)
  cluster_metadata <- cluster_metadata[ , -1]
  
  # Reorder cluster_metadata to match the order of cluster_counts columns
  cluster_metadata <- cluster_metadata[colnames(cluster_counts), ]
  
  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                colData = cluster_metadata, 
                                design = ~ sex)
  
  # Perform differential expression analysis
  dds <- DESeq(dds, test="LRT", reduced = ~ 1, useT=TRUE, minmu = 1e-6, minReplicatesForReplace=Inf, fitType = "glmGamPoi")
  
  # Extract results for the comparison of interest
  res <- results(dds)
  
  res <- res %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>%
    as_tibble()
  res$cluster_id <- cluster
  
  sigLRT_genes <- res %>% 
    filter(padj < 0.05)
  
  # Save results to a CSV file
  #write.csv(as.data.frame(res), file = paste0("All_pseudobulk_DE_results_cluster_", cluster, ".csv"))
  write.csv(as.data.frame(sigLRT_genes), file = paste0("Significant_sex_specific_DE_results_cluster_", cluster, ".csv"))
  return(res)
}

results_list_sex <- map(celltypes, perform_DE_sex_analysis)

## Get sex-specific genes found across all clusters

extract_significant_genes <- function(res) {
    res %>% 
        filter(padj < 0.05) %>% 
        pull(gene)
}
# Get the list of significant genes for each cluster
significant_genes_list <- lapply(results_list_sex, extract_significant_genes)
all_significant_genes <- unlist(significant_genes_list)

# Get the frequency of each gene
gene_frequency <- table(all_significant_genes)

# Convert to a data frame for easier viewing
gene_frequency_df <- as.data.frame(gene_frequency)
colnames(gene_frequency_df) <- c("gene", "frequency")

# Print the gene frequencies
print(gene_frequency_df)
sex_specific_genes <- as.character(gene_frequency_df$gene[which(gene_frequency_df$frequency>2)]) # no significant se-specific genes identified



############ Main DE analysis #############
results_list <- list()
perform_DE_analysis <- function(cluster) {
    # Subset the Seurat object for the current cluster
  cluster_cells <- WhichCells(ssri, idents = cluster)
  cluster_data <- subset(ssri, cells = cluster_cells)

  sub_df <- as.data.frame(cluster_data@meta.data)
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
    
    cluster_data <- subset(cluster_data, subset= barcode%in%sub_sampled$barcode)
  
    print(paste0("Running cluster ID ", cluster))
    
    # Extract raw counts
    raw_counts <- GetAssayData(cluster_data, assay="RNA", layer = "counts")
    
    raw_counts <- raw_counts[!rownames(raw_counts) %in% c("MALAT1"), ]
    
    # Create a new Seurat object with the filtered counts
    new_seurat_object <- CreateSeuratObject(counts = raw_counts, meta.data = cluster_data@meta.data)
    
    # Aggregate counts by replicate
    cluster_counts <- AggregateExpression(new_seurat_object, group.by = c("Condition", "assignment"), assays = "RNA", slot = "counts")$RNA
    
    # Create metadata for DESeq2
    cluster_metadata <- new_seurat_object@meta.data %>%
        group_by(assignment) %>%
        summarize(group_id = unique(Condition), sex = unique(sex))
    
    # Ensure row names of metadata match column names of counts
    cluster_metadata <- as.data.frame(cluster_metadata)
    rownames(cluster_metadata) <- paste0(cluster_metadata$group_id, "_", cluster_metadata$assignment)
    cluster_metadata <- cluster_metadata[ , -1]
    
    # Reorder cluster_metadata to match the order of cluster_counts columns
    cluster_metadata <- cluster_metadata[colnames(cluster_counts), ]
    
    # Create DESeq2 dataset
    dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                  colData = cluster_metadata, 
                                  design = ~ group_id)
    
    # Perform differential expression analysis
    dds <- DESeq(dds, test="LRT", reduced = ~ 1, useT=TRUE, minmu = 1e-6, minReplicatesForReplace=Inf, fitType = "glmGamPoi")
    
    # Extract results for the comparison of interest
    res <- results(dds, name="group_idlit")
    
    res <- res %>%
        data.frame() %>%
        rownames_to_column(var="gene") %>%
        as_tibble()
    res$cluster_id <- cluster
    
    sigLRT_genes <- res %>% 
        filter(padj < 0.05)
    
    # Save results to a CSV file
    write.csv(as.data.frame(res), file = paste0("./DESeq2/lrt/All_pseudobulk_DE_results_cluster_", cluster, ".csv"))
    write.csv(as.data.frame(sigLRT_genes), file = paste0("./DESeq2/lrt/Significant_pseudobulk_DE_results_cluster_", cluster, ".csv"))
    
    return(res)
}

results_list <- map(celltypes, perform_DE_analysis)

# Combine results into a single data frame
combined_results <- bind_rows(results_list, .id = "cluster")

write.xlsx(combined_results, "./DESeq2/lrt/Sert_All_DE_combined_results.xlsx")

combined_results_sig <- dplyr::filter(combined_results, padj<0.05)
write.xlsx(combined_results_sig, "./DESeq2/lrt/Sert_Significant_DE_combined_results.xlsx")
