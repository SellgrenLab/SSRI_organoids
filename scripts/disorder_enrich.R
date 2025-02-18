library(GeneOverlap)
library(ggplot2)
library(dplyr)

Geneset_with_disease_riskgenes <- readRDS("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/Internal data/C1Q/data/Geneset_with_disease_riskgenes.rds")
Trubetskoy_scz_genes <- readRDS("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/Internal data/C1Q/Bipolar/Trubetskoy_scz_genes.rds")
asd<-  readRDS("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/Internal data/C1Q/Bipolar/ASD_riskgenes.rds") #https://doi-org.proxy.kib.ki.se/10.1073/pnas.2215632120
adhd<-  readRDS("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/Internal data/C1Q/Bipolar/ADHD_riskgenes.rds") #demontis et al
mdd<- read.xlsx("/Users/susmita.malwade/MDD.xlsx") #https://www-nature-com.proxy.kib.ki.se/articles/s41588-023-01596-4

## Read in the pseudo time-related DEGenes
pde <- read.xlsx("./sertraline/Results/Sertraline_PseudoDE_Results.xlsx")
pde<- dplyr::filter(pde, padj<0.01)
genes <- pde$gene

#get genome size from sc dataset (all features that were tested) and test the overlap statistically

print(testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = genes, listB = asd))) ## New asd genes with de novo variants as well
print(testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = genes, listB = Geneset_with_disease_riskgenes$ASD_Satterstrom$Gene)))
print(testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = genes, listB = Geneset_with_disease_riskgenes$ASD_Ruzzo$`HGNC gene symbol`)))
print(testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = genes, listB = Geneset_with_disease_riskgenes$ASD_Sanders$RefSeqGeneName)))
print(testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = genes, listB = Geneset_with_disease_riskgenes$DDD$symbol)))
print(testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = genes, listB = adhd)))


## Create a data frame for plotting
df<- data.frame("Disorder" = c("ASD(Cirnigliaro)", "ASD(Satterstom)", "ASD(Ruzzo)", "ASD (Sanders)", "DDD", "ADHD", "MDD"))
df$pval <- c(2.4e-03, 0.31,0.025, 1.6e-03,0.019,1, 0.25)
df$OR <- c(12.1, 2.8, 8.5, 13.9, 4.1, 0, 2.1)
df$padj <- p.adjust(df$pval, method = "BH")
df$plot <- -log10(df$padj)
## Plot the results

plot1<- ggplot(data= df, mapping = aes(x=plot, y=Disorder, fill=OR)) + geom_bar(stat="identity") + scale_fill_viridis() + xlab(label = "-log(p-adjusted)") + theme_classic() + geom_vline(xintercept = -log10(0.05)) + ggtitle(label = "Disorder risk enrichment in PseudoDEG")


### DEGs
## read in the celltype-specific DEGs
results_de <- read.xlsx("./sertraline/201024_SSRI_pseudobulk_results_all.xlsx")
results_df_sig <- dplyr::filter(results_de, padj < 0.05)

print(testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = mdd$Gene.Name)))
print(testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = adhd)))
print(testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = gene_vector)))
print(testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = Geneset_with_disease_riskgenes$ASD_Ruzzo$`HGNC gene symbol`)))
print(testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = Geneset_with_disease_riskgenes$ASD_Sanders$RefSeqGeneName)))
print(testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = Geneset_with_disease_riskgenes$ASD_Satterstrom$Gene)))
print(testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = Geneset_with_disease_riskgenes$DDD$symbol)))


df1<- data.frame("Disorder" = c("ASD(Cirnigliaro)", "ASD(Satterstom)", "ASD(Ruzzo)", "ASD (Sanders)", "DDD", "ADHD", "MDD"))
df1$pval <- c(1.1e-05,8.6e-05,7.4e-05,5.3e-05,8.3e-06,0.26,1.9e-03)
df1$OR <- c(10.4, 7.3,9.4, 10.1,4.8,3.5, 3.4)
df1$padj <- p.adjust(df1$pval, method = "BH")
df1$plot <- -log10(df1$padj)

## Plot the results

plot2<- ggplot(data= df1, mapping = aes(x=plot, y=Disorder, fill=OR)) + geom_bar(stat="identity") + scale_fill_viridis() + xlab(label = "-log(p-adjusted)") + theme_classic() + geom_vline(xintercept = -log10(0.05)) + ggtitle(label = "Disorder risk enrichment in DEGs")


cowplot::plot_grid(plot1,plot2,ncol = 1)



#### p cutoff set to 0.01 for degs as well
results_df_sig <- dplyr::filter(results_de, padj < 0.01)

df1$pval <- c(0.04,9.8e-04, 0.26, 2.8e-03,4e-06,0.12, 2.2e-05,0.015 )
df1$OR <- c(6.6,9.8, 3.4, 11.5,8.1,8.4, 7.6, 6.1 )
df1$padj <- p.adjust(df1$pval, method = "BH")
df1$plot <- -log10(df1$padj)

ggplot(data = df1, mapping = aes(x = plot, y = Disorder, fill = OR)) +
    geom_bar(stat = "identity") +
    scale_fill_viridis() +
    xlab(label = "-log(p-adjusted)") +
    theme_classic() +
    geom_vline(xintercept = -log10(0.01)) +
    ggtitle(label = "Disorder risk enrichment in DEGs (p<0.01)")






#### New plotting
bd<-readRDS("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/Internal data/C1Q/Bipolar/Bipolar_risk_genes_with BD1.rds")

hic <- Geneset_with_disease_riskgenes$SCZ_HiC$`High-confident Hi-C defined Schizophrenia risk genes`
schema <- c("SET1DA", "CUL1", "XPO7", "TRIO", "CACNA1G", "SP4", "GRIA3", "GRIN2A", "HERC1", "RB1CC1")
cross_disorder<- read.xlsx("/Users/susmita.malwade/Crossdisorder.xlsx")
## Genes from the DE analysis
results_de <- read.xlsx("./201024_SSRI_pseudobulk_results_all.xlsx")
results_df_sig <- dplyr::filter(results_de, padj < 0.05)

pvalues <- c(
    testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = mdd$Gene.Name))@pval,
    testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = adhd))@pval,
    testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = asd))@pval,
    testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = Geneset_with_disease_riskgenes$ASD_Ruzzo$`HGNC gene symbol`))@pval,
    testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = Geneset_with_disease_riskgenes$ASD_Sanders$RefSeqGeneName))@pval,
    testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = Geneset_with_disease_riskgenes$ASD_Satterstrom$Gene))@pval,
    #testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = Geneset_with_disease_riskgenes$DDD$symbol))@pval,
    testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = Trubetskoy_scz_genes$prioritized_genes_trubetskoy))@pval,
    testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = Trubetskoy_scz_genes$all_scz_trubetskoy))@pval,
    testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = schema))@pval,
    testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = cross_disorder$Gene.updated.to.2025))@pval,
    testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = bd[[2]]))@pval
    #testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = hic))@pval
    #testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = asd))@pval
)

odds_ratios <- c(
    testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = mdd$Gene.Name))@odds.ratio,
    testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = adhd))@odds.ratio,
    testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = asd))@odds.ratio,
    testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = Geneset_with_disease_riskgenes$ASD_Ruzzo$`HGNC gene symbol`))@odds.ratio,
    testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = Geneset_with_disease_riskgenes$ASD_Sanders$RefSeqGeneName))@odds.ratio,
    testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = Geneset_with_disease_riskgenes$ASD_Satterstrom$Gene))@odds.ratio,
    #testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = Geneset_with_disease_riskgenes$DDD$symbol))@odds.ratio,
    testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = Trubetskoy_scz_genes$prioritized_genes_trubetskoy))@odds.ratio,
    testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = Trubetskoy_scz_genes$all_scz_trubetskoy))@odds.ratio,
    testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = schema))@odds.ratio,
    testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = cross_disorder$Gene.updated.to.2025))@odds.ratio,
    testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = bd[[2]]))@odds.ratio
    #testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = hic))@odds.ratio
    #testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = asd))@odds.ratio
)

# Adjust p-values using the Benjamini-Hochberg method
adjusted_pvalues <- p.adjust(pvalues, method = "BH")

# Create a data frame for plotting
test_names <- c("MDD", "ADHD", "ASD Cirnigliaro", "ASD Ruzzo", "ASD Sanders", "ASD Satterstrom", "SCZ Prioritized (Trub)", "SCZ All (Trub)", "SCZ-rare Schema", "Cross Disorder", "BD")
pvalue_df <- data.frame(Test = test_names, AdjustedPValue = adjusted_pvalues, OddsRatio = odds_ratios)
#pvalue_df <- data.frame(Test = test_names, PValue = pvalues, OddsRatio = odds_ratios)


# Get the intersecting genes for each test
intersecting_genes <- list(
    MDD = getIntersection(testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = mdd$Gene.Name))),
    ADHD = getIntersection(testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = adhd))),
    `ASD Cirnigliaro` = getIntersection(testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = asd))),
    `ASD Ruzzo` = getIntersection(testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = Geneset_with_disease_riskgenes$ASD_Ruzzo$`HGNC gene symbol`))),
    `ASD Sanders` = getIntersection(testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = Geneset_with_disease_riskgenes$ASD_Sanders$RefSeqGeneName))),
    `ASD Satterstrom` = getIntersection(testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = Geneset_with_disease_riskgenes$ASD_Satterstrom$Gene))),
    `SCZ Prioritized (Trub)` = getIntersection(testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = Trubetskoy_scz_genes$prioritized_genes_trubetskoy))),
    `SCZ All (Trub)` = getIntersection(testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = Trubetskoy_scz_genes$all_scz_trubetskoy))),
    `SCZ-rare Schema` = getIntersection(testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = schema))),
    `Cross Disorder` = getIntersection(testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = cross_disorder$Gene.updated.to.2025))),
    BD = getIntersection(testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = bd[[2]])))
    #`SCZ HiC` = getIntersection(testGeneOverlap(newGeneOverlap(genome.size = 29458, listA = unique(results_df_sig$feature), listB = hic)))
)

# Add intersecting genes to the data frame
pvalue_df$IntersectingGenes <- sapply(intersecting_genes, function(x) paste(x, collapse = ", "))

# Plot the -log10(adjusted p-values) with intersecting genes
ggplot(pvalue_df, aes(x = Test, y = -log10(AdjustedPValue), fill = OddsRatio)) +
    geom_bar(stat = "identity") +
    scale_fill_viridis() +
    xlab(label = "-log(p-val)") +
    theme_classic() +
    geom_hline(yintercept = -log10(0.05)) + coord_flip() +
    ggtitle(label = "Disorder risk enrichment in DEGs (p<0.01)") +
    geom_text(aes(label = IntersectingGenes), hjust = 1, size = 3, vjust = 0.5)

ggsave("./Results/Sertraline_DEG_disorder_enrichment_with_genes.pdf", width = 10, height = 5, dpi = 300)

