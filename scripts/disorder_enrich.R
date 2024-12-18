library(GeneOverlap)
library(ggplot2)
library(dplyr)

Geneset_with_disease_riskgenes <- readRDS("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/Internal data/C1Q/data/Geneset_with_disease_riskgenes.rds")
Trubetskoy_scz_genes <- readRDS("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/Internal data/C1Q/Bipolar/Trubetskoy_scz_genes.rds")
asd<-  readRDS("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/Internal data/C1Q/Bipolar/ASD_riskgenes.rds") #https://doi-org.proxy.kib.ki.se/10.1073/pnas.2215632120
adhd<-  readRDS("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/Internal data/C1Q/Bipolar/ADHD_riskgenes.rds") #demontis et al
mdd<- read.xlsx("/Users/susmita.malwade/MDD.xlsx") #https://www-nature-com.proxy.kib.ki.se/articles/s41588-023-01596-4


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
results_de <- read.xlsx("./sertraline/201024_SSRI_pseudobulk_results_all.xlsx")
results_df_sig <- dplyr::filter(results_de, padj < 0.05)

df1<- data.frame("Disorder" = c("ASD(Cirnigliaro)", "ASD(Satterstom)", "ASD(Ruzzo)", "ASD (Sanders)", "DDD", "ADHD", "MDD"))
df1$pval <- c(1.1e-05,8.6e-05,7.4e-05,5.3e-05,8.3e-06,0.26,1.9e-03)
df1$OR <- c(10.4, 7.3,9.4, 10.1,4.8,3.5, 3.4)
df1$padj <- p.adjust(df1$pval, method = "BH")
df1$plot <- -log10(df1$padj)
## Plot the results

plot2<- ggplot(data= df1, mapping = aes(x=plot, y=Disorder, fill=OR)) + geom_bar(stat="identity") + scale_fill_viridis() + xlab(label = "-log(p-adjusted)") + theme_classic() + geom_vline(xintercept = -log10(0.05)) + ggtitle(label = "Disorder risk enrichment in DEGs")


cowplot::plot_grid(plot1,plot2,ncol = 1)
