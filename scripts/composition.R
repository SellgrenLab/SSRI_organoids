library(Seurat)
library(cacoa)
library(ggplot2)

sample.groups <- c(rep("Control", 4), rep("Sertraline", 4))
names(sample.groups) <- c("C0", "C1", "C2", "C3", "S0", "S1", "S2", "S3")
project$replicate <- NA
for(i in 1:length(project$Condition)){ project$replicate[i] <- ifelse(project$Condition[i]=="ctrl", paste0("C", project$assignment[i]), paste0("S", project$assignment[i])) }
cell.groups <- project$celltype
names(cell.groups) <- rownames(project@meta.data)
sample.per.cell <- project$replicate
names(sample.per.cell) <- rownames(project@meta.data)
ref.level <- "Control"
target.level <- "Sertraline"

cao$plot.params <- list(size=0.1, alpha=0.1, font.size=c(2, 3))
cao$plot.theme <- cao$plot.theme + theme(legend.background=element_blank())


cao <- Cacoa$new(project[["SCT"]]@data, sample.groups=sample.groups, cell.groups=cell.groups, sample.per.cell=sample.per.cell, 
                 ref.level=ref.level, target.level=target.level, embedding=project@reductions$umap@cell.embeddings)cao$estimateCellLoadings()

cao$estimateCellLoadings()
cao$estimateExpressionShiftMagnitudes()

cao$plotCellLoadings(show.pvals=FALSE)
cao$plotExpressionShiftMagnitudes()


######

library(miloR)
embryo_milo <- Milo(as.SingleCellExperiment(project))
embryo_milo <- buildGraph(embryo_milo, k = 40, d = 30) #thumb rule for k is no. of samples*5 = 8*5
embryo_milo <- makeNhoods(embryo_milo, prop = 0.1, k = 40, d=30, refined = TRUE)
embryo_milo <- countCells(embryo_milo, meta.data = data.frame(colData(embryo_milo)), sample="replicate")
plotNhoodSizeHist(embryo_milo)
embryo_milo <- calcNhoodDistance(embryo_milo, d=30)
da_results <- testNhoods(embryo_milo, design = ~ Condition, design.df = embryo_design)
embryo_milo <- buildNhoodGraph(embryo_milo)

plotUMAP(embryo_milo) + plotNhoodGraphDA(embryo_milo, da_results, alpha=0.05) +
  plot_layout(guides="collect")



