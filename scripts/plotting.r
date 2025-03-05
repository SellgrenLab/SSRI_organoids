
## SSRI project - Pathway Heatmap ##
library(dplyr)
library(openxlsx)
library(ggplot2)
library(reshape2)
library(pals)

data <- read.xlsx("SSRI_pathways_selected.xlsx", sheet = 1)
#head(data)

# Reshape the data for heatmap 
data_melt <- melt(data, id.vars = c("Pathway", "Cell.type", "Group"), measure.vars = "NES")
data_cast <- dcast(data_melt, Pathway + Group +  Cell.type ~ variable, value.var = "value")
data_melt <- melt(data, id.vars = c("Pathway", "Group"), measure.vars = "NES")
data_cast <- dcast(data_melt, Pathway + Group ~ variable, value.var = "value")


# Create the heatmap
ggplot(data_cast, aes(x = Cell.type, y = Pathway, fill = NES)) +
    geom_tile() +
    facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient2(low = "#1c1c6d", high = "#b4050584", mid = "white", midpoint = 0) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = "Pathway Heatmap", x = "Cell Type", y = "Pathway", fill = "NES")

# Create the heatmap with viridis color scale
ggplot(data_cast, aes(x = Cell.type, y = Pathway, fill = NES)) +
        geom_tile() +
        facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
        scale_fill_viridis_c(option = "viridis") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(title = "Pathway Heatmap", x = "Cell Type", y = "Pathway", fill = "NES")

   
pal.bands(ocean.balance, ocean.delta, ocean.curl, main = "cmocean")   
# Create the heatmap with ocean.curl color scale
ggplot(data_cast, aes(x = Cell.type, y = Pathway, fill = NES)) +
    geom_tile() +
    facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradientn(colors = ocean.balance(100)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = "Pathway Heatmap", x = "Cell Type", y = "Pathway", fill = "NES")

    # Save the last plot as a PDF
ggsave("SSRI_selected_pathway_heatmap.pdf", width = 10, height = 12)



ssri <- readRDS("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/Internal data/CMV & SSRI/scRNAseq/SSRI/SSRI/sertraline/201024_SSRI_final_clusters.rds")

VlnPlot2(ssri, features=c("ARID1B", "MBD5", "TRIO", "RFX3"), cells=colnames(ssri)[ssri$types %in% c("RG")], stat.method = "wilcox.test", group.by="Condition") 
ggsave("SSRI_Disease-genes_RG_violins.pdf", width = 10, height = 12)

DotPlot2(ssri, features=c("ARID1B", "MBD5", "TRIO", "RFX3"), group.by="types", split.by="Condition", color_scheme = "BuRd")
ggsave("SSRI_Disease-genes_all_dot.pdf", width = 10, height = 5)


gene_groups <- list("ExN"= c("STMN2","MAP2","NEUROD6", "NRP1", "UNC5D", "KCNQ3", "ELAVL4", "FGF12", "SATB2", "MEF2C", "CRYM", "FEZF2", "NPR3"), "IN"= c( "ADARB2","GAD2", "DLX1", "DLX6-AS1"), "RG"= c( "VIM", "GLI3", "FABP7", "NES", "HOPX", "PAX6", "LIX1", "HMGA2"),"Astrocyte"=c("AQP4", "GFAP", "S100B"), "Glioblast"=c("EGFR", "DLK1", "ALDOA", "NTRK2", "PTN", "EDNRB", "SPARCL1") )

toplot2 <- CalcStats(ssri, features = unlist(gene_groups), method = "zscore", order = "p")
Heatmap(toplot2, lab_fill = "zscore") +
  theme(axis.text.y = element_blank())

ssri <- AddModuleScore(ssri, features = gene_groups[[1]], name = "ExN")
ssri <- AddModuleScore(ssri, features = gene_groups[[2]], name = "IN")
ssri <- AddModuleScore(ssri, features = gene_groups[[3]], name = "RG")
ssri <- AddModuleScore(ssri, features = gene_groups[[4]], name = "Ast")
ssri <- AddModuleScore(ssri, features = gene_groups[[5]], name = "GPC")
ssri <- AddModuleScore(ssri, features = c(gene_groups[[1]], gene_groups[[2]]), name = "Neu")

DimPlot2(
  ssri,
  features = "ExN1",
  cols = "Spectral-rev",
  theme = NoAxes())
ggsave("SSRI_ExN1_DimPlot.pdf", width = 5, height = 5)


DimPlot2(
  ssri,
  features = "IN1",
  cols = "Spectral-rev",
  theme = NoAxes())
ggsave("SSRI_IN_DimPlot.pdf", width = 5, height = 5)

DimPlot2(
  ssri,
  features = "RG1",
  cols = "Spectral-rev",
  theme = NoAxes())
ggsave("SSRI_RG_DimPlot.pdf", width = 5, height = 5)


DimPlot2(
  ssri,
  features = "Ast1",
  cols = "Spectral-rev",
  theme = NoAxes())
ggsave("SSRI_Ast_DimPlot.pdf", width = 5, height = 5)

DimPlot2(
  ssri,
  features = "GPC1",
  cols = "Spectral-rev",
  theme = NoAxes())
ggsave("SSRI_Gliob_DimPlot.pdf", width = 5, height = 5)


DimPlot2(
  ssri,
  features = "Neu1",
  cols = "Spectral-rev",
  theme = NoAxes())
ggsave("SSRI_Neurons_DimPlot.pdf", width = 5, height = 5)



### Other plots ####

DotPlot2(ssri, group = "Condition", color = "group", size = "n", alpha = 0.5, 
          xlab = "SSRI", ylab = "Group", main = "SSRI by Group")

syn<- c("GRIA2", "GRIA1", "GRIA4", "VAMP2", "CALM2", "DLG2")
DotPlot2(ssri, features=syn, group.by="types", split.by="Condition", color_scheme = "BuRd")


VlnPlot2(ssri , features = syn, group.by = "types", split.by = "Condition", cells=colnames(ssri)[ssri$types %in% c("Excitatory (Mature)", "Excitatory (newborn)", "Inhibitory")], stat.method = "wilcox.test", 
         hide.ns = TRUE)
ggsave("ssri_synaptic genes in neurons.pdf", width = 10, height = 10)


