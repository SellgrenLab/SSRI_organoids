## Assess proportion of neuronal subtypes very crudely ##

library(openxlsx)

df <- read.xlsx("201024_SSRI_pseudobulk_results_all.xlsx")

head(df)

genes_of_interest <- c("GRIA1", "GRIA2", "GRIA4", "DLG2", "NKAIN3", "ZEB2", "NRGN", "HSP90AA1", "HSP90AB1", "LGI2", "VAMP2", "CALM2", "CREB5", "SYT4")


df_filtered <- df[df$feature %in% genes_of_interest & 
                  df$celltype %in% c("Excitatory (newborn)", "Excitatory (Mature)", "Inhibitory"), -12 ]
head(df_filtered)

dim(df_filtered)


library(ggplot2)
df_filtered$signif <- ifelse(df_filtered$padj < 0.05, "*", "")

# Add asterisk labels for significant points
library(ggrepel)
ggplot(df_filtered, aes(x = feature, y = celltype, color = log_fc, size = abs(log_fc))) +
    geom_point() +
    # Add ring around significant points
    geom_point(
        data = subset(df_filtered, padj < 0.05),
        aes(x = feature, y = celltype, size = abs(log_fc)),
        shape = 21, fill = NA, color = "black", stroke = 1.0, inherit.aes = FALSE
    ) +
    # Add box around the plot area
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)
    ) +
    scale_color_gradient2(low = "darkblue", mid = "white", high = "darkred", midpoint = 0) +
    labs(x = "Gene", y = "Cell Type", color = "logFC", size = "|logFC|")

ggsave("genes_of_interest_dotplot.pdf", width = 11, height = 3.5, units = "in", device = "pdf")



gene_lists <- c(
    "DLX6-AS1,LHX6,ST8SIA5,PLS3,DLX5,THRB,BRINP2,NXPH1,DLX2,PDZRN3",
    "MYO5B,PPP1R1B,SLC24A2,NPR3,OPCML,NEUROD6,NEUROD2,SLC17A7,KIAA0319,RGS6",
    "MYO5B,PPP1R1B,SLC24A2,NPR3,OPCML,NEUROD6,NEUROD2,SLC17A7,KIAA0319,RGS6",
    "MYO5B,PPP1R1B,SLC24A2,NPR3,OPCML,NEUROD6,NEUROD2,SLC17A7,KIAA0319,RGS6",
    "MYO5B,PPP1R1B,SLC24A2,NPR3,OPCML,NEUROD6,NEUROD2,SLC17A7,KIAA0319,RGS6",
    "MYO5B,PPP1R1B,SLC24A2,NPR3,OPCML,NEUROD6,NEUROD2,SLC17A7,KIAA0319,RGS6",
    "MYO5B,PPP1R1B,SLC24A2,NPR3,OPCML,NEUROD6,NEUROD2,SLC17A7,KIAA0319,RGS6",
    "MYO5B,PPP1R1B,SLC24A2,NPR3,OPCML,NEUROD6,NEUROD2,SLC17A7,KIAA0319,RGS6",
    "MYO5B,PPP1R1B,SLC24A2,NPR3,OPCML,NEUROD6,NEUROD2,SLC17A7,KIAA0319,RGS6",
    "AL512590.3,SLA,SATB2,NEUROD2,SNCB,NEUROD6,EML6,AC068931.1,DPY19L1,TG",
    # Added gene lists below
    "PTH2,GNG8,HPAT5,GBX2,RNF220,TCF7L2,RASGEF1B,LINC01776,LHX9,SPOCK1",
    "PTH2,GNG8,HPAT5,GBX2,RNF220,TCF7L2,RASGEF1B,LINC01776,LHX9,SPOCK1",
    "PTH2,GNG8,HPAT5,GBX2,RNF220,TCF7L2,RASGEF1B,LINC01776,LHX9,SPOCK1",
    "PTH2,GNG8,HPAT5,GBX2,RNF220,TCF7L2,RASGEF1B,LINC01776,LHX9,SPOCK1",
    "SLC6A5,LAMP5,PAX8,PAX2,SST,CRHBP,LBX1,HOXB-AS3,HOXD3,OTP",
    "MYO5B,PPP1R1B,SLC24A2,NPR3,OPCML,NEUROD6,NEUROD2,SLC17A7,KIAA0319,RGS6",
    "MYO5B,PPP1R1B,SLC24A2,NPR3,OPCML,NEUROD6,NEUROD2,SLC17A7,KIAA0319,RGS6",
    "MYO5B,PPP1R1B,SLC24A2,NPR3,OPCML,NEUROD6,NEUROD2,SLC17A7,KIAA0319,RGS6",
    "MYO5B,PPP1R1B,SLC24A2,NPR3,OPCML,NEUROD6,NEUROD2,SLC17A7,KIAA0319,RGS6",
    "MYO5B,PPP1R1B,SLC24A2,NPR3,OPCML,NEUROD6,NEUROD2,SLC17A7,KIAA0319,RGS6",
    "MYO5B,PPP1R1B,SLC24A2,NPR3,OPCML,NEUROD6,NEUROD2,SLC17A7,KIAA0319,RGS6",
    "MYO5B,PPP1R1B,SLC24A2,NPR3,OPCML,NEUROD6,NEUROD2,SLC17A7,KIAA0319,RGS6",
    "MYO5B,PPP1R1B,SLC24A2,NPR3,OPCML,NEUROD6,NEUROD2,SLC17A7,KIAA0319,RGS6"
)

all_genes <- unlist(strsplit(gene_lists, ","))
unique_genes <- unique(c(all_genes, c("SLC17A7", "SLC17A6", "CAMK2A", "SATB2")))
print(unique_genes)




gene_lists2 <- c(
    "DLX6-AS1,LHX6,ST8SIA5,PLS3,DLX5,THRB,BRINP2,NXPH1,DLX2,PDZRN3",
    "DLX6-AS1,LHX6,ST8SIA5,PLS3,DLX5,THRB,BRINP2,NXPH1,DLX2,PDZRN3",
    "DLX6-AS1,LHX6,ST8SIA5,PLS3,DLX5,THRB,BRINP2,NXPH1,DLX2,PDZRN3",
    "LHX8,DLX5,DLX6-AS1,DLX2,DLX1,DLX6,GBX1,LHX6,ARX,SP9",
    "DLX6-AS1,LHX6,ST8SIA5,PLS3,DLX5,THRB,BRINP2,NXPH1,DLX2,PDZRN3",
    "DLX6-AS1,LHX6,ST8SIA5,PLS3,DLX5,THRB,BRINP2,NXPH1,DLX2,PDZRN3",
    "DLX6-AS1,LHX6,ST8SIA5,PLS3,DLX5,THRB,BRINP2,NXPH1,DLX2,PDZRN3",
    "DLX6-AS1,LHX6,ST8SIA5,PLS3,DLX5,THRB,BRINP2,NXPH1,DLX2,PDZRN3",
    "GNRH1,AC022433.1,ISL1,DLX6-AS1,DLX6,GPR149,TAC1,LRMDA,LINC01717,SIX3",
    "GNRH1,AC022433.1,ISL1,DLX6-AS1,DLX6,GPR149,TAC1,LRMDA,LINC01717,SIX3",
    "GNRH1,AC022433.1,ISL1,DLX6-AS1,DLX6,GPR149,TAC1,LRMDA,LINC01717,SIX3",
    "GNRH1,AC022433.1,ISL1,DLX6-AS1,DLX6,GPR149,TAC1,LRMDA,LINC01717,SIX3",
    "GNRH1,AC022433.1,ISL1,DLX6-AS1,DLX6,GPR149,TAC1,LRMDA,LINC01717,SIX3",
    "GNRH1,AC022433.1,ISL1,DLX6-AS1,DLX6,GPR149,TAC1,LRMDA,LINC01717,SIX3",
    "GNRH1,AC022433.1,ISL1,DLX6-AS1,DLX6,GPR149,TAC1,LRMDA,LINC01717,SIX3",
    "GNRH1,AC022433.1,ISL1,DLX6-AS1,DLX6,GPR149,TAC1,LRMDA,LINC01717,SIX3",
    "GNRH1,AC022433.1,ISL1,DLX6-AS1,DLX6,GPR149,TAC1,LRMDA,LINC01717,SIX3",
    "LHX8,DLX5,DLX6-AS1,DLX2,DLX1,DLX6,GBX1,LHX6,ARX,SP9",
    "DLX6-AS1,LHX6,ST8SIA5,PLS3,DLX5,THRB,BRINP2,NXPH1,DLX2,PDZRN3",
    "LHX8,DLX5,DLX6-AS1,DLX2,DLX1,DLX6,GBX1,LHX6,ARX,SP9"
)

all_genes2 <- unlist(strsplit(gene_lists2, ","))
unique_genes2 <- unique(c(all_genes2, c("GAD1", "GAD2", "SLC32A1", "PVALB", "SST", "VIP")))
print(unique_genes2)

library(Seurat)
ssri <- readRDS("201024_SSRI_final_clusters.rds")

ssri <- subset(ssri, subset = types %in% c("Excitatory (newborn)", "Excitatory (Mature)", "Inhibitory"))

ssri <- AddModuleScore(ssri, features = list(unique_genes), name = "Excitatory", assay = "SCT")
ssri <- AddModuleScore(ssri, features = list(unique_genes2), name = "Inhibitory", assay = "SCT")

exc_score_col <- "Excitatory1"
inh_score_col <- "Inhibitory1"

# Set delta threshold for classification
delta_threshold <- 0.01
ssri$neuron_type <- ifelse(
  ssri[[exc_score_col]] - ssri[[inh_score_col]] > delta_threshold,
  "Excitatory",
  ifelse(ssri[[inh_score_col]] - ssri[[exc_score_col]] > delta_threshold,
         "Inhibitory", "Unclassified")
)

# View proportion of each type
proportions <- prop.table(table(ssri$neuron_type))
print(proportions)

FeaturePlot(ssri, features = c(exc_score_col, inh_score_col))
