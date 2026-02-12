############################################################
#   Pseudotime trajectory analysis comparing Control vs
#   Sertraline using slingshot, condiments, and tradeSeq
############################################################
suppressPackageStartupMessages({
  library(slingshot)
  library(SingleCellExperiment)
  library(RColorBrewer)
  library(scales)
  library(viridis)
  library(pheatmap)
  library(msigdbr)
  library(fgsea)
  library(knitr)
  library(ggplot2)
  library(gridExtra)
  library(tradeSeq)
  library(condiments)
  library(dplyr)
  library(tidyr)
  library(Matrix)
  library(Seurat)
  library(BiocParallel)
})

set.seed(3)

# Subset Seurat Object

seo <- subset(
  project1,
  subset = types %in% c(
    "apical RG", "Astroglia", "Cycling RG",
    "Excitatory (Mature)", "Excitatory (newborn)",
    "Glioblast", "Inhibitory", "IPC",
    "RG", "Neuroepithelial", "oRG"
  )
)

sce <- as.SingleCellExperiment(seo, assay = "SCT")

############################
# Imbalance Score
############################

shuffle <- sample(ncol(sce))

scores <- imbalance_score(
  rd = reducedDims(sce)$UMAP,
  cl = colData(sce)$Condition,
  k = 20,
  smooth = 40
)

grad <- viridis::plasma(10)
names(grad) <- levels(cut(scores$scaled_scores, breaks = 10))

plot(
  reducedDims(sce)$UMAP,
  col = grad[cut(scores$scaled_scores, breaks = 10)],
  asp = 1, pch = 16,
  xlab = "UMAP-1", ylab = "UMAP-2", cex = 0.8
)
legend("topleft", legend = names(grad), col = grad,
       pch = 16, bty = "n", cex = 0.7)

############################
# 4. Run Slingshot
############################

sce <- slingshot(
  sce,
  reducedDim = "UMAP",
  clusterLabels = colData(sce)$types,
  start.clus = "Neuroepithelial",
  approx_points = 150
)

############################
# Plot Lineages
############################

plot_lineage <- function(lineage_id) {
  pt <- sce[[paste0("slingPseudotime_", lineage_id)]]
  plot(
    reducedDims(sce)$UMAP[shuffle, ],
    col = hcl.colors(100, alpha = 0.5)[cut(pt, breaks = 100)][shuffle],
    asp = 1, pch = 16,
    xlab = "UMAP-1", ylab = "UMAP-2"
  )
  lines(SlingshotDataSet(sce))
}

# Lineage 1 (Inhibitory)
plot_lineage(1)

# Lineage 4 (Excitatory Mature)
plot_lineage(4)

############################
# Density Comparison
############################

plot_density_comparison <- function(lineage_name) {

  x <- na.omit(
    slingPseudotime(sce)[colData(sce)$Condition == "ctrl", lineage_name]
  )
  y <- na.omit(
    slingPseudotime(sce)[colData(sce)$Condition == "sert", lineage_name]
  )

  ds <- list(Control = density(x),
             Sertraline = density(y))

  xlim <- range(c(ds$Control$x, ds$Sertraline$x))
  ylim <- range(c(ds$Control$y, ds$Sertraline$y))

  plot(xlim, ylim, type = "n",
       xlab = "Pseudotime", ylab = "")

  polygon(
    c(min(ds$Control$x), ds$Control$x, max(ds$Control$x)),
    c(0, ds$Control$y, 0),
    col = alpha(brewer.pal(4, "Set1")[3], 0.5)
  )

  polygon(
    c(min(ds$Sertraline$x), ds$Sertraline$x, max(ds$Sertraline$x)),
    c(0, ds$Sertraline$y, 0),
    col = alpha(brewer.pal(4, "Set1")[2], 0.5)
  )

  legend("topleft",
         legend = c("Control", "Sertraline"),
         fill = alpha(brewer.pal(4, "Set1")[3:2], 0.5),
         bty = "n")
}

plot_density_comparison("Lineage1")
plot_density_comparison("Lineage4")

############################
# Progression & Differentiation Tests
############################

sds <- SlingshotDataSet(sce)

prog_res <- progressionTest(
  sds,
  conditions = colData(sce)$Condition,
  global = TRUE,
  lineages = TRUE
)

print(knitr::kable(prog_res))

set.seed(2695)

dif_res <- differentiationTest(
  sds,
  conditions = colData(sce)$Condition,
  global = FALSE,
  pairwise = TRUE
)

print(knitr::kable(dif_res))

############################
# tradeSeq Modeling
############################

BPPARAM <- bpparam()
BPPARAM$progressbar <- TRUE
BPPARAM$workers <- parallel::detectCores()

df <- data.frame(row.names = colData(sce)$barcode)
df$Condition <- factor(colData(sce)$Condition)
df$sample <- factor(colData(sce)$assignment)

U <- model.matrix(~ Condition + sample, data = df)

sce_gam <- fitGAM(
  counts = counts(sce),
  sds = SlingshotDataSet(sce),
  conditions = factor(colData(sce)$Condition),
  nknots = 6,
  U = U,
  sce = TRUE,
  parallel = TRUE,
  BPPARAM = BPPARAM
)

saveRDS(sce_gam, "Sertraline_sce_tradeseq_fitted.rds")

############################################################
# End of Script
############################################################
