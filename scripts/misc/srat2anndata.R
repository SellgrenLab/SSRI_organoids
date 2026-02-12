###### Convert seurat to anndata object for MELD analysis ########

project1 <- readRDS("./sertraline/201024_SSRI_final_annotations.rds")

library(SeuratDisk)

project2 <- project1
project2[["RNA3"]] <- as(object = project1[["RNA"]], Class = "Assay")
DefaultAssay(project2) <- "RNA3"
project2[["RNA"]] <- NULL
project2 <- RenameAssays(object = project2, RNA3 = 'RNA')

SaveH5Seurat(project2, filename = "SSRI.h5Seurat")
Convert("SSRI.h5Seurat", dest = "h5ad")
