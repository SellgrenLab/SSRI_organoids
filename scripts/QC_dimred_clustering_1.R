# SSRI Project

## Standard data analysis for the SSRI scRNAseq (10X genomics) assay using brain organoids

## Load libraries

library(Seurat)
library(dplyr)
library(magrittr)
library(scCustomize)
library(presto)
library(DElegate)
library(openxlsx)
library(fgsea)
library(msigdbr)
library(ggplot2)

setwd("/Volumes/projects/C3_Sellgren_lab/Lab Members/Susmita/Internal data/CMV & SSRI/scRNAseq/SSRI/SSRI")

## Read in cleaned count  matrices from CellBender output
ctrl<- scCustomize::Read_CellBender_h5_Mat(file_name = "./sertraline/Control_ssri_bender_output_filtered.h5")
sert<- scCustomize::Read_CellBender_h5_Mat(file_name = "./sertraline/sertraline_bender_output_filtered.h5")

## Read in demultiplexed cell assignments and remove multiplets

demux_ctrl <- read.table("./sertraline/SSRI_control_clusters.tsv", header = TRUE, sep = "\t")
demux_sert <- read.table("./sertraline/SSRI_sert_clusters.tsv", header = TRUE, sep = "\t")

## Cell barcodes are as follows: ctrl-1 & sertraline-2

colnames(sert)<-gsub(pattern = 1,replacement = 2,colnames(sert))
demux_sert$barcode <-gsub(pattern = 1,replacement = 2,demux_sert$barcode)

###merge cell matrix of 2 samples

cm<-cbind(ctrl,sert)

###Create metadata
md<-rep(c('ctrl','sert'),times=c(length(colnames(ctrl)),
                                               length(colnames(sert))))
                                               
md<-as.data.frame(md)

row.names(md)<-c(colnames(ctrl),
                 colnames(sert))

colnames(md)<-'Condition'

#md$assignment <- ifelse(rownames(md)%in% demux_ctrl$barcode, demux_ctrl$assignment[match(rownames(md), demux_ctrl$barcode)], demux_sert$assignment[match(rownames(md), demux_sert$barcode)])
#md$multiplet_status <- ifelse(rownames(md)%in% demux_ctrl$barcode, demux_ctrl$status[match(rownames(md), demux_ctrl$barcode)], demux_sert$status[match(rownames(md), demux_sert$barcode)])

##Create Seurat object
project <- CreateSeuratObject(counts = cm, project = "SSRI",meta.data = md, min.cells =3, min.features = 200)
project <- scCustomize::Add_Cell_QC_Metrics(project, species = "human")



## Plot QC metrics
QC_Plot_UMIvsGene(seurat_object = project, low_cutoff_gene = 200, high_cutoff_gene = 12000, low_cutoff_UMI = 500, high_cutoff_UMI = 100000)

median_stats <- Median_Stats(seurat_object = project, group_by_var = "Condition")
# A tibble: 3 × 7
#  Condition          Median_nCount_RNA Median_nFeature_RNA Median_percent_mito Median_percent_ribo Median_percent_mito_ribo Median_log10GenesPerUMI
#  <chr>                          <dbl>               <dbl>               <dbl>               <dbl>                    <dbl>                   <dbl>
#1 ctrl                          11812.                4199                5.06                11.5                     17.5                   0.870
#2 sert                           4114                 1669                5.96                13.7                     23.3                   0.870
#3 Totals (All Cells)             7024                 2605                5.40                12.5                     19.6                   0.870


## Set filter threshold
nFeature_lowlimit=200
nFeature_highlimit=12000
nCount_RNA_highlimit=100000
mito_cutoff=15
#ribo_cutoff=40
pcs_number=30
resolution=0.6
top_num=10

##filter cells based on QC thresholds (Filtering Round II) the first filtering was done during ambient RNA removal
project <- subset(project, subset = nFeature_RNA > nFeature_lowlimit  &   nFeature_RNA < nFeature_highlimit &
                      nCount_RNA <  nCount_RNA_highlimit & 
                      percent_mito < mito_cutoff)


## Get genotype demultiplexed data and add to the object
demux_sert$Condition <- "Sertraline"
demux_ctrl$Condition <- "Control"
demux<- rbind(demux_ctrl, demux_sert)
project$assignment <- demux$assignment[match(rownames(project@meta.data), demux$barcode)]
project$multiplet_status <- demux$status[match(rownames(project@meta.data), demux$barcode)]
project$doublet_prob <- demux$log_prob_doublet[match(rownames(project@meta.data), demux$barcode)]

## We haven't filtered out doublets yet, lets do the first round of clustering and see later if the assigned doublets stand out in clustering
table(project$Condition, project$multiplet_status)
      
#       doublet singlet unassigned
#  ctrl     357    7323         11
#  sert     258    5425         11


## Normalize data using SCTransform
#options(future.globals.maxSize= 10091289600)
project <- SCTransform(project, vars.to.regress = c("percent_mito"))
project <- RunPCA(project)

## Integrate the two batches for cell type identification
project<-RunHarmony(project,"Condition", plot_convergence = TRUE,dims.use = 1:30)

## Embed and visualize the clustering
project <- RunUMAP(object = project, reduction = "harmony",dims=1:20)

UMAPPlot(project, group.by = c("Condition", "assignment"),label=TRUE) # check if integration was efficient across condition, if not use another methods to integrate. Also check if individual-specific effects are present based on demultiplexed genotypes

## The cluster with NA assignment status does not show any relevant markers and has a high ribosomal percentage of transcripts. This cluster, while valid, probably won't be useful in the project so we only keep the assigned singlet cells, i.e., high quality cells
project1 <- subset(project, subset=multiplet_status %in% c("singlet"))

project1 <- RunUMAP(project1, reduction = "harmony",dims=1:20, min.dist = 0.2)

## Once the embedding is finalized and cluster separation does not look confounded by QC metrics, proceed with calling clusters.
project1 <- FindNeighbors(project1, reduction = "harmony", dims = 1:20)
project1 <- FindClusters(resolution = c(0.4,0.6,0.8), project1)

UMAPPlot(project1, label=T, group.by= 'SCT_snn_res.0.4')

## Get markers for clusters

#res<- wilcoxauc(project, 'SCT_snn_res.0.4')
#res <- dplyr::filter(res, padj<0.05)

markers <- DElegate::findDE(project1, group_column='SCT_snn_res.0.4', method = 'deseq', replicate_column = 'assignment')
View(markers %>%
         group_by(group1) %>%
         dplyr::filter(padj < 0.05 & log_fc > 0)  %>%
         slice_head(n = 10) %>%
         ungroup())

top10 <-  markers %>%
         group_by(group1) %>%
         dplyr::filter(padj < 0.05 & log_fc > 0)  %>%
         slice_head(n = 10) %>%
         ungroup()

early <- c("SOX2", "PAX6", "NEUROG2", "NEUROD1", "GFAP", "FABP7", "PTPRZ1", "EMX2", "HES1", "HOPX", "HES6", "PROM1", "LIFR", "ASCL1", "DCX", "GAD1", "SLC17A6","GLI3", "EGFR", "OLIG2")
neuron <- c("MAP2", "DCX", "STMN2","SLC17A7", "GAD1", "GAD2", "PROX1", "GRIK4","SLC17A6") # mostly hippocampal and thalamic neurons

## Plot results and save
FeaturePlot(project1, features=c("nCount_RNA", "nFeature_RNA", "percent_mito", "percent_oxphos", "percent_apop", "percent_dna_repair","percent_ieg", "percent_top50", "G2M.Score"), cols=c("grey100", "firebrick"))
DimPlot_scCustom(project1,group.by=c("Condition", "assignment", "SCT_snn_res.0.4", "SCT_snn_res.0.6"), label = T)

DimPlot_scCustom(project1, split_seurat = T, split.by = "Condition", figure_plot = T, colors_use = DiscretePalette_scCustomize(num_colors = 10, palette = "alphabet2"))
## Save the above plots in the QC folder

write.xlsx(top10, "./sertraline/QC/Top10_cluster_markers_res4.xlsx")

saveRDS(project1, "./sertraline/181024_SSRI_processed_object.rds")




####### Cell type Annotation ########

## Load the processed object
project <- readRDS("./sertraline/181024_SSRI_processed_object.rds")
levels(project)

project <- RenameIdents(project, '0'= "apical RG", '1'="Astroglia", '2'="oRG", '3'= "RG", '4'="Neurons-1",'5'="undefined", '6'="Glioblast", '7'="Cycling RG", '8'="Neurons-2",'9'="Neuroepithelial")

project$celltype <- Idents(project)
DimPlot_scCustom(project1, label = T)

### Check QC for clustering
project@meta.data %>%
ggplot(aes(x=Condition, fill=Condition)) +
geom_bar() +
theme_classic() + viridis::scale_fill_viridis(discrete=T)+
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
theme(plot.title = element_text(hjust=0.5, face="bold")) +
ggtitle("NCells")

QC_Plots_UMIs(seurat_object = project, low_cutoff = 2000, high_cutoff = 100000, plot_median = T)
QC_Plots_Genes(seurat_object = project, low_cutoff = 2000, high_cutoff = 12000, plot_median = T)
QC_Plots_Complexity(seurat_object = project, high_cutoff = 0.8, plot_median = T)

project@meta.data %>%
ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent_mito)) +
geom_point() +
scale_colour_gradient(low = "gray90", high = "forestgreen") +
stat_smooth(method=lm) +
scale_x_log10() +
scale_y_log10() +
theme_classic() +
geom_vline(xintercept = 1000) +
geom_hline(yintercept = 250) +
facet_wrap(~celltype)

####### Stricter thresholding for cell type identification ########

project1 <- subset(project, subset = nFeature_RNA > 1000  &
nCount_RNA > 2000 &
percent_mito < 10)
project1 <- SCTransform(project1, vars.to.regress = c("percent_mito"))
project1 <- RunPCA(project1)
project1 <- RunUMAP(object = project1,dims=1:30)
project1<-RunHarmony(project1,"Condition", plot_convergence = TRUE,dims.use = 1:30)
project1 <- RunUMAP(object = project1, reduction = "harmony",dims=1:30)
UMAPPlot(project1, group.by="celltype", cols= colors_dutch)
project1 <- FindNeighbors(project1, reduction = "harmony", dims = 1:30)
project1 <- FindClusters(resolution = c(0.5), project1)



####### Neuron subclustering ########


neurons <- subset(project, subset=celltype%in%c("Neurons-1", "Neurons-2"))
neurons <- SCTransform(neurons, vars.to.regress= c("percent_mito"))
neurons <- RunPCA(neurons)
neurons <- RunUMAP(neurons, dims=1:30)
UMAPPlot(neurons,group.by="Condition")
neurons <- FindNeighbors(neurons, dims = 1:30)
neurons <- FindClusters(resolution = c(0.2), neurons)
neuron_markers <- DElegate::findDE(neurons, method = 'deseq', replicate_column = 'assignment')
FeaturePlot_scCustom(neurons, reduction = "umap", features = c("MAP2", "DCX", "STMN2", "GAD1", "GAD2", "PROX1","CALB1", "GRIK4","SLC17A6"), figure_plot = T)
neurons <- RenameIdents(neurons, '0'='Excitatory (Mature)', '1'= 'Excitatory (newborn)', '2'='Inhibitory', '3' = 'stressed', '4' = 'oRG-derived')
neurons$subtypes <- Idents(neurons)
DimPlot_scCustom(neurons, figure_plot=T, label.box = T, colors_use= DiscretePalette_scCustomize(num_colors = 5, palette = "alphabet2"))
FeaturePlot(neurons, features=c( "percent_oxphos"), cols=c("grey100", "firebrick"))
saveRDS(neurons,"./sertraline/201024_SSRI_neurons_subclusters.rds")



####### Final Annotations for downstream processing ########


project1$types <- project1$celltype
project1$types <- as.character(project1$types)
project1$types <- ifelse(project1$SCT_snn_res.0.5=='10', "IPC", project1$types) # had clearer expression of IPC markers

neurons$subtypes <- as.character(neurons$subtypes)
project1$types <- ifelse(rownames(project1@meta.data)%in%rownames(neurons@meta.data), neurons$subtypes[match(rownames(project1@meta.data), rownames(neurons@meta.data))], project1$types)

DimPlot_scCustom(project1, figure_plot=T, label.box = T, colors_use= DiscretePalette_scCustomize(num_colors = 20, palette = "alphabet2"), group.by="types")


saveRDS(project1, "./sertraline/201024_SSRI_final_annotations.rds")
