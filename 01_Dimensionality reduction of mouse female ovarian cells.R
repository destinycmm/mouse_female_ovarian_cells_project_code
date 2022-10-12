###################################################################
### 01_Dimensionality reduction of mouse female ovarian cells   ###
### Author: Manman cui                                          ###
### Date: 2021-08-08                                            ###
###################################################################

library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)

# Read count_table @ meta_table
count_table <- read.csv("/DataC/mmdata/mm/mouse_female_ovarian_cells_count.csv",header = T,row.names = 1)
meta_table <- read.csv("/DataC/mmdata/mm/mouse_female_ovarian_cells_meta.csv",header = T,row.names = 1)
# Create object
mouse_female_ovarian_cells <- CreateSeuratObject(counts = count_table,meta.data = meta_table)
# QC
hemocyte_marker_genes <- c("Hba-a1", "Hba-a2", "Hba-x", "Hbb-y", "Hbb-bh1") 
hemocyte_marker_genes_m <- match(hemocyte_marker_genes,rownames(mouse_female_ovarian_cells@assays$RNA))
hemocyte_marker_genes <- rownames(mouse_female_ovarian_cells@assays$RNA)[hemocyte_marker_genes_m]
hemocyte_marker_genes <- hemocyte_marker_genes[!is.na(hemocyte_marker_genes)]
hemocyte_marker_genes_table <- data.frame(t(as.matrix(mouse_female_ovarian_cells@assays$RNA@counts[hemocyte_marker_genes,])))
hemocyte_marker_genes_table$sum_hemocyte_genes <- rowSums(hemocyte_marker_genes_table)
mouse_female_ovarian_cells@meta.data <- cbind(mouse_female_ovarian_cells@meta.data,hemocyte_marker_genes_table)
mouse_female_ovarian_cells <- subset(mouse_female_ovarian_cells, subset = nFeature_RNA > 2000 & nFeature_RNA < 15000 & sum_hemocyte_genes < 500)
# LogNormalize
mouse_female_ovarian_cells <- NormalizeData(mouse_female_ovarian_cells,normalization.method = "LogNormalize",scale.factor = 100000)
# Find highly variable genes
mouse_female_ovarian_cells <- FindVariableFeatures(mouse_female_ovarian_cells,selection.method = "vst",nfeatures = 2000)
# Scale data
all.genes <- rownames(mouse_female_ovarian_cells)
mouse_female_ovarian_cells <- ScaleData(mouse_female_ovarian_cells,features = all.genes)
# Dimensionality reduction
mouse_female_ovarian_cells <- RunPCA(mouse_female_ovarian_cells,features = VariableFeatures(object = mouse_female_ovarian_cells))
# Elbow plot analysis was used to ensure the most significant PC
ElbowPlot(mouse_female_ovarian_cells)
# clustering cells
mouse_female_ovarian_cells<- FindNeighbors(mouse_female_ovarian_cells,dims = 1:10)
mouse_female_ovarian_cells<- FindClusters(mouse_female_ovarian_cells,resolution = 1)
# Non-dimensionality reduction；
mouse_female_ovarian_cells <- RunTSNE(mouse_female_ovarian_cells,dims = 1:10)
DimPlot(mouse_female_ovarian_cells,reduction = "tsne",group.by = "time")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),panel.grid = element_blank())+
  scale_color_manual(values = c("E11.5"="#66C2A5","E12.5"="#ABDDA4","E13.5"="#999999",
                                "E14.5"="#999933","E15.5"="#CCCC99","E16.5"="#FFCC66","E17.5"="#FF9966","E18.5"="#CC6633"))
