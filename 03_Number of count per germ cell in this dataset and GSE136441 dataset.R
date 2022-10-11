###############################################################################
### 03_Number of count per germ cell in this dataset and GSE136441 dataset  ###
### Author: Manman cui                                                      ###
### Date: 2021-12-20                                                        ###
###############################################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(export)

# Meta.data in this dataset and GSE136441 dataset
this_dataset <- readRDS("/DataC/mmdata/mm/mouse_female_ovarian_cells.rds")
gse136441_dataset <- readRDS("/DataC/mmdata/mm/GSE136441_mouse_female_ovarian_cells.rds")
compared_germ_cells <- c("Mitotic", "Pre-meiotic", "Leptotene", "Zygotene", "E-Pachytene", "L-Pachytene", "Diplotene", "Pre-oocyte")

this_dataset_germ_cells <- subset(this_dataset,subset=celltype %in% compared_germ_cells)
this_dataset_germ_cells_meta <- this_dataset_germ_cells@meta.data
this_dataset_germ_cells_meta$platform <- "This dataset"
this_dataset_germ_cells_meta_2 <- dplyr::select(this_dataset_germ_cells_meta,celltype,nCount_RNA,nFeature_RNA,platform)
colnames(this_dataset_germ_cells_meta_2) <- c("celltype","nCount_RNA","nFeature_RNA","platform")

gse136441_dataset_germ_cells <- subset(gse136441_dataset,subset=celltype %in% compared_germ_cells)
gse136441_dataset_germ_cells_meta <- gse136441_dataset_germ_cells@meta.data
gse136441_dataset_germ_cells_meta$platform <- "GSE136441 dataset"
gse136441_dataset_germ_cells_meta_2 <- dplyr::select(gse136441_dataset_germ_cells_meta,celltype,nCount_RNA,nFeature_RNA,platform)
colnames(gse136441_dataset_germ_cells_meta_2) <- c("celltype","nCount_RNA","nFeature_RNA","platform")
merge_compared_germ_cells_meta <- rbind(this_dataset_germ_cells_meta_2,gse136441_dataset_germ_cells_meta_2)

# Mapping
germ_log10Count_RNA_boxplot <- ggplot(merge_compared_germ_cells_meta,aes(platform,log10(nCount_RNA)))+geom_violin(aes(color=platform))+
  geom_boxplot(aes(fill=platform),width=0.2,outlier.alpha = 0.1,outlier.size = 0.2)+
  geom_jitter(height = 0,size=0.5,alpha=0.1,width = 0.3)+
  theme_bw()+ggtitle("Germ cell")+
  theme(axis.title.x=element_blank(),panel.grid = element_blank(),plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values = c("This dataset" = "#ed918a", "GSE136441 dataset" = "#bfbfbf")) +
  scale_fill_manual(values =c("This dataset" = "#ed918a", "GSE136441 dataset" = "#bfbfbf"))
graph2pdf(germ_log10Count_RNA_boxplot,"/DataC/mmdata/mm/output/germ_log10Count_RNA_boxplot.pdf",width=6.99,height=6.66)

germ_nFeature_RNA_densityplot <- ggplot(merge_compared_germ_cells_meta,aes(nFeature_RNA,fill=platform))+
  geom_density(stat = "density",position = "identity")+ggtitle("Germ cell")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),panel.grid = element_blank(),legend.position = "none")+
  scale_color_manual(values = c("This dataset" = "#ed918a", "GSE136441 dataset" = "#bfbfbf")) +
  scale_fill_manual(values =c("This dataset" = "#ed918a", "GSE136441 dataset" = "#bfbfbf"))
mean(this_dataset_germ_cells_meta$nFeature_RNA)
mean(gse136441_dataset_germ_cells_meta$nFeature_RNA)
graph2pdf(germ_nFeature_RNA_densityplot,"/DataC/mmdata/mm/output/germ_nFeature_RNA_densityplot.pdf",width=6.99,height=6.66)

