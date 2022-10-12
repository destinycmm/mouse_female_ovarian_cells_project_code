##################################################################
### 04_Estimation of dropout event                             ###
### Author: Manman cui                                         ###
### Date: 2021-12-27                                           ###
##################################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(export)

overlap_gene_list <-read.table("C:/Users/Oo/Desktop/雌性/202106-/0001-作图文件/6-作图代码和材料/Fig1_Fig2/Fig1H&FigS3D/overlap-genes-19477.txt",header = F)
# Read overlapped genes list
overlap_gene_list <- read.table("/DataC/mmdata/mm/output/overlap-genes-19477.txt",header = F)
names(overlap_gene_list)="gene"
# Creat dropout list
overlap_gene_dropout_ratio_table <- data.frame(matrix(ncol = 3,nrow = length(overlap_gene_list$gene)))
colnames(overlap_gene_dropout_ratio_table)=c("gene","this_dataset","gse136441_dataset")
rownames(overlap_gene_dropout_ratio_table)=overlap_gene_list$gene
overlap_gene_dropout_ratio_table$gene=rownames(overlap_gene_dropout_ratio_table)
# Data_table in this dataset and GSE136441 dataset
this_dataset <- readRDS("/DataC/mmdata/mm/mouse_female_ovarian_cells.rds")
gse136441_dataset <- readRDS("/DataC/mmdata/mm/GSE136441_mouse_female_ovarian_cells.rds")
compared_germ_cells <- c("Mitotic", "Pre-meiotic", "Leptotene", "Zygotene", "E-Pachytene", "L-Pachytene", "Diplotene", "Pre-oocyte")
this_dataset_germ_cells <- subset(this_dataset,subset=celltype %in% compared_germ_cells)
this_dataset_germ_cells_data_table <- data.frame(this_dataset_germ_cells@assays$RNA@data)
gse136441_dataset_germ_cells <- subset(gse136441_dataset,subset=celltype %in% compared_germ_cells)
gse136441_dataset_germ_cells_data_table <- data.frame(gse136441_dataset_germ_cells@assays$RNA@data)
# Dropout ratio of each gene in this dataset and GSE136441 dataset
for (gene in overlap_gene_dropout_ratio_table$gene) {
  this_dataset_dropout_ratio <- sum(this_dataset_germ_cells_data_table[gene,]==0)/ncol(this_dataset_germ_cells_data_table)
  gse136441_dropout_ratio <- sum(gse136441_dataset_germ_cells_data_table[gene,]==0)/ncol(gse136441_dataset_germ_cells_data_table)
  overlap_gene_dropout_ratio_table[gene,"this_dataset"]=this_dataset_dropout_ratio
  overlap_gene_dropout_ratio_table[gene,"gse136441_dataset"]=gse136441_dropout_ratio
}
write.table(overlap_gene_dropout_ratio_table,"/DataC/mmdata/mm/output/overlap_gene_dropout_ratio_table.txt",sep = "\t",quote = F,row.names = F)
# mapping
genes <- c("Gadl1","Lgr4","Kit","Slc25a1","Ddx4")
mark_gene <- dplyr::filter(overlap_gene_dropout_ratio_table,gene %in% genes)
overlap_gene_dropout_ratio_plot <- ggplot()+geom_point(overlap_gene_dropout_ratio_table,mapping = aes(this_dataset,gse136441_dataset),size=0.5,alpha=0.1)+
  geom_point(mark_gene,mapping = aes(this_dataset,gse136441_dataset),color="#ed918a")+
  geom_text(mark_gene,mapping = aes(this_dataset,gse136441_dataset,label=gene),hjust=-0.2)+theme_bw()+theme(panel.grid = element_blank())
graph2pdf(overlap_gene_dropout_ratio_plot,"/DataC/mmdata/mm-study/output/作图代码和材料/Fig1/Fig1J/overlap_gene_dropout_ratio_plot.pdf",width=4.52,height=3.51)


