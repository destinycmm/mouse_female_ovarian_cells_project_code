##################################################################
### 05_Saturation analysis                                     ###
### Author: Manman cui                                         ###
### Date: 2021-12-28                                           ###
##################################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(export)

# Count_table in this dataset and GSE136441 dataset
this_dataset <- readRDS("/DataC/mmdata/mm/mouse_female_ovarian_cells.rds")
gse136441_dataset <- readRDS("/DataC/mmdata/mm/GSE136441_mouse_female_ovarian_cells.rds")
compared_germ_cells <- c("Mitotic", "Pre-meiotic", "Leptotene", "Zygotene", "E-Pachytene", "L-Pachytene", "Diplotene", "Pre-oocyte")
this_dataset_germ_cells <- subset(this_dataset,subset=celltype %in% compared_germcell)
this_dataset_germ_cells_count_table <- data.frame(this_dataset_germ_cells@assays$RNA@counts)
gse136441_dataset_germ_cells <- subset(gse136441_dataset,subset=celltype %in% compared_germcell)
gse136441_dataset_germ_cells_count_table <- data.frame(gse136441_dataset_germ_cells@assays$RNA@counts)
# Number of captured genes in 100 germ cells in this dataset
out_dir <- "/DataC/mmdata/mm/output/Saturation analysis/this_dataset/" 
for (sample in colnames(this_dataset_germ_cells_count_table)) {
  exp_gene = rownames(dplyr::select(this_dataset_germ_cells_count_table,all_of(sample)) %>% dplyr::filter(!!as.name(sample)> 0))
  out_file = paste0(out_dir,"/",sample,"_exp_gene.txt")
  write.table(exp_gene,out_file,quote = F,row.names = F,col.names = F)
}
this_dataset_germ_cells_file <- list.files(out_dir)
this_dataset_germ_cells_file_100 <- sample(this_dataset_germ_cells_file,100)
this_dataset_germ_cells_100_gene_list <- c()
this_dataset_germ_cells_100_sat_table <- data.frame(matrix(ncol = 2,nrow = length(this_dataset_germ_cells_file_100)))
colnames(this_dataset_germ_cells_100_sat_table) <- c("cell","gene")
for (i in 1:length(this_dataset_germ_cells_file_100)) {
  file_i <- paste0(out_dir,"/",this_dataset_germ_cells_file_100[i])
  table_i <- read.table(file_i,header = F)
  gene_i <- table_i$V1
  this_dataset_germ_cells_100_gene_list <- c(this_dataset_germ_cells_100_gene_list,gene_i)
  this_dataset_germ_cells_100_gene_list <- c(unique(this_dataset_germ_cells_100_gene_list))
  this_dataset_germ_cells_100_sat_table[i,1]=i
  this_dataset_germ_cells_100_sat_table[i,2]=length(this_dataset_germ_cells_100_gene_list)
}
this_dataset_germ_cells_100_sat_table["sample"]="this_dataset"
write.csv(this_dataset_germ_cells_100_sat_table,"/DataC/mmdata/mm/output/Saturation analysis/this_dataset/this_dataset_germ_cells_100_sat_table.csv")
# Number of captured genes in 600 germ cells in gse136441 dataset
out_dir <- "/DataC/mmdata/mm/output/Saturation analysis/gse136441_dataset/" 
for (sample in colnames(gse136441_dataset_germ_cells_count_table)) {
  exp_gene = rownames(dplyr::select(gse136441_dataset_germ_cells_count_table,all_of(sample)) %>% dplyr::filter(!!as.name(sample)> 0))
  out_file = paste0(out_dir,"/",sample,"_exp_gene.txt")
  write.table(exp_gene,out_file,quote = F,row.names = F,col.names = F)
}
gse136441_dataset_germ_cells_file <- list.files(out_dir)
gse136441_dataset_germ_cells_file_600 <- sample(gse136441_dataset_germ_cells_file,600)
gse136441_dataset_germ_cells_600_gene_list <- c()
gse136441_dataset_germ_cells_600_sat_table <- data.frame(matrix(ncol = 2,nrow = length(gse136441_dataset_germ_cells_file_600)))
colnames(gse136441_dataset_germ_cells_600_sat_table) <- c("cell","gene")
for (i in 1:length(gse136441_dataset_germ_cells_file_600)) {
  file_i <- paste0(out_dir,"/",gse136441_dataset_germ_cells_file_600[i])
  table_i <- read.table(file_i,header = F)
  gene_i <- table_i$V1
  gse136441_dataset_germ_cells_600_gene_list <- c(gse136441_dataset_germ_cells_600_gene_list,gene_i)
  gse136441_dataset_germ_cells_600_gene_list <- c(unique(gse136441_dataset_germ_cells_600_gene_list))
  gse136441_dataset_germ_cells_600_sat_table[i,1]=i
  gse136441_dataset_germ_cells_600_sat_table[i,2]=length(gse136441_dataset_germ_cells_600_gene_list)
}
gse136441_dataset_germ_cells_600_sat_table["sample"]="GSE136441_dataset"
write.csv(gse136441_dataset_germ_cells_600_sat_table,"/DataC/mmdata/mm/output/Saturation analysis/gse136441_dataset_germ_cells_600_sat_table.csv")
# mapping
merge_sat_table_100_600 <- rbind(this_dataset_germ_cells_100_sat_table,gse136441_dataset_germ_cells_600_sat_table)
merge_sat_table_100_600_plot <- ggplot(merge_sat_table_100_600,aes(cell,gene,color=sample))+geom_point(size=0.5)+geom_line()+
  theme_bw()+scale_color_manual(values = c("this_dataset" = "#ed918a","GSE136441_dataset" = "#ABB2B9"))
graph2pdf(merge_sat_table_100_600_plot,"/DataC/mmdata/mm/output/Saturation analysis/merge_sat_table_100_600_plot.pdf",width=7.99,height=6.66)



