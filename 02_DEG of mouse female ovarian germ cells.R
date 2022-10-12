##############################################################
### 02_DEG of mouse female ovarian germ cells              ###
### Author: Manman cui                                     ###
### Date: 2021-09-25                                       ###
##############################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(export)

# Find all DEG
mouse_female_ovarian_cells <- readRDS("/DataC/mmdata/mm/mouse_female_ovarian_cells.rds")
mouse_female_ovarian_germ_cells <- subset(mouse_female_ovarian_cells,
                                          subset=celltype %in% c("Mitotic", "Pre-meiotic", "Leptotene", "Zygotene", 
                                                                 "E-Pachytene", "L-Pachytene","Diplotene", "Pre-oocyte","Aberrant cell"))
allDEG_germ_cells <-  FindAllMarkers(mouse_female_ovarian_germ_cells,min.pct = 0.25,only.pos = T,test.use = "wilcox",logfc.threshold = 0.25)
allDEG_germ_cells_2 <- filter(allDEG_germ_cells,avg_log2FC>1,p_val_adj<0.01)
# DEG_pheatmap
germ_cell_meta_table <- mouse_female_ovarian_germ_cells@meta.data
germ_cell_meta_scaledata_table <-data.frame(mouse_female_ovarian_germ_cells@assays$RNA@scale.data)

c1_data <- germ_cell_meta_scaledata_table[rownames(germ_cell_meta_scaledata_table)%in% allDEG_germ_cells_2[allDEG_germ_cells_2$cluster=="Mitotic","gene"],
                                     filter(germ_cell_meta_table,celltype10=="Mitotic")$sample]
c1_data2 <- c1_data[order(rowMeans(c1_data),decreasing = T),order(colnames(c1_data),decreasing = T)]

c2_data <- germ_cell_meta_scaledata_table[rownames(germ_cell_meta_scaledata_table) %in% allDEG_germ_cells_2[allDEG_germ_cells_2$cluster=="Pre-meiotic","gene"],
                                     filter(germ_cell_meta_table,celltype10=="Pre-meiotic")$sample]
c2_data2 <- c2_data[order(rowMeans(c2_data),decreasing = T),order(colMeans(c2_data),decreasing = T)]

c3_data <- germ_cell_meta_scaledata_table[rownames(germ_cell_meta_scaledata_table) %in% allDEG_germ_cells_2[allDEG_germ_cells_2$cluster=="Leptotene","gene"],
                                     filter(germ_cell_meta_table,celltype10=="Leptotene")$sample]
c3_data2 <- c3_data[order(rowMeans(c3_data),decreasing = T),order(colMeans(c3_data),decreasing = T)]

c4_data <- germ_cell_meta_scaledata_table[rownames(germ_cell_meta_scaledata_table) %in% allDEG_germ_cells_2[allDEG_germ_cells_2$cluster=="Zygotene","gene"],
                                     filter(germ_cell_meta_table,celltype10=="Zygotene")$sample]
c4_data2 <- c4_data[order(rowMeans(c4_data),decreasing = T),order(colMeans(c4_data),decreasing = T)]

c5_data <- germ_cell_meta_scaledata_table[rownames(germ_cell_meta_scaledata_table)%in% allDEG_germ_cells_2[allDEG_germ_cells_2$cluster=="E-Pachytene","gene"],
                                     filter(germ_cell_meta_table,celltype10=="E-Pachytene")$sample]
c5_data2 <- c5_data[order(rowMeans(c5_data),decreasing = T),order(colMeans(c5_data),decreasing = T)]

c6_data <- germ_cell_meta_scaledata_table[rownames(germ_cell_meta_scaledata_table)%in% allDEG_germ_cells_2[allDEG_germ_cells_2$cluster=="L-Pachytene","gene"],
                                     filter(germ_cell_meta_table,celltype10=="L-Pachytene")$sample]
c6_data2 <- c6_data[order(rowMeans(c6_data),decreasing = T),order(colMeans(c6_data),decreasing = T)]

c7_data <- germ_cell_meta_scaledata_table[rownames(germ_cell_meta_scaledata_table)%in% allDEG_germ_cells_2[allDEG_germ_cells_2$cluster=="Diplotene","gene"],
                                     filter(germ_cell_meta_table,celltype10=="Diplotene")$sample]
c7_data2 <- c7_data[order(rowMeans(c7_data),decreasing = T),order(colMeans(c7_data),decreasing = T)]

c8_data <- germ_cell_meta_scaledata_table[rownames(germ_cell_meta_scaledata_table)%in% allDEG_germ_cells_2[allDEG_germ_cells_2$cluster=="Pre-oocyte","gene"],
                                     filter(germ_cell_meta_table,celltype10=="Pre-oocyte")$sample]
c8_data2 <- c8_data[order(rowMeans(c8_data),decreasing = T),order(colMeans(c8_data),decreasing = T)]

c9_data <- germ_cell_meta_scaledata_table[rownames(germ_cell_meta_scaledata_table)%in% allDEG_germ_cells_2[allDEG_germ_cells_2$cluster=="Aberrant cell","gene"],
                                     filter(germ_cell_meta_table,celltype10=="Aberrant cell")$sample]
c9_data2 <- c9_data[order(rowMeans(c9_data),decreasing = T),order(colMeans(c9_data),decreasing = T)]

germ_cell_meta_scaledata_table_2 <- germ_cell_meta_scaledata_table[c(rownames(c1_data2),rownames(c2_data2),rownames(c3_data2),rownames(c4_data2),rownames(c5_data2),
                                                          rownames(c6_data2),rownames(c7_data2),rownames(c8_data2),rownames(c9_data2)), 
                                                        c(colnames(c1_data2),colnames(c2_data2),colnames(c3_data2),colnames(c4_data2),colnames(c5_data2),
                                                          colnames(c6_data2),colnames(c7_data2),colnames(c8_data2),colnames(c9_data2))]
germ_cell_meta_scaledata_table_2[germ_cell_meta_scaledata_table_2>=2]=2
germ_cell_meta_scaledata_table_2[germ_cell_meta_scaledata_table_2<=-2]=-2
# annotation colname
annotation_col <- data.frame(
  cells = c(rep("Mitotic",339),rep("Pre-meiotic",157),rep("Leptotene",132),rep("Zygotene",135),rep("E-Pachytene",164),
            rep("L-Pachytene",200),rep("Diplotene",72),rep("Pre-oocyte",84),rep("Aberrant cell",54)))
rownames(annotation_col) <- colnames(germ_cell_meta_scaledata_table_2)
annotation_col$cells <- factor(annotation_col$cells, levels = c("Mitotic", "Pre-meiotic", "Leptotene", "Zygotene", 
                                                                "E-Pachytene", "L-Pachytene",
                                                                "Diplotene", "Pre-oocyte","Aberrant cell"))
# annotation rowname
annotation_row <- data.frame(
  genes = c(rep("Mitotic",913),rep("Pre-meiotic",109),rep("Leptotene",47),rep("Zygotene",194),rep("E-Pachytene",616),
            rep("L-Pachytene",496),rep("Diplotene",377),rep("Pre-oocyte",705),rep("Aberrant cell",813)))
rownames(annotation_row) <- rownames(germ_cell_meta_scaledata_table_2)
annotation_row$genes <- factor(annotation_row$genes, levels = c("Mitotic", "Pre-meiotic", "Leptotene", "Zygotene", 
                                                                "E-Pachytene", "L-Pachytene",
                                                                "Diplotene", "Pre-oocyte","Aberrant cell"))
# mapping
pheatmap(germ_cell_meta_scaledata_table_2,
         cluster_rows = F,cluster_cols = F,
         show_colnames = F,show_rownames = F,
         annotation_row = annotation_row,
         annotation_col = annotation_col)
graph2jpg(file="/DataC/mmdata/mm/output/C1_C9_DEG.jpg",width=8.51,height=7.72)

