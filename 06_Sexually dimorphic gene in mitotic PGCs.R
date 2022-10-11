##################################################################
### 06_Sexually dimorphic gene in mitotic PGCs                 ###
### Author: Manman cui                                         ###
### Date: 2022-02-12                                           ###
##################################################################

library(Seurat)
library(dplyr)
library(openxlsx)

mitotic_PGCs <- readRDS("/DataC/mmdata/mm/XX_XY_mitotic_PGCs.rds")
# XX_undiff.ed_PGCs VS XY_undiff.ed_PGCs
undiffed_mitotic <- subset(mitotic_PGCs,subset= celltype=="Undiff.ed mitotic PGCs")
undiffed_mitotic <- SetIdent(mitotic_PGCs,value = undiffed_mitotic@meta.data$sex)
undiff_XX_VS_XY_DEG <- FindMarkers(undiffed_mitotic,ident.1 = "Female",ident.2 = "Male",
                               min.pct = 0.25,test.use = "wilcox",logfc.threshold = 0.25)
undiff_XX_VS_XY_DEG$cluster <- as.factor(ifelse(undiff_XX_VS_XY_DEG$p_val_adj<0.05,ifelse(undiff_XX_VS_XY_DEG$avg_log2FC>0,"XX_undiff.ed_PGCs","XY_undiff.ed_PGCs"),"unsig"))
undiff_XX_VS_XY_DEG$gene <- rownames(undiff_XX_VS_XY_DEG)
write.xlsx(undiff_XX_VS_XY_DEG,"/DataC/mmdata/mm/output/undiff_XX_VS_XY_DEG.xlsx")

# F-PGC VS M-PGC 
F_M_PGC <- subset(mitotic_PGCs,subset=celltype3 %in% c("Female mitotic PGCs","Male mitotic PGCs"))
F_M_PGC_DEG <- FindMarkers(F_M_PGC,ident.1 ="F-PGCs",ident.2 ="M-PGCs",
                           min.pct = 0.25,test.use = "wilcox",logfc.threshold = 0.25)
F_M_PGC_DEG$cluster <- as.factor(ifelse(F_M_PGC_DEG$p_val_adj<0.05,ifelse(F_M_PGC_DEG$avg_log2FC>0,"F-PGCs","M-PGCs"),"unsig"))
F_M_PGC_DEG$gene <- rownames(F_M_PGC_DEG)
write.xlsx(F_M_PGC_DEG,"/DataC/mmdata/mm/output/F_M_PGC_DEG.xlsx")
