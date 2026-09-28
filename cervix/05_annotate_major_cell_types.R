# annotate major cell types

library(dplyr)
library(ggplot2)
library(Seurat)
library(scCustomize)
library(table1)
packageVersion("Seurat")
# Load data

cervix.integrated <- readRDS("/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_all_filtered_normalized_integrated_clustered_prepMarker.RDS")

DefaultAssay(object = cervix.integrated) <- "SCT"
# UMAP
DimPlot(cervix.integrated, reduction = "umap", group.by ="integrated_snn_res.1.8", label = TRUE)

DimPlot(cervix.integrated, reduction = "umap", group.by ="integrated_snn_res.2", label = TRUE)
# Define major cell type markers
T.cells <- c("CD3D", "CD3E", "CD3G", "CD27", "TRAC", "TRBC1")
CD4.CD8.T <- c("SELL","CD4", "LDHB", "IL7R", "CD8A", "CD8B")
CD4.Treg1 <- c("IL7R", "IL2RA",  "IL2RB","FOXP3", "SATB1", "LGALS3")
CD4.Treg2 <- c("CTLA4", "ENTPD1", "NT5E", "TNFRSF18", "LAG3")
T.sub <- c("CCR7", "ID2", "CCR7", "CD27")

MAIT1 <- c("KLRB1", "IL18R1", "CCR5", "CCR6", "CXCR6", "TRAV1-2")

MAIT2 <- c( "SLC4A10", "MAF", "PRSS35", "CXCR4")

gd.T <- c("TRGC2","TRDC", "TRGC1")
NK <- c("NKG7", "NKTR","GNLY", "KLRB1","KLRD1", "KLRG1")
B <- c("CD79A", "CD79B", "MS4A1", "CD19")
Plasma <- c("IRF4", "XBP1", "CXCR4", "SDC1","KLF4", "TNFRSF17")
Plasma2 <- c("CD27", "CD19", "CD38", "CXCR4", "BCMA")
B.plasma <- c("IGHD","IGHM", "IGLC1", "IGKC")
Myeloid1 <- c("CD14", "FCGR3B", "FCGR3A", "CCR2", "CCR5", "HLA-DRA")
Myeloid2 <- c( "S100A8", "S100A10", "CD86", "LYZ", "CD74", "CD68")
Progenitor <- c("PTPRC", "CD34")

# T cells
FeaturePlot(object = cervix.integrated, features = T.cells)

## CD4 vs CD8 T cells
FeaturePlot(object = cervix.integrated, features = CD4.CD8.T)
## gd T cells
FeaturePlot(object = cervix.integrated, features = gd.T)

# NK / NK T cells T cells
FeaturePlot(object = cervix.integrated, features = NK)
# B/plasma cells
## B cells
FeaturePlot(object = cervix.integrated, features = B)

FeaturePlot(object = cervix.integrated, features = B.plasma)

## Plasma cells
FeaturePlot(object = cervix.integrated, features = Plasma)
FeaturePlot(object = cervix.integrated, features = Plasma2)
## B/plasma cells
# Original notebook chunk used eval=FALSE; active here to compute the required result.
FeaturePlot(object = cervix.integrated, features = B.plasma)

# CD14+ and CD16+ Myeloid cells
FeaturePlot(object = cervix.integrated, features = Myeloid1)

FeaturePlot(object = cervix.integrated, features = Myeloid2)

# Progenitors
FeaturePlot(object = cervix.integrated, features = Progenitor)

# T cell subtypes
## CD4 T reg
FeaturePlot(object = cervix.integrated, features = CD4.Treg1)

FeaturePlot(object = cervix.integrated, features = CD4.Treg2)

## CM vs. EM vs. Naive
FeaturePlot(object = cervix.integrated, features = T.sub)

## MAIT
FeaturePlot(object = cervix.integrated, features = MAIT1)
FeaturePlot(object = cervix.integrated, features = MAIT2)
# Define cell types
DimPlot(cervix.integrated, reduction = "umap", group.by ="integrated_snn_res.2", label = TRUE) + ggtitle("integrated_snn_res.2")

cervix.integrated$major_celltypes <- "Undecided"
cervix.integrated$major_celltypes[cervix.integrated$integrated_snn_res.2 %in% c(4,16)] <- "B"
cervix.integrated$major_celltypes[cervix.integrated$integrated_snn_res.2 %in% c(5, 7, 11, 13, 18)] <- "CD8_T"
cervix.integrated$major_celltypes[cervix.integrated$integrated_snn_res.2 %in% c(0, 1, 3, 6, 8, 9, 17, 20)] <- "CD4_T"
cervix.integrated$major_celltypes[cervix.integrated$integrated_snn_res.2 %in% c(10, 14)] <- "NK"
cervix.integrated$major_celltypes[cervix.integrated$integrated_snn_res.2 %in% c(15, 19)] <- "CD14_Mono"
cervix.integrated$major_celltypes[cervix.integrated$integrated_snn_res.2 %in% c(12)] <- "CD16_Mono"
cervix.integrated$major_celltypes[cervix.integrated$integrated_snn_res.2 %in% c(2)] <- "MAIT"

DimPlot(cervix.integrated, reduction = "umap", group.by ="major_celltypes", label = TRUE) + ggtitle("Major Cell Types")
# Save data
saveRDS(cervix.integrated,
        "/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_all_filtered_normalized_integrated_clustered_prepMarker_majorCT.RDS")
