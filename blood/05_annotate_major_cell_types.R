# annotate major cell types

library(dplyr)
library(ggplot2)
library(Seurat)
library(scCustomize)
library(table1)
packageVersion("Seurat")
# Load data

blood.integrated <- readRDS("/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/PBMC/SCTransformV2_Oct26/blood_all_filtered_normalized_integrated_clustered_prepMarker.RDS")

DefaultAssay(object = blood.integrated) <- "SCT"
# UMAP
DimPlot(blood.integrated, reduction = "umap", group.by ="integrated_snn_res.1.8", label = TRUE)

DimPlot(blood.integrated, reduction = "umap", group.by ="integrated_snn_res.2", label = TRUE)
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

# T cells
FeaturePlot(object = blood.integrated, features = T.cells)

## CD4 vs CD8 T cells
FeaturePlot(object = blood.integrated, features = CD4.CD8.T)
## gd T cells
FeaturePlot(object = blood.integrated, features = gd.T)

# NK / NK T cells T cells
FeaturePlot(object = blood.integrated, features = NK)
# B/plasma cells
## B cells
FeaturePlot(object = blood.integrated, features = B)

FeaturePlot(object = blood.integrated, features = B.plasma)

## Plasma cells
FeaturePlot(object = blood.integrated, features = Plasma)
FeaturePlot(object = blood.integrated, features = Plasma2)
# T cell subtypes
## CD4 T reg
FeaturePlot(object = blood.integrated, features = CD4.Treg1)

FeaturePlot(object = blood.integrated, features = CD4.Treg2)

## CM vs. EM vs. Naive
FeaturePlot(object = blood.integrated, features = T.sub)

## MAIT
FeaturePlot(object = blood.integrated, features = MAIT1)
FeaturePlot(object = blood.integrated, features = MAIT2)
# Define cell types
DimPlot(blood.integrated, reduction = "umap", group.by ="integrated_snn_res.2", label = TRUE) + ggtitle("integrated_snn_res.2")

blood.integrated$major_celltypes <- "Undecided"

blood.integrated$major_celltypes[blood.integrated$integrated_snn_res.2 %in% c(15,1,7,11,10,0,19,20,24,3,16,4,2)] <- "CD4_T"

blood.integrated$major_celltypes[blood.integrated$integrated_snn_res.2 %in% c(9,8,12,17,5,22,18,21,14)] <- "CD8_T"

blood.integrated$major_celltypes[blood.integrated$integrated_snn_res.2 %in% c(6)] <- "NK_T"

blood.integrated$major_celltypes[blood.integrated$integrated_snn_res.2 %in% c(23,13)] <- "gd_T"

DimPlot(blood.integrated, reduction = "umap", group.by ="major_celltypes", label = TRUE) + ggtitle("Major Cell Types")
# Save data
saveRDS(blood.integrated,
        "/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/PBMC/SCTransformV2_Oct26/blood_all_filtered_normalized_integrated_clustered_prepMarker_majorCT.RDS")
