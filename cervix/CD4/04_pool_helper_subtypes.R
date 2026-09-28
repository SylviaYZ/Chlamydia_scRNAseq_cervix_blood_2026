# pool helper subtypes

library(dplyr)
library(ggplot2)
library(Seurat)
library(scCustomize)
library(table1)
library(SingleR)
packageVersion("Seurat")
# Load data
cervix.CD4T.integrated <- readRDS("/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_CD4T_normalized_integrated_clustered_prepMarker_SingleR.RDS")

DefaultAssay(object = cervix.CD4T.integrated) <- "SCT"

table(cervix.CD4T.integrated$SingleR.pruned, cervix.CD4T.integrated$ID)

table(cervix.CD4T.integrated$ID, cervix.CD4T.integrated$condition)

DimPlot(cervix.CD4T.integrated, reduction = "umap", group.by = "integrated_snn_res.2", label = T, repel = T)

DimPlot(cervix.CD4T.integrated, reduction = "umap", group.by = "integrated_snn_res.1.2", label = T, repel = T)
# Aggregate SingleR result and visualization
cervix.CD4T.integrated$SingleR.pruned.aggregated <- ifelse(cervix.CD4T.integrated$SingleR.pruned %in% c("Th1 cells", "Th1/Th17 cells", "Th17 cells", "Th2 cells"), "T-helper cells", cervix.CD4T.integrated$SingleR.pruned)

DimPlot(cervix.CD4T.integrated, reduction = "umap", group.by = "SingleR.pruned.aggregated",
        label = FALSE)

DimPlot(cervix.CD4T.integrated, reduction = "umap", group.by = "SingleR.pruned.aggregated",
        label = TRUE)

# Save Data
saveRDS(cervix.CD4T.integrated, file = "/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_CD4T_normalized_integrated_clustered_prepMarker_SingleR_manual.RDS")

