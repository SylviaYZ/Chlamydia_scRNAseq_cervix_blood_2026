# SingleR annotation

library(dplyr)
library(ggplot2)
library(Seurat)
library(scCustomize)
library(table1)
library(SingleR)
packageVersion("Seurat")

# Load data
cervix.CD8T.integrated <- readRDS("/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_CD8T_normalized_integrated_clustered_prepMarker.RDS")

# Using Monaco data as reference for SingleR
monaco.ref <- celldex::MonacoImmuneData()

table(monaco.ref$label.fine, monaco.ref$label.main)

# Select CD4 T cell subtypes
monaco.CD8 <- monaco.ref[, monaco.ref$label.main == "CD8+ T cells"]

table(monaco.CD8$label.fine, monaco.CD8$label.main)
# Run SingleR on CD8 T cells
cervix.CD8.SingleR <- SingleR( GetAssayData(cervix.CD8T.integrated, assay = "RNA", layer = "counts"),
                          ref = monaco.CD8, labels = monaco.CD8$label.fine)

plotScoreHeatmap(cervix.CD8.SingleR)
# UMAP visualization
all.equal(cervix.CD8.SingleR@rownames, colnames(cervix.CD8T.integrated))

cervix.CD8T.integrated$SingleR.pruned <- cervix.CD8.SingleR$pruned.labels

DimPlot(cervix.CD8T.integrated, reduction = "umap", group.by = "SingleR.pruned",
        label = FALSE)

DimPlot(cervix.CD8T.integrated, reduction = "umap", group.by = "SingleR.pruned",
        label = TRUE)

# Save data
saveRDS(cervix.CD8T.integrated, file = "/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_CD8T_normalized_integrated_clustered_prepMarker_SingleR.RDS")
