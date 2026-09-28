# SingleR annotation

library(dplyr)
library(ggplot2)
library(Seurat)
library(scCustomize)
library(table1)
library(SingleR)
packageVersion("Seurat")
# Load data
cervix.CD4T.integrated <- readRDS("/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_CD4T_normalized_integrated_clustered_prepMarker.RDS")

# Using Monaco data as reference for SingleR
monaco.ref <- celldex::MonacoImmuneData()

table(monaco.ref$label.fine, monaco.ref$label.main)

# Select CD4 T cell subtypes
monaco.CD4 <- monaco.ref[, monaco.ref$label.main == "CD4+ T cells"]

table(monaco.CD4$label.fine, monaco.CD4$label.main)
# Run SingleR on CD4 T cells
cervix.CD4.SingleR <- SingleR( GetAssayData(cervix.CD4T.integrated, assay = "RNA", layer = "counts"),
                          ref = monaco.CD4, labels = monaco.CD4$label.fine)

plotScoreHeatmap(cervix.CD4.SingleR)
# UMAP visualization
all.equal(cervix.CD4.SingleR@rownames, colnames(cervix.CD4T.integrated))

cervix.CD4T.integrated$SingleR.pruned <- cervix.CD4.SingleR$pruned.labels

DimPlot(cervix.CD4T.integrated, reduction = "umap", group.by = "SingleR.pruned",
        label = FALSE)

DimPlot(cervix.CD4T.integrated, reduction = "umap", group.by = "SingleR.pruned",
        label = TRUE)

# Save data
saveRDS(cervix.CD4T.integrated, file = "/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_CD4T_normalized_integrated_clustered_prepMarker_SingleR.RDS")
