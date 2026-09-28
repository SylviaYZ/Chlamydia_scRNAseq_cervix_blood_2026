# cluster markers

library(dplyr)
library(ggplot2)
library(Seurat)
library(scCustomize)
library(table1)
library(SingleR)

packageVersion("Seurat")
# Load data
cervix.CD8T.integrated <- readRDS("/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_CD8T_normalized_integrated_clustered_prepMarker_SingleR.RDS")

DefaultAssay(object = cervix.CD8T.integrated) <- "SCT"

DimPlot(cervix.CD8T.integrated, reduction = "umap", group.by = "integrated_snn_res.2", label = T, repel = T)

# Check Number of DE genes
for (r in seq(0.2,2,0.2)){
  Idents(object = cervix.CD8T.integrated) <- paste0("integrated_snn_res.", r)

  print(paste0("integrated_snn_res.", r ," - number of cells in each cluster:"))

  print(table(cervix.CD8T.integrated@active.ident))

  All.Markers <- FindAllMarkers(cervix.CD8T.integrated, assay = "SCT",  logfc.threshold = 1, base = 2, verbose = FALSE)

  All.Markers <- All.Markers %>% filter(p_val_adj < 0.01)

  print(paste0("integrated_snn_res.", r ," - number of DE genes w/ abs. avg_log2fc > 1 and p_val_adj < 0.01 for each cluster:"))

  print(table(All.Markers$cluster))

  write.csv(All.Markers, file = paste0("Wilcoxon_CD8T_Cluster_Markers_Res_",r,".csv"))
}
