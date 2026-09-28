# cluster markers

library(dplyr)
library(ggplot2)
library(Seurat)
library(scCustomize)
library(table1)
packageVersion("Seurat")
# Load data
blood.integrated <- readRDS("/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/PBMC/SCTransformV2_Oct26/blood_all_filtered_normalized_integrated_clustered_prepMarker.RDS")

# Check Number of DE genes
for (r in seq(0.2,2,0.2)){
  Idents(object = blood.integrated) <- paste0("integrated_snn_res.", r)

  print(paste0("integrated_snn_res.", r ," - number of cells in each cluster:"))

  print(table(blood.integrated@active.ident))

  All.Markers <- FindAllMarkers(blood.integrated, assay = "SCT",  logfc.threshold = 1, base = 2, verbose = FALSE)

  All.Markers <- All.Markers %>% filter(p_val_adj < 0.01)

  print(paste0("integrated_snn_res.", r ," - number of DE genes w/ abs. avg_log2fc > 1 and p_val_adj < 0.01 for each cluster:"))

  print(table(All.Markers$cluster))

  write.csv(All.Markers, file = paste0("/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Blood/Blood_check_major_clusters/Wilcoxon_All_Cluster_Markers_Res_",r,".csv"))
}
