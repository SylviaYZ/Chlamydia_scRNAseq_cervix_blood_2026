# subset normalize integrate

library(dplyr)
library(ggplot2)
library(Seurat)
library(scCustomize)
library(table1)
packageVersion("Seurat")
# Load data and select CD4T
cervix.integrated <- readRDS("/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_all_filtered_normalized_integrated_clustered_prepMarker_majorCT.RDS")

cervix.CD4T <- cervix.integrated[,cervix.integrated$major_celltypes== "CD4_T"]

cervix.CD4T

table( cervix.CD4T$ID,cervix.CD4T$orig.ident, cervix.CD4T$major_celltypes)

cervix.CD4T <- cervix.CD4T[, ! cervix.CD4T$orig.ident %in% c("aggr_cx1", "cx3_1")]

table(cervix.CD4T$major_celltypes, cervix.CD4T$orig.ident)

# Normalization & Integration
cervix.CD4T.list <- SplitObject(cervix.CD4T, split.by = "orig.ident")

normalize.SCT <- function(dt){

  if(length(unique(dt$orig.ident)) == 1){

    dt <- SCTransform(dt, vst.flavor = "v2", vars.to.regress = c('percent.mt', 'nCount_RNA', 'nFeature_RNA', 'S.Score', 'G2M.Score'), verbose = FALSE) %>%
    RunPCA(npacs = 30, verbose = FALSE)

  }
  else if (length(unique(dt$ID)) > 1){

    dt <- SCTransform(dt, vst.flavor = "v2", vars.to.regress = c('percent.mt', 'nCount_RNA', 'nFeature_RNA', 'S.Score', 'G2M.Score', 'ID'), verbose = FALSE) %>%
    RunPCA(npacs = 30, verbose = FALSE) }

  return(dt)
}

cervix.CD4T.list <- lapply(cervix.CD4T.list,
                           function(dt) normalize.SCT(dt))

CD4.features  <- SelectIntegrationFeatures(cervix.CD4T.list, nfeatures = 3000)

# Check MT genes
CD4.features[grepl("^MT-", CD4.features)]
# Remove MT genes
CD4.features <- CD4.features[!grepl("^MT-", CD4.features)]

# Check RB genes
CD4.features[grepl("^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", CD4.features)]
# Remove RB genes
CD4.features <- CD4.features[!grepl("^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", CD4.features)]

# Check TCR genes
CD4.features[grepl("^TRAV|^TRBV|^TRDV|^TRGV|^TRAC|^TRBC", CD4.features)]
# Remove TCF genes
CD4.features <- CD4.features[!grepl("^TRAV|^TRBV|^TRDV|^TRGV|^TRAC|^TRBC", CD4.features)]

cervix.CD4T.list <- PrepSCTIntegration(cervix.CD4T.list, anchor.features = CD4.features)

CD4.anchors <- FindIntegrationAnchors(cervix.CD4T.list, normalization.method = "SCT", anchor.features = CD4.features)

cervix.CD4T.integrated <- IntegrateData(anchorset = CD4.anchors, normalization.method = "SCT")

cervix.CD4T.integrated <- RunPCA(cervix.CD4T.integrated)

ElbowPlot(cervix.CD4T.integrated, ndims = 30, reduction = "pca")

cervix.CD4T.integrated <- RunUMAP(cervix.CD4T.integrated, reduction = "pca", dims = 1:30)

cervix.CD4T.integrated <- FindNeighbors(cervix.CD4T.integrated, dims = 1:30)

for (r in seq(0.2,2,0.2)){
    cervix.CD4T.integrated <- FindClusters(cervix.CD4T.integrated, resolution = r)
    print(
      DimPlot(cervix.CD4T.integrated, reduction = "umap", group.by = paste0("integrated_snn_res.", r)) + ggtitle(paste0("integrated_snn_res.", r))
    )
}

# UMAP visualizaiton
DimPlot(cervix.CD4T.integrated, reduction = "umap", group.by = "ID")

DimPlot(cervix.CD4T.integrated, reduction = "umap", group.by = "orig.ident")

DimPlot(cervix.CD4T.integrated, reduction = "umap", group.by = "condition")
# PrepSCTFindMarker
cervix.CD4T.integrated <- PrepSCTFindMarkers(cervix.CD4T.integrated)

# Save data
# Original notebook chunk used eval=FALSE; active here to compute the required result.
saveRDS(cervix.CD4T.integrated, file = "/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_CD4T_normalized_integrated_clustered_prepMarker.RDS")
