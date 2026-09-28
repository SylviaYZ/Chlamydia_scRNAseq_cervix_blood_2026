# subset normalize integrate

library(dplyr)
library(ggplot2)
library(Seurat)
library(scCustomize)
library(table1)
packageVersion("Seurat")
# Load data and select CD8 T
cervix.integrated <- readRDS("/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_all_filtered_normalized_integrated_clustered_prepMarker_majorCT.RDS")

cervix.CD8T <- cervix.integrated[,cervix.integrated$major_celltypes== "CD8_T"]

cervix.CD8T

table( cervix.CD8T$ID,cervix.CD8T$orig.ident, cervix.CD8T$major_celltypes)

cervix.CD8T <- cervix.CD8T[, ! cervix.CD8T$orig.ident %in% c("aggr_cx1", "cx3_1")]

table(cervix.CD8T$major_celltypes, cervix.CD8T$orig.ident)

# Normalization & Integration
cervix.CD8T.list <- SplitObject(cervix.CD8T, split.by = "orig.ident")

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

cervix.CD8T.list <- lapply(cervix.CD8T.list,
                           function(dt) normalize.SCT(dt))

CD8.features  <- SelectIntegrationFeatures(cervix.CD8T.list, nfeatures = 3000)

# Check MT genes
CD8.features[grepl("^MT-", CD8.features)]
# Remove MT genes
CD8.features <- CD8.features[!grepl("^MT-", CD8.features)]

# Check RB genes
CD8.features[grepl("^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", CD8.features)]
# Remove RB genes
CD8.features <- CD8.features[!grepl("^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", CD8.features)]

# Check TCR genes
CD8.features[grepl("^TRAV|^TRBV|^TRDV|^TRGV|^TRAC|^TRBC", CD8.features)]
# Remove TCF genes
CD8.features <- CD8.features[!grepl("^TRAV|^TRBV|^TRDV|^TRGV|^TRAC|^TRBC", CD8.features)]

cervix.CD8T.list <- PrepSCTIntegration(cervix.CD8T.list, anchor.features = CD8.features)

CD8.anchors <- FindIntegrationAnchors(cervix.CD8T.list, normalization.method = "SCT", anchor.features = CD8.features,
                                         k.filter = NA)

cervix.CD8T.integrated <- IntegrateData(anchorset = CD8.anchors, normalization.method = "SCT", k.weight = 50)

cervix.CD8T.integrated <- RunPCA(cervix.CD8T.integrated)

ElbowPlot(cervix.CD8T.integrated, ndims = 30, reduction = "pca")

cervix.CD8T.integrated <- RunUMAP(cervix.CD8T.integrated, reduction = "pca", dims = 1:30)

cervix.CD8T.integrated <- FindNeighbors(cervix.CD8T.integrated, dims = 1:30)

for (r in seq(0.2,2,0.2)){
    cervix.CD8T.integrated <- FindClusters(cervix.CD8T.integrated, resolution = r)
    print(
      DimPlot(cervix.CD8T.integrated, reduction = "umap", group.by = paste0("integrated_snn_res.", r)) + ggtitle(paste0("integrated_snn_res.", r))
    )
}

# UMAP visualizaiton
DimPlot(cervix.CD8T.integrated, reduction = "umap", group.by = "ID")

DimPlot(cervix.CD8T.integrated, reduction = "umap", group.by = "orig.ident")

DimPlot(cervix.CD8T.integrated, reduction = "umap", group.by = "condition")
# PrepSCTFindMarker
cervix.CD8T.integrated <- PrepSCTFindMarkers(cervix.CD8T.integrated)

# Save data
saveRDS(cervix.CD8T.integrated, file = "/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_CD8T_normalized_integrated_clustered_prepMarker.RDS")
