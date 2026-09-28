# normalize integrate cluster

library(dplyr)
library(ggplot2)
library(Seurat)
library(scCustomize)
library(table1)
packageVersion("Seurat")
# Load filtered data
cervix.all <- readRDS(file = "/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_all_filtered.RDS")
# Normalization using SCTransform V2 and PCA
## Load file for cell cycle
exp.mat <- read.table(file = "/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/nestorawa_forcellcycle_expressionMatrix.txt",
    header = TRUE, as.is = TRUE, row.names = 1)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Create list for normalization
cervix.list <- SplitObject(object = cervix.all, split.by = "ident")

normalize.SCT <- function(dt){

  print(dt)

  if(length(unique(dt$orig.ident)) == 1){

    dt <- NormalizeData(dt, verbose = FALSE)
    dt <- CellCycleScoring(dt, g2m.features=g2m.genes, s.features=s.genes)
    dt <- ScaleData(dt, features = rownames(dt))

    print( FeaturePlot_scCustom(seurat_object = RunPCA(dt, features = c(s.genes, g2m.genes), verbose = FALSE), features = c('S.Score', 'G2M.Score'),  na_cutoff = 0) )

     print(paste0(unique(dt$orig.ident), ": S and G2M score correlation = ", round(cor(dt@meta.data$S.Score, dt@meta.data$G2M.Score), 5)))

    dt <- SCTransform(dt, vst.flavor = "v2", vars.to.regress = c('percent.mt', 'nCount_RNA', 'nFeature_RNA', 'S.Score', 'G2M.Score'), verbose = FALSE) %>%
    RunPCA(npacs = 30, verbose = FALSE)

    print( FeaturePlot_scCustom(seurat_object = RunPCA(dt, features = c(s.genes, g2m.genes), verbose = FALSE), features = c('S.Score', 'G2M.Score'), na_cutoff = 0) )

  }
  else if (length(unique(dt$ID)) > 1){

    dt <- NormalizeData(dt, verbose = FALSE)
    dt <- CellCycleScoring(dt, g2m.features=g2m.genes, s.features=s.genes)
    dt <- ScaleData(dt, features = rownames(dt))

    print( FeaturePlot_scCustom(seurat_object = RunPCA(dt, features = c(s.genes, g2m.genes), verbose = FALSE), features = c('S.Score', 'G2M.Score'),  na_cutoff = 0) )

     print(paste0(unique(dt$orig.ident), ": S and G2M score correlation = ", round(cor(dt@meta.data$S.Score, dt@meta.data$G2M.Score), 5)))

    dt <- SCTransform(dt, vst.flavor = "v2", vars.to.regress = c('percent.mt', 'nCount_RNA', 'nFeature_RNA', 'S.Score', 'G2M.Score', 'ID'), verbose = FALSE) %>%
    RunPCA(npacs = 30, verbose = FALSE)

    print( FeaturePlot_scCustom(seurat_object = RunPCA(dt, features = c(s.genes, g2m.genes), verbose = FALSE), features = c('S.Score', 'G2M.Score'), na_cutoff = 0) )
    }

  return(dt)
}

# SCT V2 normalization
cervix.list <- lapply(cervix.list, function(dt) normalize.SCT(dt) )

## Integration

cervix.list <- cervix.list[c("aggr_cx1", "cx2", "cx3_1", "cx3_2", "aggr_cx4", "aggr_cx5", "cx6_1")]

top.features <- SelectIntegrationFeatures(object.list = cervix.list, nfeatures = 3000)

# Check Ribosomal genes in Top Features: 81 RB genes
top.features[grepl("^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", top.features)]

# Check Mitachondrial genes in Top Features: 13 MT genes
top.features[grepl("^MT-", top.features)]

# Check TCR in Top Features
top.features[grepl("^TRAV|^TRBV|^TRDV|^TRGV|^TRAC|^TRBC", top.features)]

# Check TCR in Top Features
top.features[grepl("^IG[HKL]V|^IG[KL]J|^IG[KL]C|^IGH[ADEGM]", top.features)]

# Remove MT and RB genes
top.features <- top.features[!grepl("^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", top.features) & !grepl("^MT-", top.features)]

length(top.features)

cervix.list <- PrepSCTIntegration(object.list = cervix.list, anchor.features = top.features)

cervix.anchors <- FindIntegrationAnchors(object.list = cervix.list,
                                         normalization.method = "SCT",
                                         anchor.features = top.features,
                                         k.filter = NA)

cervix.integrated <- IntegrateData(anchorset = cervix.anchors, normalization.method = "SCT")

## Dimension reduction and clustering on integrated data
DefaultAssay(object = cervix.integrated) <- "integrated"

cervix.integrated <- RunPCA(cervix.integrated, dims = 1:30, verbose = FALSE)

ElbowPlot(cervix.integrated, ndims = 30, reduction = "pca")

print(cervix.integrated[["pca"]], dims = 1:30, nfeatures = 10)

cervix.integrated <- RunUMAP(cervix.integrated, reduction = "pca", dims = 1:30, verbose = FALSE)

cervix.integrated <- FindNeighbors(cervix.integrated, reduction = "pca", dims = 1:30)

for (r in seq(0.2,2,0.2)){
    cervix.integrated <- FindClusters(cervix.integrated, resolution = r)
    print(
      DimPlot(cervix.integrated, reduction = "umap", group.by = paste0("integrated_snn_res.", r)) + ggtitle(paste0("integrated_snn_res.", r))
    )
}

# UMAP visualization
DimPlot(cervix.integrated, reduction = "umap", group.by = "ID")

DimPlot(cervix.integrated, reduction = "umap", group.by = "orig.ident")

DimPlot(cervix.integrated, reduction = "umap", group.by = "condition")

# PrepSCTFindMarker
cervix.integrated <- PrepSCTFindMarkers(cervix.integrated)

# Save data
saveRDS(cervix.integrated, file = "/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_all_filtered_normalized_integrated_clustered_prepMarker.RDS")

