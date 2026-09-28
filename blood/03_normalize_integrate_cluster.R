# normalize integrate cluster

library(dplyr)
library(ggplot2)
library(Seurat)
library(scCustomize)
library(table1)
packageVersion("Seurat")
# Load filtered data
blood.all <- readRDS(file = "/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/PBMC/SCTransformV2_Oct26/blood_all_filtered.RDS")
# Normalization using SCTransform V2 and PCA
## Load file for cell cycle
exp.mat <- read.table(file = "/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/nestorawa_forcellcycle_expressionMatrix.txt",
    header = TRUE, as.is = TRUE, row.names = 1)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Create list for normalization
blood.list <- SplitObject(object = blood.all, split.by = "ident")

normalize.SCT <- function(dt){

  print(dt)

    dt <- NormalizeData(dt, verbose = FALSE)
    dt <- CellCycleScoring(dt, g2m.features=g2m.genes, s.features=s.genes)
    dt <- ScaleData(dt, features = rownames(dt))

    print( FeaturePlot_scCustom(seurat_object = RunPCA(dt, features = c(s.genes, g2m.genes), verbose = FALSE), features = c('S.Score', 'G2M.Score'),  na_cutoff = 0) )

     print(paste0(unique(dt$orig.ident), ": S and G2M score correlation = ", round(cor(dt@meta.data$S.Score, dt@meta.data$G2M.Score), 5)))

    dt <- SCTransform(dt, vst.flavor = "v2", vars.to.regress = c('percent.mt', 'nCount_RNA', 'nFeature_RNA', 'S.Score', 'G2M.Score'), verbose = FALSE) %>%
    RunPCA(verbose = FALSE)

    print( FeaturePlot_scCustom(seurat_object = RunPCA(dt, features = c(s.genes, g2m.genes), verbose = FALSE), features = c('S.Score', 'G2M.Score'), na_cutoff = 0) )

  return(dt)
}

# SCT V2 normalization
blood.list <- lapply(blood.list, function(dt) normalize.SCT(dt) )

## Integration

blood.list <- blood.list[c("aggr_p1", "p2_1","p2_2")]

top.features <- SelectIntegrationFeatures(object.list = blood.list, nfeatures = 3000)

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

blood.list <- PrepSCTIntegration(object.list = blood.list, anchor.features = top.features)

blood.anchors <- FindIntegrationAnchors(object.list = blood.list,
                                         normalization.method = "SCT",
                                         anchor.features = top.features)

blood.integrated <- IntegrateData(anchorset = blood.anchors, normalization.method = "SCT")

## Dimension reduction and clustering on integrated data
DefaultAssay(object = blood.integrated) <- "integrated"

blood.integrated <- RunPCA(blood.integrated,  verbose = FALSE)

ElbowPlot(blood.integrated, ndims = 30, reduction = "pca")

print(blood.integrated[["pca"]], dims = 1:30, nfeatures = 10)

blood.integrated <- RunUMAP(blood.integrated, reduction = "pca", dims = 1:30, verbose = FALSE)

blood.integrated <- FindNeighbors(blood.integrated, reduction = "pca", dims = 1:30)

for (r in seq(0.2,2,0.2)){
    blood.integrated <- FindClusters(blood.integrated, resolution = r)
    print(
      DimPlot(blood.integrated, reduction = "umap", group.by = paste0("integrated_snn_res.", r)) + ggtitle(paste0("integrated_snn_res.", r))
    )
}

# UMAP visualization
DimPlot(blood.integrated, reduction = "umap", group.by = "ID")

DimPlot(blood.integrated, reduction = "umap", group.by = "orig.ident")

DimPlot(blood.integrated, reduction = "umap", group.by = "condition")

FeaturePlot_scCustom(blood.integrated, reduction = "umap", features = "S.Score")

FeaturePlot_scCustom(blood.integrated, reduction = "umap", features = "G2M.Score")

# PrepSCTFindMarker
blood.integrated <- PrepSCTFindMarkers(blood.integrated)

# Save data
saveRDS(blood.integrated, file = "/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/PBMC/SCTransformV2_Oct26/blood_all_filtered_normalized_integrated_clustered_prepMarker.RDS")

