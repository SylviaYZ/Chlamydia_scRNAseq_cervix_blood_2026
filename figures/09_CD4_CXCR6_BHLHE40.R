# CD4 CXCR6 and BHLHE40 expression

library(dplyr)
library(ggplot2)
library(Seurat)
library(scCustomize)
library(table1)
library(pheatmap)
library(patchwork)
library(cowplot)
library(reshape2)
library(ggpubr)
library(grid)
library(tidyr)
packageVersion("Seurat")

# load data from Jing's project space

# cervix.CD4T.integrated <- readRDS("/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_CD4T_normalized_integrated_clustered_prepMarker_SingleR_manual.RDS") 
cervix.CD4T.integrated <- readRDS("/proj/xiaojinzlab/projects/Chlamydia_scRNA-seq_processed_Seurat/cervix_CD4T_normalized_integrated_clustered_prepMarker_SingleR_manual.RDS") 

cervix.CD4T.integrated$condition <- factor(cervix.CD4T.integrated$condition,
                                           levels=c("uninfected", "infected"),
                                           labels=c("CT-","CT+"))

cervix.CD4T.integrated$SingleR.pruned <- factor(cervix.CD4T.integrated$SingleR.pruned,
                                                levels = c("Naive CD4 T cells",
                                                           "Follicular helper T cells",
                                                           "Th1 cells",
                                                           "Th1/Th17 cells",
                                                           "Th17 cells",
                                                           "Th2 cells",
                                                           "Terminal effector CD4 T cells",
                                                           "T regulatory cells"),
                                                labels = c("TN",
                                                           "Tfh",
                                                           "Th1",
                                                           "Th1/17",
                                                           "Th17",
                                                           "Th2",
                                                           "TTE",
                                                           "Treg"))

DefaultAssay(cervix.CD4T.integrated) <- "SCT"

# Subset of CT+ participants
table(cervix.CD4T.integrated$condition)

CD4_inf <- subset(
  cervix.CD4T.integrated,
  subset = condition == "CT+"
)

CD4_inf

# Gene visualization on existing FeaturePlot on UMAP
FeaturePlot(
  CD4_inf,
  features = c("CXCR6", "BHLHE40"),
  reduction = "umap",
  ncol = 2
)

# violin plot by cell type
VlnPlot(
  CD4_inf,
  features = c("CXCR6", "BHLHE40"),
  group.by = "SingleR.pruned",
  stack = TRUE,
  flip = TRUE,
  pt.size = 0
) +
  NoLegend()

# violin plot by Seurat Clusters
VlnPlot(
  CD4_inf,
  features = c("CXCR6", "BHLHE40"),
  group.by = "integrated_snn_res.2",
  stack = TRUE,
  flip = TRUE,
  pt.size = 0
) +
  NoLegend()

table(CD4_inf$integrated_snn_res.2,CD4_inf$SingleR.pruned)

# scatterplot 
scatter_data <- FetchData(
  CD4_inf,
  vars = c("CXCR6", "BHLHE40", "integrated_snn_res.2", "SingleR.pruned")
)

ggplot(
  scatter_data,
  aes(
    x = CXCR6,
    y = BHLHE40,
    color = factor(integrated_snn_res.2)
  )
) +
  geom_point(
    position = position_jitter(
      width = 0.1,
      height = 0.1,
      seed = 123
    ),
    size = 0.5,
    alpha = 0.5
  ) +
  labs(
    x = "CXCR6 normalized RNA expression",
    y = "BHLHE40 normalized RNA expression",
    color = "Cluster"
  ) +
  theme_classic()


ggplot(
  scatter_data,
  aes(
    x = CXCR6,
    y = BHLHE40,
    color = factor(SingleR.pruned)
  )
) +
  geom_point(
    position = position_jitter(
      width = 0.1,
      height = 0.1,
      seed = 123
    ),
    size = 0.5,
    alpha = 0.5
  ) +
  labs(
    x = "CXCR6 normalized RNA expression",
    y = "BHLHE40 normalized RNA expression",
    color = "Subtype"
  ) +
  theme_classic()
