# trajectory condition DE

library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(glmGamPoi)
library(knitr)
library(here)
library(slingshot)
library(gridExtra)
library(tradeSeq)
library(ComplexHeatmap)
library(tidyr)
library(ggpubr)
packageVersion("Seurat")
# Load data
cervix.CD8T.integrated <- readRDS("/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_CD8T_normalized_integrated_clustered_prepMarker_SingleR.RDS")

DefaultAssay(object = cervix.CD8T.integrated) <- "SCT"

table(cervix.CD8T.integrated$SingleR.pruned, cervix.CD8T.integrated$ID)

table(cervix.CD8T.integrated$ID, cervix.CD8T.integrated$condition)

DimPlot(cervix.CD8T.integrated, reduction = "umap", group.by = "SingleR.pruned",
        label = T, repel = T) +
  ggtitle("SingleR for CD8 T")

cervix.CD8T.integrated.SCE <- as.SingleCellExperiment(cervix.CD8T.integrated, assay = "SCT")

reducedDim(cervix.CD8T.integrated.SCE, "umap") <- Embeddings(cervix.CD8T.integrated, "umap")
reducedDim(cervix.CD8T.integrated.SCE, "pca") <- Embeddings(cervix.CD8T.integrated, "pca")

# Slingshot + PCA

cervix.CD8T.integrated.SCE <- slingshot(cervix.CD8T.integrated.SCE,
                                        reducedDim = 'pca',
                                        clusterLabels = colData(cervix.CD8T.integrated.SCE)$SingleR.pruned,

                                        start.clus = "Naive CD8 T cells",
                                       )

cervix.CD8T.embed <- embedCurves(cervix.CD8T.integrated.SCE, reducedDim(cervix.CD8T.integrated.SCE, "umap"))

cervix.CD8T.curve <- slingCurves(cervix.CD8T.embed, as.df = TRUE)
cervix.CD8T.mst <- slingMST(cervix.CD8T.embed, as.df = TRUE) %>% mutate(shape = ifelse(Cluster =="Naive CD8 T cells", 17, 19))
cervix.CD8T.mst %>% dplyr::select(Lineage,Order, Cluster)

df.CD8 <- cbind(reducedDims(cervix.CD8T.integrated.SCE)$umap, cervix.CD8T.integrated.SCE$SingleR.pruned) %>%
  as.data.frame()
colnames(df.CD8)[3] <- "SingleR.pruned"
df.CD8$umap_1 <- as.numeric(df.CD8$umap_1)
df.CD8$umap_2 <- as.numeric(df.CD8$umap_2)

p.CD8 <-  ggplot(df.CD8, aes(x = umap_1, y = umap_2)) +
    geom_point(aes(color= SingleR.pruned),  size = 0.5) +
    theme_classic() +
  labs(color= "Cell type")

p.CD8 + geom_path(data = cervix.CD8T.curve %>% arrange(Order),
              aes(group = Lineage))+
  ggtitle("Slingshot trajectory on CD8 T cells when given cell types") +
  guides(color = guide_legend(override.aes = list(size=3)))

# Lineage visualization

celltype.color <- c('Naive CD8 T cells' = '#FF0000',
                 'Central memory CD8 T cells' = '#0000FF',
                 'Effector memory CD8 T cells' = '#228B22',
                 'Terminal effector CD8 T cells' = '#00FFFF')

Trajectory.annotation.data.frame <- data.frame(cell_type = cervix.CD8T.integrated.SCE$SingleR.pruned,
                       condition = cervix.CD8T.integrated.SCE$condition,
                       pseudotime = as.numeric(cervix.CD8T.integrated.SCE$slingPseudotime_1)) %>% arrange(pseudotime)

Trajectory.annotation <- HeatmapAnnotation(cell_type = Trajectory.annotation.data.frame$cell_type,
                       pseudotime = Trajectory.annotation.data.frame$pseudotime,
                       col = list(cell_type = celltype.color))

plot(Trajectory.annotation)

# DE analysis using TradeSeq
## Remove MT and RB genes, and AC genes

genes.remove <- readRDS(file = "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Cervix/Gene_Remove_DE.rds")

# Remove long integrative non-coding genes
cervix.CD8T.integrated.SCE.DE <- cervix.CD8T.integrated.SCE[!rownames(cervix.CD8T.integrated.SCE) %in% genes.remove, ]

cervix.CD8T.integrated.SCE.DE
## select knots for model
# Original notebook chunk used eval=FALSE; active here to compute the required result.
set.seed(425)

knots <- evaluateK(counts = assays(cervix.CD8T.integrated.SCE.DE)$counts, sds = cervix.CD8T.embed, k = 3:10,
                   nGenes = 500, verbose = T)

## Fit the model
# Original notebook chunk used eval=FALSE; active here to compute the required result.
set.seed(425)
sls.model.CD8 <- fitGAM(counts = assays(cervix.CD8T.integrated.SCE.DE)$counts, sds = cervix.CD8T.embed,
              conditions = factor(cervix.CD8T.integrated.SCE.DE$condition),
              nknots = 6, verbose = T)

saveRDS(sls.model.CD8, "/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_CD8_slingshot_tradeSeq_model.RDS")

mean(rowData(sls.model.CD8)$tradeSeq$converged)

sls.model.CD8 <- readRDS(file = "/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_CD8_slingshot_tradeSeq_model.RDS")
## Condition Test
# Original notebook chunk used eval=FALSE; active here to compute the required result.
condRes <- conditionTest(sls.model.CD8, l2fc = log2(2))
condRes$padj <- p.adjust(condRes$pvalue, "fdr")
saveRDS(condRes, "cervix_CD8lineage_condTest.RDS")

condRes <- readRDS("cervix_CD8lineage_condTest.RDS")

condRes <- condRes %>% filter(padj < 0.05) %>% arrange(padj)

sum(condRes$padj <= 0.05, na.rm = TRUE)

### heatmap after condition test

yhatSmooth2 <- predictSmooth(sls.model.CD8, gene = rownames(condRes), nPoints = 50, tidy = FALSE)

# Scale data
yhatSmooth2Scaled <- t(scale(t(yhatSmooth2)))

# Heatmap for temporary visualization
heatSmooth2_infected <- heatmap(yhatSmooth2Scaled[, 1:50],
scale = "none", keep.dendro=TRUE,
 Rowv = TRUE, Colv = FALSE
)

# Get cluster from heatmap
heatSmooth2_infected_row.clusters <-  as.hclust(heatSmooth2_infected$Rowv)

infected.DE.Traj.cluster <- cutree(heatSmooth2_infected_row.clusters, k=5)

infected.DE.Traj.cluster.df <- data.frame(gene = names(infected.DE.Traj.cluster),
                                          cluster = infected.DE.Traj.cluster)

saveRDS(infected.DE.Traj.cluster.df, "infected_DE_Traj_cluster.RDS")

all.equal(names(infected.DE.Traj.cluster), rownames(yhatSmooth2Scaled))

# generated infected heatmap
infected_yhatSmooth2Scaled_long <- as.data.frame(yhatSmooth2Scaled[, 1:50])

colnames(infected_yhatSmooth2Scaled_long) <- c(1:50)

infected_yhatSmooth2Scaled_long$gene <- rownames(infected_yhatSmooth2Scaled_long)

infected_yhatSmooth2Scaled_long <- gather(infected_yhatSmooth2Scaled_long, key = "pseudotime", value = "scaled_predicted_expression", -gene) %>% mutate(pseudotime = as.numeric(pseudotime))

infected_yhatSmooth2Scaled_long <- left_join(infected_yhatSmooth2Scaled_long, infected.DE.Traj.cluster.df,
                                             by = "gene")

infected_heatmap_plot <- ggplot(infected_yhatSmooth2Scaled_long,
       aes(pseudotime, y=gene, fill=scaled_predicted_expression))  +
  geom_tile() +
      scale_fill_gradient(low="white", high="darkgreen",
                          limits=c(-3,7)) +
  facet_grid(rows=vars(cluster), space="free", scale="free") +
  ggtitle("     Infected") +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 18, face = "bold"))

print(infected_heatmap_plot)

# Generated uninfected heatmap, gene group is from infected data cluster

uninfected_yhatSmooth2Scaled_long <- as.data.frame(yhatSmooth2Scaled[, 51:100])

colnames(uninfected_yhatSmooth2Scaled_long) <- c(1:50)

uninfected_yhatSmooth2Scaled_long$gene <- rownames(uninfected_yhatSmooth2Scaled_long)

uninfected_yhatSmooth2Scaled_long <- gather(uninfected_yhatSmooth2Scaled_long, key = "pseudotime", value = "scaled_predicted_expression", -gene) %>% mutate(pseudotime = as.numeric(pseudotime))

uninfected_yhatSmooth2Scaled_long <- left_join(uninfected_yhatSmooth2Scaled_long, infected.DE.Traj.cluster.df,
                                             by = "gene")

uninfected_heatmap_plot <- ggplot(uninfected_yhatSmooth2Scaled_long,
       aes(pseudotime, y=gene, fill=scaled_predicted_expression))  +
  geom_tile() +
  scale_fill_gradient(low="white", high="darkgreen",
                          limits=c(-3,7)) +
  facet_grid(rows=vars(cluster), space="free", scale="free")+
  ggtitle("    Uninfected") +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 18, face = "bold"))
print(uninfected_heatmap_plot)

# plot two heatmap together
ggarrange(infected_heatmap_plot, uninfected_heatmap_plot, ncol = 2,
          common.legend = TRUE,
          legend="right")

