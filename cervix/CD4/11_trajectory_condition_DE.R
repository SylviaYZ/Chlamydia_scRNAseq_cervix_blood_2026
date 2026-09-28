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
cervix.CD4T.integrated <- readRDS("/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_CD4T_normalized_integrated_clustered_prepMarker_SingleR_manual.RDS")

table(cervix.CD4T.integrated$SingleR.pruned.aggregated, cervix.CD4T.integrated$ID)

table(cervix.CD4T.integrated$ID, cervix.CD4T.integrated$condition)

celltype.color <- c('TN' = "#F8766D",
                 'Tfh' = "#CD9600",
                 'Th1 cells' = '#228B22',
                 'Th17 cells' = '#800080',
                 'Th1/Th17 cells' = '#FFA500',
                 'Th2 cells' = '#FFFF00',
                 'T regulatory cells' = '#FF00FF',
                 'TTE' = "#C77CFF" ,
                 'Th' = "#00BFC4",
                 "Th1,Th1/17,Th17,Th2" = "#00BFC4")

DimPlot(cervix.CD4T.integrated, reduction = "umap", group.by = "SingleR.pruned.aggregated",
        label = T, repel = T,
        cols = celltype.color) +
  ggtitle("SingleR for CD4 T")

# Remove T reg
cervix.CD4T.integrated <- cervix.CD4T.integrated[,cervix.CD4T.integrated$SingleR.pruned.aggregated != "T regulatory cells"]

cervix.CD4T.integrated.SCE <- as.SingleCellExperiment(cervix.CD4T.integrated, assay = "SCT")

reducedDim(cervix.CD4T.integrated.SCE, "umap") <- Embeddings(cervix.CD4T.integrated, "umap")
reducedDim(cervix.CD4T.integrated.SCE, "pca") <- Embeddings(cervix.CD4T.integrated, "pca")

# Slingshot + PCA

cervix.CD4T.integrated.SCE <- slingshot(cervix.CD4T.integrated.SCE,
                                        reducedDim = 'pca',
                                        clusterLabels = colData(cervix.CD4T.integrated.SCE)$SingleR.pruned.aggregated,

                                        start.clus = "Naive CD4 T cells",
                                       )

cervix.CD4T.embed <- embedCurves(cervix.CD4T.integrated.SCE, reducedDim(cervix.CD4T.integrated.SCE, "umap"))

cervix.CD4T.curve <- slingCurves(cervix.CD4T.embed, as.df = TRUE)
cervix.CD4T.mst <- slingMST(cervix.CD4T.embed, as.df = TRUE) %>% mutate(shape = ifelse(Cluster =="Naive CD4 T cells", 17, 19))
cervix.CD4T.mst %>% dplyr::select(Lineage,Order, Cluster)

df.CD4 <- cbind(reducedDims(cervix.CD4T.integrated.SCE)$umap, cervix.CD4T.integrated.SCE$SingleR.pruned.aggregated) %>%
  as.data.frame()
colnames(df.CD4)[3] <- "SingleR.pruned.aggregated"
df.CD4$umap_1 <- as.numeric(df.CD4$umap_1)
df.CD4$umap_2 <- as.numeric(df.CD4$umap_2)

df.CD4$SingleR.pruned.aggregated <- factor(df.CD4$SingleR.pruned.aggregated,
                                                           levels = c("Naive CD4 T cells",
                                                                      "Follicular helper T cells",
                                                                      "T-helper cells",
                                                                      "Terminal effector CD4 T cells"),
                                                           labels = c("TN", "Tfh", "Th1,Th1/17,Th17,Th2", "TTE"))

p.CD4 <-  ggplot(df.CD4, aes(x = umap_1, y = umap_2)) +
    geom_point(aes(color= SingleR.pruned.aggregated),  size = 0.5) +
    theme_classic() +
  labs(color= "Cell type") +
  scale_color_manual(values=celltype.color)

# p.CD4  + geom_point(data = cervix.CD4T.mst ,
#                     shape = cervix.CD4T.mst$shape,size = 3) +
#   geom_path(data = cervix.CD4T.mst  %>% arrange(Order), aes(group = Lineage), linewidth = 1) +
#   ggtitle("Slingshot spanning trees on CD4 T cells w/o Th2, Treg")+
#   guides(color = guide_legend(override.aes = list(size=3)))

plot_traj <- p.CD4 + geom_path(data = cervix.CD4T.curve %>% arrange(Order),
              aes(group = Lineage))+
  guides(color = guide_legend(title = NULL, override.aes = list(size = 3)))+
  theme(
    text = element_text(size = 18),
    legend.text = element_text(size = 18))
     # #legend.position = "bottom",
     #  legend.position = c(0.8, 0.2))

plot_traj
ggsave("/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Cervix/Cervix_CD4T_analysis/Cervix_CD4T_trajectory/Cervix_CD4T_PCA_THelper/Visualization/CD4_traj.png",
       plot_traj,
       width = 8, height = 4, units = "in",
       dpi = 300)
# Lineage visualization

Trajectory.annotation.data.frame <- data.frame(cell_type = cervix.CD4T.integrated.SCE$SingleR.pruned.aggregated,
                       condition = cervix.CD4T.integrated.SCE$condition,
                       pseudotime = as.numeric(cervix.CD4T.integrated.SCE$slingPseudotime_1)) %>% arrange(pseudotime)

Trajectory.annotation <- HeatmapAnnotation(cell_type = Trajectory.annotation.data.frame$cell_type,
                       pseudotime = Trajectory.annotation.data.frame$pseudotime,
                       col = list(cell_type = celltype.color))

plot(Trajectory.annotation)

# DE analysis using TradeSeq
## Remove MT and RB genes, and AC genes

genes.remove <- readRDS(file = "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Cervix/Gene_Remove_DE.rds")

# Remove long integrative non-coding genes
cervix.CD4T.integrated.SCE.DE <- cervix.CD4T.integrated.SCE[!rownames(cervix.CD4T.integrated.SCE) %in% genes.remove, ]

cervix.CD4T.integrated.SCE.DE
## select knots for model
# Original notebook chunk used eval=FALSE; active here to compute the required result.
set.seed(1249)

knots <- evaluateK(counts = assays(cervix.CD4T.integrated.SCE.DE)$counts, sds = cervix.CD4T.embed, k = 3:10,
                   nGenes = 500, verbose = T)

## Fit the model
# Original notebook chunk used eval=FALSE; active here to compute the required result.
set.seed(1249)
sls.model.CD4 <- fitGAM(counts = assays(cervix.CD4T.integrated.SCE.DE)$counts, sds = cervix.CD4T.embed,
              conditions = factor(cervix.CD4T.integrated.SCE.DE$condition),
              nknots = 7, verbose = T)

saveRDS(sls.model.CD4, "/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_CD4_slingshot_tradeSeq_model.RDS")

mean(rowData(sls.model.CD4)$tradeSeq$converged)

## Condition Test
# Original notebook chunk used eval=FALSE; active here to compute the required result.
condRes <- conditionTest(sls.model.CD4, l2fc = log2(2))
condRes$padj <- p.adjust(condRes$pvalue, "fdr")
saveRDS(condRes, "cervix_CD4lineage_condTest.RDS")

condRes <- readRDS("cervix_CD4lineage_condTest.RDS")

condRes <- condRes %>% filter(padj < 0.05) %>% arrange(padj)

sum(condRes$padj <= 0.05, na.rm = TRUE)

### heatmap after condition test

yhatSmooth2 <- predictSmooth(sls.model.CD4, gene = rownames(condRes), nPoints = 50, tidy = FALSE)

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

infected.DE.Traj.cluster.df.output <- left_join(infected.DE.Traj.cluster.df,
                                                condRes %>% mutate(gene = rownames(condRes)),
                                                by = "gene")

saveRDS(infected.DE.Traj.cluster.df.output, "infected_DE_Traj_cluster.RDS")

write.csv(infected.DE.Traj.cluster.df.output, "cervix_CD4_infected_DE_Traj_cluster_output.csv")

all.equal(names(infected.DE.Traj.cluster), rownames(yhatSmooth2Scaled))

# generated infected heatmap
infected_yhatSmooth2Scaled_long <- as.data.frame(yhatSmooth2Scaled[, 1:50])

colnames(infected_yhatSmooth2Scaled_long) <- c(1:50)

infected_yhatSmooth2Scaled_long$gene <- rownames(infected_yhatSmooth2Scaled_long)

infected_yhatSmooth2Scaled_long <- gather(infected_yhatSmooth2Scaled_long, key = "pseudotime", value = "scaled_predicted_expression", -gene) %>% mutate(pseudotime = as.numeric(pseudotime))

infected_yhatSmooth2Scaled_long <- left_join(infected_yhatSmooth2Scaled_long, infected.DE.Traj.cluster.df,
                                             by = "gene")

infected_yhatSmooth2Scaled_long$cluster <- factor(infected_yhatSmooth2Scaled_long$cluster, levels=c(2,5,3,1,4))

infected_heatmap_plot <- ggplot(infected_yhatSmooth2Scaled_long,
       aes(pseudotime, y=gene, fill=scaled_predicted_expression))  +
  geom_tile() +
      scale_fill_gradient(name = "Predicted\nExpression", low="lightyellow", high="red",
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

uninfected_yhatSmooth2Scaled_long$cluster <- factor(uninfected_yhatSmooth2Scaled_long$cluster, levels=c(2,5,3,1,4))

uninfected_heatmap_plot <- ggplot(uninfected_yhatSmooth2Scaled_long,
       aes(pseudotime, y=gene, fill=scaled_predicted_expression))  +
  geom_tile() +
  scale_fill_gradient(name = "Predicted\nExpression", low="lightyellow", high="red",
                          limits=c(-3,7)) +
  facet_grid(rows=vars(cluster), space="free", scale="free")+
  ggtitle("    Uninfected") +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 18, face = "bold")
)
print(uninfected_heatmap_plot)

# plot two heatmap together
p_save <- ggarrange(infected_heatmap_plot + ggtitle("  CT+"),
  uninfected_heatmap_plot + ggtitle("  CT-"), ncol = 2,
          common.legend = TRUE,
          legend="right")

p_save
ggsave(filename = "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Cervix/Cervix_CD4T_analysis/Cervix_CD4T_trajectory/Cervix_CD4T_PCA_THelper/Visualization/combined_heatmap.png", plot = p_save, width = 8, height = 6, dpi = 300)
