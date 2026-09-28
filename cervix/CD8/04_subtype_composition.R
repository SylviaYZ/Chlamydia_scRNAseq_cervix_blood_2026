# subtype composition

library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(knitr)
library(table1)
library(reshape2)
library(cowplot)
library(ggpubr)
packageVersion("Seurat")
# Load data
cervix.CD8T.integrated <- readRDS("/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_CD8T_normalized_integrated_clustered_prepMarker_SingleR.RDS")

# Observed cell counts

cell_counts <- function(cluster){
  test_counts <- table( ifelse(Idents(cervix.CD8T.integrated) == cluster, paste0("InCluster"), paste0("NotInCluster" )), cervix.CD8T.integrated$condition)
  test.data <- as.data.frame(test_counts)
   rownames(test.data) <- paste0(test.data$Var1,"_",test.data$Var2)
  test.data["Cluster",] <- c(NA,NA, cluster)

  df_ret <- as.data.frame(test.data$Freq)
  rownames(df_ret) <- rownames(test.data)
  return(t(df_ret))
}

## UMAP visialization of INF vs UNINF within each cell type

table(cervix.CD8T.integrated$SingleR.pruned, cervix.CD8T.integrated$condition)

DimPlot(cervix.CD8T.integrated, reduction = "umap", split.by = "condition", label = T, pt.size = 0.1, group.by = "SingleR.pruned")

## Cell composition for INF vs UNINF within each cell type
Idents(cervix.CD8T.integrated) <- factor(cervix.CD8T.integrated$SingleR.pruned)

Chi2_CD8T <- t(unlist(sapply(unique(cervix.CD8T.integrated$SingleR.pruned),
                          function(i) cell_counts(cluster = i)))) %>% as.data.frame()

colnames(Chi2_CD8T) <- c("InCluster_infected", "NotInCluster_infected",
                               "InCluster_uninfected", "NotInCluster_uninfected",
                               "cell_type")

Chi2_CD8T <- Chi2_CD8T %>% dplyr::select(cell_type, InCluster_infected,
                                              InCluster_uninfected, NotInCluster_infected,
                                              NotInCluster_uninfected) %>%
  mutate(InCluster_infected = as.numeric(InCluster_infected),
         InCluster_uninfected = as.numeric(InCluster_uninfected),
         NotInCluster_infected = as.numeric(NotInCluster_infected),
          NotInCluster_uninfected = as.numeric( NotInCluster_uninfected),
    inf_total = sum(as.numeric(InCluster_infected, NotInCluster_infected)), uninf_total = sum(as.numeric(InCluster_uninfected, NotInCluster_uninfected)),
    inf_proportion = InCluster_infected / inf_total,
    uninf_proportion = InCluster_uninfected / uninf_total)

Chi2_CD8T

write.csv(Chi2_CD8T, file = "Cervix_INFvsUNINF_CD8T_SingleR_Chi2Enrich.csv")

## Proporiton INF vs UNINF within each cell type visualization
inf_byCelltype_plot <- Chi2_CD8T %>%
  mutate(infected = InCluster_infected/inf_total,
         uninfected = InCluster_uninfected / uninf_total) %>%
  dplyr::select(cell_type, infected, uninfected)

inf_byCelltype_plot <- melt(inf_byCelltype_plot, id.vars = "cell_type")
inf_byCelltype_plot$value <- as.numeric(inf_byCelltype_plot$value)
inf_byCelltype_plot$posting <- as.numeric(inf_byCelltype_plot$value)*100

ggplot(inf_byCelltype_plot, aes(y=variable, x=cell_type)) +
      geom_tile(aes(fill=value)) +
      scale_fill_gradient(low="white", high="red",
                          limits=c(0,1)) +
      geom_text(aes(label=round(posting, digits = 1)), size = 7) +
      labs(title=str_wrap("Cervix %infected vs. %uninfected cells within each cell type", 60),
           x = "",
           y = "") +
      theme(text = element_text(size = 15),
            axis.text.x = element_text(size = 15,angle = 45, hjust=1,vjust=1),
            axis.text.y = element_text(size = 15)) +
      coord_cartesian(xlim = c(1,4), ylim = c(1,2), clip = "off")
