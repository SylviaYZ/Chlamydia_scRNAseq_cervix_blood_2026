# collectri viper inference

library(dplyr)
library(ggplot2)
library(Seurat)
library(scCustomize)
library(table1)
library(ggrepel)
library(decoupleR)
library(OmnipathR)
library(tibble)
library(tidyr)
library(patchwork)
library(pheatmap)
library(EnhancedVolcano)
library(xml2)
library(purrr)
library(viper)
packageVersion("Seurat")

###################################
# Load data
##################################

cervix.integrated <- readRDS("/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_all_filtered_normalized_integrated_clustered_prepMarker_majorCT.RDS")

DefaultAssay(object = cervix.integrated) <- "SCT"

cervix.integrated$major_celltypes[which(cervix.integrated$major_celltypes%in% c("CD14_Mono","CD16_Mono"))] <- "Monocyte"

Idents(cervix.integrated) <- "major_celltypes"

###################################
# Load CollecTRI network
##################################

net <- get_collectri(organism='human', split_complexes=FALSE)

mat <- as.matrix(cervix.integrated@assays[["SCT"]]@data)

net_summarize <-  net %>% filter(target %in% rownames(mat)) %>% group_by(source) %>%
  summarise(
    positive_targets = if(any(mor > 0))
      paste(sort(unique(target[mor > 0])), collapse = "///")
    else NA_character_,
    negative_targets = if(any(mor < 0))
      paste(sort(unique(target[mor < 0])), collapse = "///")
    else NA_character_
  ) %>%
  ungroup()

#################################
# Run VIPER
################################

acts <- run_viper(mat=mat, net=net, .source='source', .target='target', nes = TRUE,
                .mor='mor', minsize = 5)

saveRDS(acts,
        file="/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Cervix/Cervix_major_TF_VIPER/Cervix_major_TF_VIPER.rds")
