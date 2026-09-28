# subtype DE

for (folder in c("ALL_INFvsUNINF/", "NEG_INFvsUNINF/", "POS_INFvsUNINF/")) dir.create(folder, recursive = TRUE, showWarnings = FALSE)

library(dplyr)
library(ggplot2)
library(Seurat)
library(scCustomize)
library(table1)
library(SingleR)
packageVersion("Seurat")
# Load data
cervix.CD8T.integrated <- readRDS("/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_CD8T_normalized_integrated_clustered_prepMarker_SingleR.RDS")

DefaultAssay(object = cervix.CD8T.integrated) <- "SCT"

# Remove genes that do not need to be included for DE analysis
genes.remove <- readRDS(file = "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Cervix/Gene_Remove_DE.rds")

cervix.CD8T.integrated <- cervix.CD8T.integrated[!rownames(cervix.CD8T.integrated) %in% genes.remove,]

cervix.CD8T.integrated

table(cervix.CD8T.integrated$SingleR.pruned, cervix.CD8T.integrated$condition)

table(cervix.CD8T.integrated$SingleR.pruned, cervix.CD8T.integrated$ID)

table(cervix.CD8T.integrated$ID, cervix.CD8T.integrated$condition)
# Perform Wilcoxon Rank Sum test
DE.folder <- "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Cervix/Cervix_CD8T_Analysis/Cervix_CD8T_DE_INFvsUNINF"

Idents(cervix.CD8T.integrated) <- "SingleR.pruned"

for (ct in c("Central memory CD8 T cells", "Effector memory CD8 T cells", "Terminal effector CD8 T cells")) {
  Wilcox.DE <- FindMarkers(cervix.CD8T.integrated,
                           ident.1 = "infected",
                           ident.2 = "uninfected",
                           group.by = "condition",
                           subset.ident = ct,
                           recorrect_umi = FALSE,
                           logfc.threshold = 1,
                           assay = "SCT",
                           slot = "data",
                           min.pct = 0.1)

  Wilcox.DE.pos <- Wilcox.DE %>% filter(avg_log2FC >0 & p_val_adj < 0.05 ) %>%
    arrange(-avg_log2FC) %>% dplyr::select(-c("pct.1","pct.2"))
  Wilcox.DE.neg <- Wilcox.DE %>% filter(avg_log2FC <= 0 & p_val_adj < 0.05) %>%
    arrange(avg_log2FC) %>% dplyr::select(-c("pct.1","pct.2"))

  write.csv(Wilcox.DE,
            paste0("ALL_INFvsUNINF/",ct,"_INFvsUNINF_All_Wilcox.csv"))

  write.csv(Wilcox.DE.pos,
            paste0("POS_INFvsUNINF/",ct,"_Up_DE_Wilcox.csv"))

   write.csv(Wilcox.DE.neg,
            paste0("NEG_INFvsUNINF/",ct,"_Down_DE_Wilcox.csv"))
}

