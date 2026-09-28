# helper subtype DE

for (folder in c("ALL_INFvsUNINF/", "NEG_INFvsUNINF/", "POS_INFvsUNINF/")) dir.create(folder, recursive = TRUE, showWarnings = FALSE)

library(dplyr)
library(ggplot2)
library(Seurat)
library(scCustomize)
library(table1)
library(SingleR)
packageVersion("Seurat")
# Load Data

cervix.CD4T.integrated <- readRDS("/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_CD4T_normalized_integrated_clustered_prepMarker_SingleR_manual.RDS")

cervix.CD4T.integrated <- cervix.CD4T.integrated[, cervix.CD4T.integrated$SingleR.pruned.aggregated == "T-helper cells"]

Idents(cervix.CD4T.integrated) <- "SingleR.pruned"

cervix.CD4T.integrated$SingleR.pruned[cervix.CD4T.integrated$SingleR.pruned =="Th1/Th17 cells"] <- "Th1_17 cells"
# Remove genes that do not need to be included for DE analysis
genes.remove <- readRDS(file = "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Cervix/Gene_Remove_DE.rds")

cervix.CD4T.integrated <- cervix.CD4T.integrated[!rownames(cervix.CD4T.integrated) %in% genes.remove,]

cervix.CD4T.integrated

table(cervix.CD4T.integrated$SingleR.pruned, cervix.CD4T.integrated$condition)

table(cervix.CD4T.integrated$SingleR.pruned, cervix.CD4T.integrated$ID)

# Perform Wilcoxon Rank Sum test
DE.folder <- "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Cervix/Cervix_CD4T_analysis/Cervix_CD4Thelper_DE_INFvsUNINF/"

Idents(cervix.CD4T.integrated) <- "SingleR.pruned"

for (ct in c("Th17 cells", "Th1_17 cells", "Th1 cells")) {
  Wilcox.DE <- FindMarkers(cervix.CD4T.integrated,
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
            paste0("ALL_INFvsUNINF/",ct,"_All_DE_Wilcox.csv"))

  write.csv(Wilcox.DE.pos,
            paste0("POS_INFvsUNINF/",ct,"_Up_DE_Wilcox.csv"))

   write.csv(Wilcox.DE.neg,
            paste0("NEG_INFvsUNINF/",ct,"_Down_DE_Wilcox.csv"))
}

