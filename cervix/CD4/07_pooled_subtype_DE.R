# pooled subtype DE
# Supplementary Table 9: exports differential expression results for pooled cervical CD4 subtypes.


library(dplyr)
library(ggplot2)
library(Seurat)
library(scCustomize)
library(table1)
library(SingleR)
packageVersion("Seurat")
# Load Data

# Run cervix/CD4/04_pool_helper_subtypes.R first; it generates this input file.
cervix.CD4T.integrated <- readRDS("/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_CD4T_normalized_integrated_clustered_prepMarker_SingleR_manual.RDS")

# Remove genes that do not need to be included for DE analysis
# Run cervix/07_differential_expression.R first; it generates this input file.
genes.remove <- readRDS(file = "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Cervix/Gene_Remove_DE.rds")

cervix.CD4T.integrated <- cervix.CD4T.integrated[!rownames(cervix.CD4T.integrated) %in% genes.remove,]

cervix.CD4T.integrated

table(cervix.CD4T.integrated$SingleR.pruned.aggregated, cervix.CD4T.integrated$condition)

table(cervix.CD4T.integrated$SingleR.pruned.aggregated, cervix.CD4T.integrated$ID)

table(cervix.CD4T.integrated$ID, cervix.CD4T.integrated$condition)
# Perform Wilcoxon Rank Sum test
DE.folder <- "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Cervix/Cervix_CD4T_analysis/Cervix_CD4T_DE_INFvsUNINF"
for (folder in c("ALL_INFvsUNINF/", "POS_INFvsUNINF/", "NEG_INFvsUNINF/")) {
  dir.create(file.path(DE.folder, folder), recursive = TRUE, showWarnings = FALSE)
}


Idents(cervix.CD4T.integrated) <- "SingleR.pruned.aggregated"

for (ct in c("Follicular helper T cells", "T-helper cells", "T regulatory cells")) {
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

  # Supplementary Table 9: pooled CD4 subtype DE; significant genes have p_val_adj < 0.05.
  write.csv(Wilcox.DE,
            paste0(DE.folder, "/", "ALL_INFvsUNINF/",ct,"_INFvsUNINF_All_Wilcox.csv"))

  write.csv(Wilcox.DE.pos,
            paste0(DE.folder, "/", "POS_INFvsUNINF/",ct,"_INFvsUNINF_Pos_Wilcox.csv"))

   write.csv(Wilcox.DE.neg,
            paste0(DE.folder, "/", "NEG_INFvsUNINF/",ct,"_INFvsUNINF_Neg_Wilcox.csv"))
}

