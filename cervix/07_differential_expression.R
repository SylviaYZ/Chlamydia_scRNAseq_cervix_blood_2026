# differential expression
# Supplementary Table 3: exports broad cervical cell-type differential expression results.

# Create output folders after updating the paths for your computer.
dir.create("/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Cervix", recursive = TRUE, showWarnings = FALSE)


library(dplyr)
library(ggplot2)
library(Seurat)
library(scCustomize)
library(table1)
library(SingleR)
packageVersion("Seurat")
# Load data

# Run cervix/05_annotate_major_cell_types.R first; it generates this input file.
cervix.integrated <- readRDS("/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_all_filtered_normalized_integrated_clustered_prepMarker_majorCT.RDS")

DefaultAssay(object = cervix.integrated) <- "SCT"

cervix.integrated$major_celltypes[which(cervix.integrated$major_celltypes%in% c("CD14_Mono","CD16_Mono"))] <- "Monocyte"

table(cervix.integrated$major_celltypes, cervix.integrated$condition)
# Remove MT and RB genes, and AC genes

# Pseudogenes
rownames(cervix.integrated)[grepl("^A[CFLP][0-9.]+", rownames(cervix.integrated))]

rownames(cervix.integrated)[grepl("^Z[[:digit:]]", rownames(cervix.integrated))]

rownames(cervix.integrated)[grepl("^AUXG|^BX|^U[0-9.]+", rownames(cervix.integrated))]

rownames(cervix.integrated)[grepl("LINC|^LINC-", rownames(cervix.integrated))]

rownames(cervix.integrated)[grepl("^FP[0-9.]+|^CR[0-9.]+", rownames(cervix.integrated))]

rownames(cervix.integrated)[grepl("^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", rownames(cervix.integrated))]

rownames(cervix.integrated)[grepl("^MT-", rownames(cervix.integrated))]

# Defined genes to remove
genes.remove <- rownames(cervix.integrated)[grepl("^A[CFLP][0-9.]+|^Z[[:digit:]]|^AUXG|^BX|^U[0-9.]+|LINC|^LINC-|^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA|^MT-|^FP[0-9.]+|^CR[0-9.]+", rownames(cervix.integrated))]
genes.remove

# Add back some genes begins with AP[digit][letter] and U[digit][letter], also some FP and CR
genes.remove <- genes.remove[!grepl("AP.*[A-Za-z]|U.*[A-Za-z]|FP.*[A-Za-z]|CR.*[A-Za-z]", genes.remove)]
genes.remove

saveRDS(genes.remove, file = "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Cervix/Gene_Remove_DE.rds")

# Remove long integrative non-coding genes
cervix.integrated <- cervix.integrated[!rownames(cervix.integrated) %in% genes.remove, ]

cervix.integrated

# Wilcoxon rank sum test of infected cells vs. uninfected cell in each celltype
DE.folder <- "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Cervix/Cervix_DE_INFvsUNINF"
for (folder in c("All_INFvsUNINF/", "POS_INFvsUNINF/", "NEG_INFvsUNINF/")) {
  dir.create(file.path(DE.folder, folder), recursive = TRUE, showWarnings = FALSE)
}


Idents(cervix.integrated) <- "major_celltypes"

for (ct in c("CD4_T", "CD8_T", "MAIT", "Monocyte", "B", "NK")) {
  Wilcox.DE <- FindMarkers(cervix.integrated,
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

  # Supplementary Table 3: broad cervical cell-type DE; significant genes have p_val_adj < 0.05.
  write.csv(Wilcox.DE,
            paste0(DE.folder, "/", "All_INFvsUNINF/",ct,"_INFvsUNINF_All_Wilcox.csv"))

  write.csv(Wilcox.DE.pos,
            paste0(DE.folder, "/", "POS_INFvsUNINF/",ct,"_INFvsUNINF_Pos_Wilcox.csv"))

   write.csv(Wilcox.DE.neg,
            paste0(DE.folder, "/", "NEG_INFvsUNINF/",ct,"_INFvsUNINF_Neg_Wilcox.csv"))
}

