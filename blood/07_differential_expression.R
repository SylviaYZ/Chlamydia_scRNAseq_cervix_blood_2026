# differential expression

# Create output folders after updating the paths for your computer.
dir.create("/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Blood", recursive = TRUE, showWarnings = FALSE)


library(dplyr)
library(ggplot2)
library(Seurat)
library(scCustomize)
library(table1)
library(SingleR)
packageVersion("Seurat")
# Load data
# Run blood/05_annotate_major_cell_types.R first; it generates this input file.
blood.integrated <- readRDS("/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/PBMC/SCTransformV2_Oct26/blood_all_filtered_normalized_integrated_clustered_prepMarker_majorCT.RDS")

blood.integrated <- blood.integrated[, blood.integrated$major_celltypes != "Undecided"]

DefaultAssay(object = blood.integrated) <- "SCT"

Idents(blood.integrated) <- "major_celltypes"

blood.integrated

table(blood.integrated$major_celltypes, blood.integrated$condition)
# Remove MT and RB genes, and AC genes

# Pseudogenes
rownames(blood.integrated)[grepl("^A[CFLP][0-9.]+", rownames(blood.integrated))]

rownames(blood.integrated)[grepl("^Z[[:digit:]]", rownames(blood.integrated))]

rownames(blood.integrated)[grepl("^AUXG|^BX|^U[0-9.]+", rownames(blood.integrated))]

rownames(blood.integrated)[grepl("LINC|^LINC-", rownames(blood.integrated))]

rownames(blood.integrated)[grepl("^FP[0-9.]+|^CR[0-9.]+", rownames(blood.integrated))]

rownames(blood.integrated)[grepl("^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", rownames(blood.integrated))]

rownames(blood.integrated)[grepl("^MT-", rownames(blood.integrated))]

rownames(blood.integrated)[grepl("^CR[0-9.]+", rownames(blood.integrated))]

# Defined genes to remove
genes.remove <- rownames(blood.integrated)[grepl("^A[CFLP][0-9.]+|^Z[[:digit:]]|^AUXG|^BX|^U[0-9.]+|LINC|^LINC-|^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA|^MT-|^FP[0-9.]+|^CR[0-9.]+", rownames(blood.integrated))]

# Add back some genes begins with AP[digit][letter] and U[digit][letter], also some FP and CR
genes.remove <- genes.remove[!grepl("AP.*[A-Za-z]|U.*[A-Za-z]|FP.*[A-Za-z]|CR.*[A-Za-z]", genes.remove)]

genes.remove <- genes.remove[!genes.remove %in% c("CR1", "CR2")]
genes.remove

saveRDS(genes.remove, file = "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Blood/Gene_Remove_DE.rds")

# Remove long integrative non-coding genes
blood.integrated <- blood.integrated[!rownames(blood.integrated) %in% genes.remove, ]

blood.integrated

# Wilcoxon rank sum test of infected cells vs. uninfected cell in each celltype
DE.folder <- "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Blood/Blood_majorCT_DE_INFvsUNINF"
for (folder in c("All_INFvsUNINF/", "POS_INFvsUNINF/", "NEG_INFvsUNINF/")) {
  dir.create(file.path(DE.folder, folder), recursive = TRUE, showWarnings = FALSE)
}


Idents(blood.integrated) <- "major_celltypes"

for (ct in unique(blood.integrated$major_celltypes)) {
  Wilcox.DE <- FindMarkers(blood.integrated,
                           ident.1 = "infected",
                           ident.2 = "uninfected",
                           group.by = "condition",
                           subset.ident = ct,
                           recorrect_umi = FALSE,
                           logfc.threshold = 1,
                           assay = "SCT",
                           slot = "data",
                           min.pct = 0.1)

  Wilcox.DE.pos <- Wilcox.DE %>% filter(avg_log2FC >0 )
  Wilcox.DE.neg <- Wilcox.DE %>% filter(avg_log2FC <= 0 )

  write.csv(Wilcox.DE,
            paste0(DE.folder, "/", "All_INFvsUNINF/",ct,"_INFvsUNINF_All_Wilcox.csv"))

  write.csv(Wilcox.DE.pos,
            paste0(DE.folder, "/", "POS_INFvsUNINF/",ct,"_INFvsUNINF_Pos_Wilcox.csv"))

   write.csv(Wilcox.DE.neg,
            paste0(DE.folder, "/", "NEG_INFvsUNINF/",ct,"_INFvsUNINF_Neg_Wilcox.csv"))
}

