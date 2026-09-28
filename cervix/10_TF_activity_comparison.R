# TF activity comparison

library(dplyr)
library(ggplot2)
library(Seurat)
library(table1)
library(ggrepel)
library(tibble)
library(tidyr)
library(patchwork)
library(pheatmap)
library(EnhancedVolcano)
library(xml2)
library(purrr)
library(viper)
library(Matrix)

cervix.integrated <- readRDS("/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_all_filtered_normalized_integrated_clustered_prepMarker_majorCT.RDS")

DefaultAssay(object = cervix.integrated) <- "SCT"

cervix.integrated$major_celltypes[which(cervix.integrated$major_celltypes%in% c("CD14_Mono","CD16_Mono"))] <- "Monocyte"

Idents(cervix.integrated) <- "major_celltypes"

acts <- readRDS(file="/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Cervix/Cervix_major_TF_VIPER/Cervix_major_TF_VIPER.rds")

cols <- colnames(cervix.integrated)

acts_long <- acts %>%
  filter(condition %in% cols) %>%                # keep only cells present in Seurat object
  select(source, condition, score) %>%           # TF, cell barcode, activity score
  distinct()                                     # guard against duplicates

genes <- sort(unique(acts_long$source))
cell_index <- match(acts_long$condition, cols)
gene_index <- match(acts_long$source, genes)

VIPER_mat <- sparseMatrix(
  i = gene_index, j = cell_index, x = acts_long$score,
  dims = c(length(genes), length(cols)),
  dimnames = list(genes, cols)
)

setequal(unique(acts_long$condition), colnames(cervix.integrated))

VIPER_assay <- CreateAssayObject(counts = VIPER_mat)
VIPER_assay <- SetAssayData(VIPER_assay, layer = "data", new.data = VIPER_mat)
VIPER_assay <- SetAssayData(VIPER_assay, layer = "scale.data", new.data = as.matrix(VIPER_mat))

cervix.integrated[["VIPER"]] <- VIPER_assay
DefaultAssay(cervix.integrated) <- "VIPER"

cts <- sort(unique(cervix.integrated$major_celltypes))
markers_by_ct <- lapply(cts, function(ct) {
  obj <- subset(cervix.integrated, subset = major_celltypes == ct)
  DefaultAssay(obj) <- "VIPER"
  Idents(obj) <- "condition"
  FindMarkers(
    obj,
    ident.1 = "infected",
    ident.2 = "uninfected",
    min.pct = 0,
    logfc.threshold = 0,
    assay = "VIPER",
    slot = "scale.data",
    only.pos = TRUE,
    mean.fxn = Matrix::rowMeans,
    fc.name = "avg_diff"
  ) %>%
    rownames_to_column("regulator") %>%
    mutate(celltype = ct)
}) %>% bind_rows()

saveRDS(markers_by_ct, "markers_viper_by_celltype.rds")

CD4 <- markers_by_ct %>% filter(celltype == "CD4_T")

quantile(CD4$pct.1)

quantile(CD4$avg_diff)

write.csv(CD4,
          "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Cervix/Cervix_major_TF_VIPER/TF_VIPER_Seurat/VIPER_Seurat_CD4T.csv")

write.csv(CD4 %>% filter(pct.1 > 0.5, avg_diff > quantile(CD4$avg_diff, 0.9)),
          "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Cervix/Cervix_major_TF_VIPER/TF_VIPER_Seurat/VIPER_Seurat_CD4T_selected.csv")

CD8 <- markers_by_ct %>% filter(celltype == "CD8_T")

quantile(CD8$pct.1)

quantile(CD8$avg_diff)

write.csv(CD8,
          "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Cervix/Cervix_major_TF_VIPER/TF_VIPER_Seurat/VIPER_Seurat_CD8T.csv")

write.csv(CD8 %>% filter(pct.1 > 0.5, avg_diff > quantile(CD8$avg_diff, 0.9)),
          "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Cervix/Cervix_major_TF_VIPER/TF_VIPER_Seurat/VIPER_Seurat_CD8T_selected.csv")
