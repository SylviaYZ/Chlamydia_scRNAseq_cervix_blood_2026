#  S02 tissue QC markers

library(dplyr)
library(ggplot2)
library(Seurat)
library(svglite)
library(writexl)
library(forcats)
library(patchwork)
library(ragg)

output_dir <- "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Figures/SuppFig1"

theme_jci <- theme_minimal(base_size = 8, base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 8),
    axis.title  = element_text(size = 8),
    plot.title  = element_blank(),
    legend.text  = element_text(size = 8),
    legend.title = element_text(size = 8)
  )

save_svg <- function(plot, file, width_in, height_in, bg = "white") {
  ggsave(
    filename = file,
    plot = plot,
    device  = svglite::svglite,
    width   = width_in, height = height_in, units = "in",
    bg = bg, limitsize = FALSE
  )
}

# =========================================================
# Shared distinct ID palette
# =========================================================
# 1-13 are your cervix IDs
# 14-24 are additional blood-only IDs
# ID 9 is preserved exactly as requested
id_cols_all <- c(
  "1"  = "#E41A1C",  # red
  "2"  = "#377EB8",  # blue
  "3"  = "#4DAF4A",  # green
  "4"  = "#984EA3",  # purple
  "5"  = "#FF7F00",  # orange
  "6"  = "#FFFF33",  # yellow
  "7"  = "#A65628",  # brown
  "8"  = "#F781BF",  # pink
  "9"  = "#17BECF",  # cyan
  "10" = "#BCBD22",  # olive
  "11" = "#1F77B4",  # dark blue
  "12" = "#2CA02C",  # dark green
  "13" = "#FF1493",  # deep pink
  "14" = "#8C564B",
  "15" = "#9467BD",
  "16" = "#D62728",
  "17" = "#7F7F7F",
  "18" = "#AEC7E8",
  "19" = "#98DF8A",
  "20" = "#FFBB78",
  "21" = "#C5B0D5",
  "22" = "#C49C94",
  "23" = "#F7B6D2",
  "24" = "#9EDAE5",
  "NA" = "grey75"
)

# =========================================================
# Cervix
# =========================================================
cervix.integrated <- readRDS(
  "/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_all_filtered_normalized_integrated_clustered_prepMarker_majorCT.RDS"
)

# -------------------------
# Mask patient ID
# -------------------------
set.seed(20251028)
ids <- cervix.integrated$ID
u_no_na <- unique(ids[!is.na(ids)])
k <- length(u_no_na)

mapping <- setNames(sample.int(k, size = k, replace = FALSE), u_no_na)
cervix.integrated$ID_mask <- as.integer(unname(mapping[as.character(ids)]))

# Key table
key_tbl <- cervix.integrated@meta.data %>%
  dplyr::select(ID, ID_mask) %>%
  distinct() %>%
  arrange(ID_mask)

# Stable factor for metadata
cervix.integrated$ID_mask <- factor(cervix.integrated$ID_mask, levels = key_tbl$ID_mask)

# Explicit plotting factor including NA
cervix.integrated$ID_mask_plot <- as.character(cervix.integrated$ID_mask)
cervix.integrated$ID_mask_plot <- forcats::fct_na_value_to_level(cervix.integrated$ID_mask_plot)

levs_num <- sort(unique(as.integer(as.character(
  cervix.integrated$ID_mask_plot[cervix.integrated$ID_mask_plot != "NA"]
))))
levs <- c(as.character(levs_num), "NA")
cervix.integrated$ID_mask_plot <- factor(cervix.integrated$ID_mask_plot, levels = levs)

# -------------------------
# Other metadata
# -------------------------
DefaultAssay(cervix.integrated) <- "SCT"
cervix.integrated$major_celltypes_updated <- cervix.integrated$major_celltypes
cervix.integrated$major_celltypes_updated[
  which(cervix.integrated$major_celltypes %in% c("CD14_Mono", "CD16_Mono"))
] <- "Monocytes"

cervix.integrated$major_celltypes_updated <- factor(
  cervix.integrated$major_celltypes_updated,
  levels = c("B", "CD4_T", "CD8_T", "MAIT", "NK", "Monocytes"),
  labels = c("B", "CD4 T", "CD8 T", "MAIT/NK T", "NK", "Monocyte")
)

cervix.integrated$condition <- factor(
  cervix.integrated$condition,
  levels = c( "infected", "uninfected"),
  labels = c( "CT+", "CT-")
)

Idents(cervix.integrated) <- "major_celltypes_updated"

# -------------------------
# Cervix ID UMAP with distinct colors
# -------------------------
cervix_id_cols <- id_cols_all[c(levels(cervix.integrated$ID_mask_plot))]
cervix_id_cols <- cervix_id_cols[!is.na(cervix_id_cols)]

Cervix_p1 <- DimPlot(
  cervix.integrated,
  reduction = "umap",
  group.by  = "ID_mask_plot",
  pt.size   = 3,
  raster    = TRUE,
  raster.dpi = c(600, 600),
  alpha     = 0.5
) +
  theme_jci +
  theme(
    legend.text  = element_text(size = 6),
    legend.title = element_text(size = 6)
  ) +
  scale_color_manual(values = cervix_id_cols, drop = TRUE) +
  guides(color = guide_legend(ncol = 2, byrow = TRUE, override.aes = list(size = 3)))

# Write cervix key with colors
key_tbl_out <- key_tbl %>%
  mutate(Color = unname(id_cols_all[as.character(ID_mask)]))

writexl::write_xlsx(
  list("ID_mask_key" = key_tbl_out),
  path = file.path(output_dir, "Cervix_ID_Mask_Convert.xlsx")
)

Cervix_p2 <- DimPlot(
  cervix.integrated, reduction = "umap", group.by = "orig.ident",
  pt.size = 3, raster = TRUE, raster.dpi = c(600, 600), alpha = 0.5
) +
  theme_jci +
  theme(
    legend.text  = element_text(size = 6),
    legend.title = element_text(size = 6, face = "plain")
  )

Cervix_p3 <- DimPlot(
  cervix.integrated, reduction = "umap", group.by = "condition",
  pt.size = 3, raster = TRUE, raster.dpi = c(600, 600), alpha = 0.5
) +
  theme_jci +
  theme(
    legend.text  = element_text(size = 6),
    legend.title = element_text(size = 6, face = "plain")
  )

Cervix_p4 <- DimPlot(
  cervix.integrated, reduction = "umap", group.by = "integrated_snn_res.2",
  pt.size = 3, raster = TRUE, raster.dpi = c(600, 600), alpha = 0.5,
  label = TRUE, label.size = 2
) +
  theme_jci +
  theme(
    legend.text  = element_text(size = 6),
    legend.title = element_text(size = 6, face = "plain")
  ) +
  guides(color = guide_legend(ncol = 3, byrow = TRUE, override.aes = list(size = 3)))

save_svg(
  Cervix_p1,
  file.path(output_dir, "SuppFig1A_UMAP_ID.svg"),
  width_in = 3, height_in = 2
)

save_svg(
  Cervix_p2,
  file.path(output_dir, "SuppFig1B_UMAP_Original_Ident.svg"),
  width_in = 3, height_in = 2
)

save_svg(
  Cervix_p3,
  file.path(output_dir, "SuppFig1C_UMAP_Original_condition.svg"),
  width_in = 3, height_in = 2
)

save_svg(
  Cervix_p4,
  file.path(output_dir, "SuppFig1D_UMAP_Original_cluster.svg"),
  width_in = 4, height_in = 2
)

features.cervix <- unique(c(
  "CD79A", "CD79B", "MS4A1", "CD19", "IGHD", "IGHM", "IGLC1", "IGKC",
  "CD3D", "CD3E", "CD3G", "CD27", "TRAC", "TRBC1", "TRDC", "TRGC1", "TRGC2", "CD4",
  "CD8A", "CD8B", "SLC4A10", "TRAV1-2", "KLRB1", "GNLY", "NKG7",
  "KLRD1", "KLRG1", "CD14", "FCGR3A", "FCGR3B",
  "HLA-DRA", "LYZ", "CD68", "CD86", "S100A8"
))

p.list <- FeaturePlot(
  cervix.integrated,
  features = features.cervix,
  combine  = FALSE,
  pt.size  = 0.005,
  order    = TRUE
)

p.list <- lapply(p.list, function(p) {
  p +
    theme_void(base_size = 7) +
    coord_fixed(1, clip = "off") +
    theme(
      plot.title      = element_text(size = 7, face = "bold", hjust = 0.5, margin = margin(b = 1)),
      plot.margin     = margin(1, 1, 1, 1),
      panel.border    = element_blank(),
      legend.position = "none"
    )
})

Cervix_p5_small <- wrap_plots(p.list, ncol = 10) +
  plot_annotation(theme = theme(plot.margin = margin(2, 2, 2, 2)))

svglite::svglite(file.path(output_dir, "SuppFig1E_features.svg"),
                 width = 6, height = 2.5)
print(Cervix_p5_small)
dev.off()

ragg::agg_png(file.path(output_dir, "SuppFig1E_features.png"),
              width = 6, height = 2.5, units = "in", res = 600)
print(Cervix_p5_small)
dev.off()

ggplot2::ggsave(
  filename = file.path(output_dir, "SuppFig1E_features_alt.svg"),
  plot     = Cervix_p5_small,
  width    = 6, height = 2.5, units = "in",
  device   = svglite::svglite
)

# =========================================================
# Blood
# =========================================================
blood.integrated <- readRDS(
  "/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/PBMC/SCTransformV2_Oct26/blood_all_filtered_normalized_integrated_clustered_prepMarker_majorCT.RDS"
)

# Start from cervix key so shared IDs keep same mask
set.seed(2026)

map_cervix <- setNames(key_tbl$ID_mask, as.character(key_tbl$ID))
next_start <- if (length(map_cervix)) max(as.integer(map_cervix), na.rm = TRUE) + 1L else 1L

blood_ids <- blood.integrated$ID
u_blood   <- unique(blood_ids)
u_blood   <- u_blood[!is.na(u_blood)]
new_only  <- setdiff(as.character(u_blood), names(map_cervix))

k_new   <- length(new_only)
pool    <- seq.int(from = next_start, length.out = k_new)

map_new <- if (k_new > 0L) {
  setNames(sample(pool, size = k_new, replace = FALSE), new_only)
} else {
  setNames(integer(0), character(0))
}

map_all <- c(map_cervix, map_new)

blood.integrated$ID_mask <- unname(map_all[as.character(blood_ids)])
blood.integrated$ID_mask <- as.integer(blood.integrated$ID_mask)

# sanity checks
overlap_ids <- intersect(as.character(u_blood), names(map_cervix))
if (length(overlap_ids)) {
  stopifnot(all(map_all[overlap_ids] == map_cervix[overlap_ids]))
}
unmapped <- unique(blood_ids[!is.na(blood_ids) & is.na(blood.integrated$ID_mask)])
if (length(unmapped)) {
  stop(sprintf("Some non-missing blood IDs did not get a mask: %s",
               paste(unmapped, collapse = ", ")))
}

combined_key <- tibble(
  ID = names(map_all),
  ID_mask = unname(map_all)
) %>%
  arrange(ID_mask) %>%
  mutate(Color = unname(id_cols_all[as.character(ID_mask)]))

writexl::write_xlsx(
  list("ID_mask_key" = combined_key),
  path = file.path(output_dir, "Combined_ID_Mask_Convert.xlsx")
)

# -------------------------
# Blood metadata
# -------------------------
table(blood.integrated$ID, blood.integrated$condition)

blood.integrated$condition <- factor(
  blood.integrated$condition,
  levels = c( "infected", "uninfected"),
  labels = c( "CT+", "CT-")
)

blood.integrated$major_celltypes <- factor(
  blood.integrated$major_celltypes,
  levels = c("CD4_T", "CD8_T", "NK_T", "gd_T", "Undecided"),
  labels = c("CD4 T", "CD8 T", "MAIT/NK T", "gd T", "Undecided")
)

blood.integrated <- blood.integrated[, blood.integrated$major_celltypes != "Undecided"]

# -------------------------
# Blood ID plotting factor
# participant 9 reuses cervix color
# 14-24 use other distinct colors
# -------------------------
blood.integrated$ID_mask_plot <- as.character(blood.integrated$ID_mask)
blood.integrated$ID_mask_plot <- forcats::fct_na_value_to_level(blood.integrated$ID_mask_plot)

blood_levs_num <- sort(unique(as.integer(as.character(
  blood.integrated$ID_mask_plot[blood.integrated$ID_mask_plot != "NA"]
))))
blood_levs <- c(as.character(blood_levs_num), "NA")
blood.integrated$ID_mask_plot <- factor(blood.integrated$ID_mask_plot, levels = blood_levs)

blood_id_cols <- id_cols_all[levels(blood.integrated$ID_mask_plot)]
blood_id_cols <- blood_id_cols[!is.na(blood_id_cols)]

# If blood has any IDs above 24, stop so you notice
missing_blood_cols <- setdiff(levels(blood.integrated$ID_mask_plot), names(id_cols_all))
missing_blood_cols <- setdiff(missing_blood_cols, "NA")
if (length(missing_blood_cols) > 0) {
  stop(sprintf(
    "No manual colors defined for blood ID_mask values: %s",
    paste(missing_blood_cols, collapse = ", ")
  ))
}

Blood_p1 <- DimPlot(
  blood.integrated,
  reduction = "umap",
  group.by  = "ID_mask_plot",
  pt.size   = 3,
  raster    = TRUE,
  raster.dpi = c(600, 600),
  alpha     = 0.5
) +
  theme_jci +
  theme(
    legend.text  = element_text(size = 6),
    legend.title = element_text(size = 6)
  ) +
  scale_color_manual(values = blood_id_cols, drop = TRUE) +
  guides(color = guide_legend(ncol = 2, byrow = TRUE, override.aes = list(size = 3)))

Blood_p2 <- DimPlot(
  blood.integrated, reduction = "umap", group.by = "orig.ident",
  pt.size = 3, raster = TRUE, raster.dpi = c(600, 600), alpha = 0.5
) +
  theme_jci +
  theme(
    legend.text  = element_text(size = 6),
    legend.title = element_text(size = 6, face = "plain")
  )

Blood_p3 <- DimPlot(
  blood.integrated, reduction = "umap", group.by = "condition",
  pt.size = 3, raster = TRUE, raster.dpi = c(600, 600), alpha = 0.5
) +
  theme_jci +
  theme(
    legend.text  = element_text(size = 6),
    legend.title = element_text(size = 6, face = "plain")
  )

Blood_p4 <- DimPlot(
  blood.integrated, reduction = "umap", group.by = "integrated_snn_res.2",
  pt.size = 3, raster = TRUE, raster.dpi = c(600, 600), alpha = 0.5,
  label = TRUE, label.size = 2
) +
  theme_jci +
  theme(
    legend.text  = element_text(size = 6),
    legend.title = element_text(size = 6, face = "plain")
  ) +
  guides(color = guide_legend(ncol = 4, byrow = TRUE, override.aes = list(size = 3)))

save_svg(
  Blood_p1,
  file.path(output_dir, "SuppFig2A_UMAP_ID.svg"),
  width_in = 3, height_in = 2
)

save_svg(
  Blood_p2,
  file.path(output_dir, "SuppFig2B_UMAP_Original_Ident.svg"),
  width_in = 3, height_in = 2
)

save_svg(
  Blood_p3,
  file.path(output_dir, "SuppFig2C_UMAP_Original_condition.svg"),
  width_in = 3, height_in = 2
)

save_svg(
  Blood_p4,
  file.path(output_dir, "SuppFig2D_UMAP_Original_cluster.svg"),
  width_in = 4, height_in = 2
)

features.blood <- unique(c(
  "CD3D", "CD3E", "CD3G", "CD27", "TRAC", "TRBC1", "CD4",
  "CD8A", "CD8B", "SLC4A10", "TRAV1-2", "KLRB1", "GNLY", "NKG7",
  "KLRD1", "KLRG1", "TRDC", "TRGC1", "TRGC2"
))

DefaultAssay(object = blood.integrated) <- "SCT"

p.list2 <- FeaturePlot(
  blood.integrated,
  features = features.blood,
  combine  = FALSE,
  pt.size  = 0.005,
  order    = TRUE
)

p.list2 <- lapply(p.list2, function(p) {
  p +
    theme_void(base_size = 7) +
    coord_fixed(1, clip = "off") +
    theme(
      plot.title      = element_text(size = 7, face = "bold", hjust = 0.5, margin = margin(b = 1)),
      plot.margin     = margin(1, 1, 1, 1),
      panel.border    = element_blank(),
      legend.position = "none"
    )
})

Blood_p5_small <- wrap_plots(p.list2, ncol = 10) +
  plot_annotation(theme = theme(plot.margin = margin(2, 2, 2, 2)))

svglite::svglite(file.path(output_dir, "SuppFig2E_features.svg"),
                 width = 6, height = 2.5)
print(Blood_p5_small)
dev.off()

ragg::agg_png(file.path(output_dir, "SuppFig2E_features.png"),
              width = 6, height = 2.5, units = "in", res = 600)
print(Blood_p5_small)
dev.off()

ggplot2::ggsave(
  filename = file.path(output_dir, "SuppFig2E_features_alt.svg"),
  plot     = Blood_p5_small,
  width    = 6, height = 2.5, units = "in",
  device   = svglite::svglite
)
