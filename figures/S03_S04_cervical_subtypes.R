#  S04 cervical subtypes

library(dplyr)
library(ggplot2)
library(Seurat)
library(table1)
library(pheatmap)
library(scCustomize)
library(RColorBrewer)
library(patchwork)
library(cowplot)
library(reshape2)
library(ggpubr)
library(grid)
library(svglite)
library(ragg)
library(readxl)
library(forcats)

packageVersion("Seurat")

output_dir <- "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Figures/SuppFig3"

font_family <- "Arial"

theme_jci <- theme_minimal(base_size = 8, base_family = font_family) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1, face = "plain"),
    axis.text.y = element_text(size = 8, face = "plain"),
    axis.title   = element_text(size = 8, face = "plain"),
    plot.title   = element_blank(),
    legend.text  = element_text(size = 8),
    legend.title = element_text(size = 8, face = "plain")
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

ID_mask <- read_xlsx(
  "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Figures/SuppFig3/Combined_ID_Mask_Convert.xlsx"
)

# helper to make explicit ID plotting factor
make_id_plot_factor <- function(x) {
  x_chr <- as.character(x)
  x_chr <- forcats::fct_na_value_to_level(x_chr, level = "NA")
  levs_num <- sort(unique(as.integer(as.character(x_chr[x_chr != "NA"]))))
  levs <- c(as.character(levs_num), "NA")
  factor(as.character(x_chr), levels = levs)
}

# =========================
# CD4 T cells
# =========================
cervix.CD4T.integrated <- readRDS(
  "/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_CD4T_normalized_integrated_clustered_prepMarker_SingleR_manual.RDS"
)

cervix.CD4T.integrated$ID_mask <- ID_mask$ID_mask[
  match(cervix.CD4T.integrated$ID, ID_mask$ID)
]

table(cervix.CD4T.integrated$ID_mask, cervix.CD4T.integrated$ID)

cervix.CD4T.integrated$ID_mask_plot <- make_id_plot_factor(cervix.CD4T.integrated$ID_mask)

cd4_id_cols <- id_cols_all[levels(cervix.CD4T.integrated$ID_mask_plot)]
cd4_id_cols <- cd4_id_cols[!is.na(cd4_id_cols)]

missing_cd4_cols <- setdiff(levels(cervix.CD4T.integrated$ID_mask_plot), names(id_cols_all))
missing_cd4_cols <- setdiff(missing_cd4_cols, "NA")
if (length(missing_cd4_cols) > 0) {
  stop(sprintf(
    "No manual colors defined for CD4 ID_mask values: %s",
    paste(missing_cd4_cols, collapse = ", ")
  ))
}

cervix.CD4T.integrated$condition <- factor(
  cervix.CD4T.integrated$condition,
  levels = c( "infected", "uninfected"),
  labels = c( "CT+", "CT-")
)

CD4_p1 <- DimPlot(
  cervix.CD4T.integrated,
  reduction = "umap",
  group.by = "ID_mask_plot",
  pt.size = 5,
  raster = TRUE,
  raster.dpi = c(600, 600),
  alpha = 0.8
) +
  theme_jci +
  theme(
    legend.text  = element_text(size = 6),
    legend.title = element_text(size = 6, face = "plain")
  ) +
  scale_color_manual(values = cd4_id_cols, drop = TRUE) +
  guides(
    color = guide_legend(ncol = 2, byrow = TRUE, override.aes = list(size = 3))
  )

CD4_p2 <- DimPlot(
  cervix.CD4T.integrated,
  reduction = "umap",
  group.by = "orig.ident",
  pt.size = 5,
  raster = TRUE,
  raster.dpi = c(600, 600),
  alpha = 0.8
) +
  theme_jci +
  theme(
    legend.text  = element_text(size = 6),
    legend.title = element_text(size = 6, face = "plain")
  )

CD4_p3 <- DimPlot(
  cervix.CD4T.integrated,
  reduction = "umap",
  group.by = "condition",
  pt.size = 5,
  raster = TRUE,
  raster.dpi = c(600, 600),
  alpha = 0.8
) +
  theme_jci +
  theme(
    legend.text  = element_text(size = 6),
    legend.title = element_text(size = 6, face = "plain")
  )

save_svg(
  CD4_p1,
  file.path(output_dir, "SuppFig3A_UMAP_ID.svg"),
  width_in = 3, height_in = 2
)

save_svg(
  CD4_p2,
  file.path(output_dir, "SuppFig3B_UMAP_Original_Ident.svg"),
  width_in = 3, height_in = 2
)

save_svg(
  CD4_p3,
  file.path(output_dir, "SuppFig3C_UMAP_Original_condition.svg"),
  width_in = 3, height_in = 2
)

# =========================
# CD8 T cells
# =========================
cervix.CD8T.integrated <- readRDS(
  "/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_CD8T_normalized_integrated_clustered_prepMarker_SingleR.RDS"
)

cervix.CD8T.integrated$ID_mask <- ID_mask$ID_mask[
  match(cervix.CD8T.integrated$ID, ID_mask$ID)
]

table(cervix.CD8T.integrated$ID_mask, cervix.CD8T.integrated$ID)

cervix.CD8T.integrated$ID_mask_plot <- make_id_plot_factor(cervix.CD8T.integrated$ID_mask)

cd8_id_cols <- id_cols_all[levels(cervix.CD8T.integrated$ID_mask_plot)]
cd8_id_cols <- cd8_id_cols[!is.na(cd8_id_cols)]

missing_cd8_cols <- setdiff(levels(cervix.CD8T.integrated$ID_mask_plot), names(id_cols_all))
missing_cd8_cols <- setdiff(missing_cd8_cols, "NA")
if (length(missing_cd8_cols) > 0) {
  stop(sprintf(
    "No manual colors defined for CD8 ID_mask values: %s",
    paste(missing_cd8_cols, collapse = ", ")
  ))
}

cervix.CD8T.integrated$condition <- factor(
  cervix.CD8T.integrated$condition,
  levels = c( "infected", "uninfected"),
  labels = c( "CT+", "CT-")
)

CD8_p1 <- DimPlot(
  cervix.CD8T.integrated,
  reduction = "umap",
  group.by = "ID_mask_plot",
  pt.size = 5,
  raster = TRUE,
  raster.dpi = c(600, 600),
  alpha = 0.8
) +
  theme_jci +
  theme(
    legend.text  = element_text(size = 6),
    legend.title = element_text(size = 6, face = "plain")
  ) +
  scale_color_manual(values = cd8_id_cols, drop = TRUE) +
  guides(
    color = guide_legend(ncol = 2, byrow = TRUE, override.aes = list(size = 3))
  )

CD8_p2 <- DimPlot(
  cervix.CD8T.integrated,
  reduction = "umap",
  group.by = "orig.ident",
  pt.size = 5,
  raster = TRUE,
  raster.dpi = c(600, 600),
  alpha = 0.8
) +
  theme_jci +
  theme(
    legend.text  = element_text(size = 6),
    legend.title = element_text(size = 6, face = "plain")
  )

CD8_p3 <- DimPlot(
  cervix.CD8T.integrated,
  reduction = "umap",
  group.by = "condition",
  pt.size = 5,
  raster = TRUE,
  raster.dpi = c(600, 600),
  alpha = 0.8
) +
  theme_jci +
  theme(
    legend.text  = element_text(size = 6),
    legend.title = element_text(size = 6, face = "plain")
  )

save_svg(
  CD8_p1,
  file.path(output_dir, "SuppFig4A_UMAP_ID.svg"),
  width_in = 3, height_in = 2
)

save_svg(
  CD8_p2,
  file.path(output_dir, "SuppFig4B_UMAP_Original_Ident.svg"),
  width_in = 3, height_in = 2
)

save_svg(
  CD8_p3,
  file.path(output_dir, "SuppFig4C_UMAP_Original_condition.svg"),
  width_in = 3, height_in = 2
)
