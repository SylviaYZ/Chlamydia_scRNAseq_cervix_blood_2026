# CD4 subtypes

library(dplyr)
library(ggplot2)
library(Seurat)
library(scCustomize)
library(table1)
library(pheatmap)
library(RColorBrewer)
library(patchwork)
library(cowplot)
library(reshape2)
library(ggpubr)
library(grid)
library(tidyr)
library(svglite)

packageVersion("Seurat")

# ============================================================
# OUTPUT DIRECTORY
# ============================================================

output_dir <- "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Figures/Fig4"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# SHARED CELL-TYPE ORDERS
# ============================================================

celltype_order <- c(
  "TN",
  "Tfh",
  "Th1",
  "Th1/17",
  "Th17",
  "Th2",
  "TTE",
  "Treg"
)

celltype_levels_raw <- c(
  "Naive CD4 T cells",
  "Follicular helper T cells",
  "Th1 cells",
  "Th1/Th17 cells",
  "Th17 cells",
  "Th2 cells",
  "Terminal effector CD4 T cells",
  "T regulatory cells"
)


# Pooled CD4 populations
celltype_order_D <- c(
  "TN",
  "Tfh",
  "Treg",
  "Th",
  "TTE"
)

celltype_levels_raw_D <- c(
  "Naive CD4 T cells",
  "Follicular helper T cells",
  "T regulatory cells",
  "T-helper cells",
  "Terminal effector CD4 T cells"
)

# Helper T-cell subsets
celltype_order_E <- c(
  "Th1",
  "Th1/17",
  "Th17",
  "Th2"
)

celltype_levels_raw_E <- c(
  "Th1 cells",
  "Th1/Th17 cells",
  "Th17 cells",
  "Th2 cells"
)

condition_order <- c("CT-", "CT+")

# ============================================================
# GENERAL PLOTTING SETTINGS
# ============================================================

font_family <- "Arial"

theme_jci <- theme_minimal(
  base_size = 8,
  base_family = font_family
) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(
      size = 8,
      angle = 45,
      hjust = 1,
      vjust = 1,
      face = "plain"
    ),
    axis.text.y = element_text(
      size = 8,
      face = "plain"
    ),
    axis.title = element_text(
      size = 8,
      face = "plain"
    ),
    plot.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.title = element_text(
      size = 8,
      face = "plain"
    )
  )

save_svg <- function(plot, file, width_in, height_in, bg = "white") {
  ggsave(
    filename = file,
    plot = plot,
    device = svglite::svglite,
    width = width_in,
    height = height_in,
    units = "in",
    bg = bg,
    limitsize = FALSE
  )
}

# ============================================================
# CELL-TYPE COLORS
# ============================================================

celltype_colors <- c(
  "TN"     = "#8DD3C7",
  "Tfh"    = "#FFFFB3",
  "Treg"   = "#BEBADA",
  "Th1"    = "#FB8072",
  "Th1/17" = "#80B1D3",
  "Th17"   = "#FDB462",
  "Th2"    = "#B3DE69",
  "TTE"    = "#FCCDE5"
)

celltype_colors <- celltype_colors[celltype_order]

# Panel D: pooled Th remains gray
color_palette_D <- c(
  "TN"   = unname(celltype_colors["TN"]),
  "Tfh"  = unname(celltype_colors["Tfh"]),
  "Treg" = unname(celltype_colors["Treg"]),
  "Th"   = "grey70",
  "TTE"  = unname(celltype_colors["TTE"])
)

color_palette_D <- color_palette_D[celltype_order_D]

# Panel E
color_palette_E <- celltype_colors[celltype_order_E]

print(celltype_colors)
print(color_palette_D)
print(color_palette_E)

# ============================================================
# LOAD SEURAT OBJECT
# ============================================================

cervix.CD4T.integrated <- readRDS(
  "/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_CD4T_normalized_integrated_clustered_prepMarker_SingleR_manual.RDS"
)

cervix.CD4T.integrated$condition <- factor(
  cervix.CD4T.integrated$condition,
  levels = c("uninfected", "infected"),
  labels = condition_order
)

cervix.CD4T.integrated$SingleR.pruned <- factor(
  cervix.CD4T.integrated$SingleR.pruned,
  levels = celltype_levels_raw,
  labels = celltype_order
)

DefaultAssay(cervix.CD4T.integrated) <- "SCT"

# ============================================================
# FIGURE 4A — UMAP BY CELL TYPE
# ============================================================

um <- Embeddings(
  cervix.CD4T.integrated,
  "umap"
)

xr <- range(um[, 1], na.rm = TRUE)
yr <- range(um[, 2], na.rm = TRUE)

stubx <- diff(xr) * 0.06
stuby <- diff(yr) * 0.06

UMAP_4A <- DimPlot(
  cervix.CD4T.integrated,
  reduction = "umap",
  group.by = "SingleR.pruned",
  cols = celltype_colors,
  label = FALSE,
  pt.size = 3,
  raster = TRUE,
  raster.dpi = c(600, 600),
  shuffle = TRUE
) +
  scale_color_manual(
    values = celltype_colors,
    breaks = celltype_order,
    limits = celltype_order,
    drop = FALSE
  ) +
  guides(
    color = guide_legend(
      ncol = 1,
      reverse = FALSE,
      override.aes = list(size = 3, alpha = 1)
    )
  ) +
  theme_jci +
  coord_fixed(1, clip = "off") +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(2, 2, 2, 2),
    legend.position = "right",
    legend.background = element_rect(
      fill = "white",
      color = NA
    ),
    legend.key.size = unit(0.22, "cm"),
    legend.text = element_text(size = 6)
  ) +
  annotate(
    "segment",
    x = xr[1],
    xend = xr[1] + stubx,
    y = yr[1],
    yend = yr[1],
    linewidth = 0.3
  ) +
  annotate(
    "segment",
    x = xr[1],
    xend = xr[1],
    y = yr[1],
    yend = yr[1] + stuby,
    linewidth = 0.3
  ) +
  annotate(
    "text",
    x = xr[1] + stubx * 1.1,
    y = yr[1] - stuby * 0.25,
    label = "UMAP_1",
    size = 7 / ggplot2::.pt
  ) +
  annotate(
    "text",
    x = xr[1] - stubx * 0.25,
    y = yr[1] + stuby * 1.1,
    label = "UMAP_2",
    angle = 90,
    vjust = 0,
    size = 7 / ggplot2::.pt
  )

save_svg(
  UMAP_4A,
  file.path(output_dir, "Fig4A_UMAP_celltype.svg"),
  width_in = 3,
  height_in = 3
)

# ============================================================
# FIGURE 4C — UMAP SPLIT BY CONDITION
# ============================================================

UMAP_4C <- DimPlot(
  cervix.CD4T.integrated,
  reduction = "umap",
  group.by = "SingleR.pruned",
  split.by = "condition",
  cols = celltype_colors,
  label = FALSE,
  combine = TRUE,
  pt.size = 3,
  raster = TRUE,
  raster.dpi = c(600, 600),
  shuffle = TRUE
) +
  theme_jci +
  labs(
    x = NULL,
    y = NULL,
    title = NULL
  ) +
  scale_x_continuous(
    breaks = NULL,
    labels = NULL
  ) +
  scale_y_continuous(
    breaks = NULL,
    labels = NULL
  ) +
  coord_fixed(1, clip = "off") +
  NoLegend() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_text(size = 7)
  )

save_svg(
  UMAP_4C,
  file.path(output_dir, "Fig4C_UMAP_celltype_condition.svg"),
  width_in = 3,
  height_in = 1.5
)

# ============================================================
# FIGURE 4B — MARKER HEATMAP
# ============================================================

Idents(cervix.CD4T.integrated) <- "SingleR.pruned"

combined_averages <- AverageExpression(
  cervix.CD4T.integrated,
  return.seurat = TRUE
)

# Explicitly match heatmap groups to the shared UMAP/bar-plot palette.
heatmap_labels <- as.character(Idents(combined_averages))
if (anyNA(heatmap_labels) ||
    any(!heatmap_labels %in% names(celltype_colors))) {
  stop("Averaged cell-type labels do not match celltype_colors; check Idents(combined_averages).")
}

combined_averages$heatmap_celltype <- droplevels(factor(
  heatmap_labels,
  levels = celltype_order
))
Idents(combined_averages) <- "heatmap_celltype"

# Order colors by the actual grouping levels, including when a type is absent.
heatmap_group_colors <- unname(celltype_colors[
  levels(combined_averages$heatmap_celltype)
])
stopifnot(!anyNA(heatmap_group_colors))

heatmap <- Seurat::DoHeatmap(
  combined_averages,
  group.by = "heatmap_celltype",
  group.bar = TRUE,
  group.colors = heatmap_group_colors,
  features = c(
    "CCR7",
    "TCF7",
    "LEF1",
    "KLF2",
    "PLAC8",
    "ICAM2",
    "MYC",

    # Tfh
    "BCL6",
    "CXCR5",
    "CD69",
    "IL21",
    "IL10",
    "SLAMF6",

    "GZMH",
    "ANKRD28",
    "ANXA1",
    "CXCR3",
    "EOMES",
    "GZMK",
    "HLA-DPB1",
    "HLA-DMA",
    "ANXA2",
    "GZMA",
    "IFNG",
    "TBX21",
    "RUNX3",
    "CCR1",
    "CCR5",

    "NKG7",
    "GZMB",
    "PRF1",
    "GNLY",
    "LAG3",
    "HAVCR2",

    "PDCD1",
    "IL10",
    "RORC",
    "IL17F",
    "IL17A",
    "CCR6",
    "ADAM12",
    "IL16",
    "TNF",
    "IL23R",

    # Th2
    "CCR4",
    "CCR7",
    "GATA3",

    # Treg
    "FOXP3",
    "IL2RA",
    "IKZF2",
    "TNFRSF9",
    "IL1R2",
    "LAIR2",
    "IL1R1"
  ),
  label = TRUE,
  draw.lines = FALSE
) +
  scale_fill_gradientn(
    name = "Standardized\nExpression",
    colors = rev(
      RColorBrewer::brewer.pal(
        n = 4,
        name = "RdBu"
      )
    )
  ) +
  guides(color = "none") +
  theme(
    axis.text.y = element_text(size = 6),
    legend.title = element_text(size = 6),
    text = element_text(size = 6),
    legend.position = "bottom",
    legend.justification = "center",
    legend.direction = "horizontal"
  )

save_svg(
  heatmap,
  file.path(output_dir, "Fig4B_heatmap.svg"),
  width_in = 2.8,
  height_in = 6
)

# ============================================================
# FIGURE 4D — % OF TOTAL CD4
# ============================================================

Chi2_5celltypes <- read.csv(
  file = "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Cervix/Cervix_CD4T_analysis/Cervix_CD4T_enrichment/Cervix_INFvsUNINF_CD4T_Thelper_Chi2Enrich.csv"
)

inf_byCelltype_plot <- Chi2_5celltypes %>%
  dplyr::mutate(
    infected = InCluster_infected / inf_total,
    uninfected = InCluster_uninfected / uninf_total
  ) %>%
  dplyr::select(
    cell_type,
    infected,
    uninfected
  )

# Treg immediately after Tfh
inf_byCelltype_plot$cell_type <- factor(
  inf_byCelltype_plot$cell_type,
  levels = celltype_levels_raw_D,
  labels = celltype_order_D
)

inf_byCelltype_plot <- reshape2::melt(
  inf_byCelltype_plot,
  id.vars = "cell_type"
)

inf_byCelltype_plot$value <- as.numeric(
  inf_byCelltype_plot$value
)

inf_byCelltype_plot$posting <- inf_byCelltype_plot$value * 100

inf_byCelltype_plot$status <- factor(
  ifelse(
    inf_byCelltype_plot$variable == "infected",
    "CT+",
    "CT-"
  ),
  levels = condition_order
)

bar_df <- inf_byCelltype_plot %>%
  tidyr::complete(
    status,
    cell_type,
    fill = list(value = 0)
  ) %>%
  dplyr::mutate(
    status = factor(status, levels = condition_order),
    cell_type = factor(
      cell_type,
      levels = celltype_order_D
    )
  )

bar_4D1 <- ggplot(
  bar_df,
  aes(
    x = status,
    y = value,
    fill = cell_type
  )
) +
  # Default stacking matches legend order from top to bottom
  geom_col(
    position = position_fill(reverse = FALSE),
    width = 0.7,
    show.legend = FALSE
  ) +
  geom_point(
    alpha = 0,
    shape = 21,
    size = 4,
    stroke = 0,
    show.legend = TRUE
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_fill_manual(
    values = color_palette_D,
    breaks = celltype_order_D,
    limits = celltype_order_D,
    drop = FALSE,
    name = NULL
  ) +
  labs(
    x = NULL,
    y = "% of Total CD4"
  ) +
  coord_cartesian(
    ylim = c(0, 1),
    expand = FALSE
  ) +
  theme_jci +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 0.5
    ),
    legend.position = "right"
  ) +
  guides(
    fill = guide_legend(
      ncol = 1,
      reverse = FALSE,
      override.aes = list(
        alpha = 1,
        shape = 21,
        size = 4,
        color = NA,
        stroke = 0
      )
    )
  )

save_svg(
  bar_4D1,
  file.path(output_dir, "Fig4D1_bar_enrichment.svg"),
  width_in = 2.5,
  height_in = 2.5
)

# ============================================================
# FIGURE 4E — % OF CD4 Th
# ============================================================

Chi2_4celltypes <- read.csv(
  file = "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Cervix/Cervix_CD4T_analysis/Cervix_CD4T_enrichment/Cervix_INFvsUNINF_CD4T_Chi2Enrich.csv"
)

inf_byCelltype_plot2 <- Chi2_4celltypes %>%
  dplyr::mutate(
    infected = InCluster_infected / inf_total,
    uninfected = InCluster_uninfected / uninf_total
  ) %>%
  dplyr::select(
    cell_type,
    infected,
    uninfected
  )

inf_byCelltype_plot2$cell_type <- factor(
  inf_byCelltype_plot2$cell_type,
  levels = celltype_levels_raw_E,
  labels = celltype_order_E
)

inf_byCelltype_plot2 <- reshape2::melt(
  inf_byCelltype_plot2,
  id.vars = "cell_type"
)

inf_byCelltype_plot2$value <- as.numeric(
  inf_byCelltype_plot2$value
)

inf_byCelltype_plot2$posting <- inf_byCelltype_plot2$value * 100

inf_byCelltype_plot2$status <- factor(
  ifelse(
    inf_byCelltype_plot2$variable == "infected",
    "CT+",
    "CT-"
  ),
  levels = condition_order
)

bar_df2 <- inf_byCelltype_plot2 %>%
  tidyr::complete(
    status,
    cell_type,
    fill = list(value = 0)
  ) %>%
  dplyr::mutate(
    status = factor(status, levels = condition_order),
    cell_type = factor(
      cell_type,
      levels = celltype_order_E
    )
  )

bar_4D2 <- ggplot(
  bar_df2,
  aes(
    x = status,
    y = value,
    fill = cell_type
  )
) +
  geom_col(
    position = position_fill(reverse = FALSE),
    width = 0.7,
    show.legend = FALSE
  ) +
  geom_point(
    alpha = 0,
    shape = 21,
    size = 4,
    stroke = 0,
    show.legend = TRUE
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_fill_manual(
    values = color_palette_E,
    breaks = celltype_order_E,
    limits = celltype_order_E,
    drop = FALSE,
    name = NULL
  ) +
  labs(
    x = NULL,
    y = "% of CD4 Th"
  ) +
  coord_cartesian(
    ylim = c(0, 1),
    expand = FALSE
  ) +
  theme_jci +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 0.5
    ),
    legend.position = "right"
  ) +
  guides(
    fill = guide_legend(
      ncol = 1,
      reverse = FALSE,
      override.aes = list(
        alpha = 1,
        shape = 21,
        size = 4,
        color = NA,
        stroke = 0
      )
    )
  )

save_svg(
  bar_4D2,
  file.path(output_dir, "Fig4D2_bar_enrichment.svg"),
  width_in = 2.5,
  height_in = 2.5
)
