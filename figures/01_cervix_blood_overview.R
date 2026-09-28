# cervix blood overview

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
packageVersion("Seurat")

output_dir <- "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Figures/Fig2"

###############################################
# Cervix
##############################################

color_palette <- brewer.pal(n = 7, name = "Set1")

color_palette <- color_palette[-6]

font_family <- "Arial"

okabe_ito <- c("#000000","#E69F00","#56B4E9","#009E73",
               "#F0E442","#0072B2","#D55E00","#CC79A7")

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

cervix.integrated <- readRDS("/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_all_filtered_normalized_integrated_clustered_prepMarker_majorCT.RDS")

DefaultAssay(object = cervix.integrated) <- "SCT"

cervix.integrated$major_celltypes_updated <- cervix.integrated$major_celltypes

cervix.integrated$major_celltypes_updated[which(cervix.integrated$major_celltypes %in% c("CD14_Mono", "CD16_Mono"))] <- "Monocytes"

cervix.integrated$major_celltypes_updated <- factor(cervix.integrated$major_celltypes_updated,
                                                    levels = c("Monocytes" , "NK", "B", "MAIT","CD4_T", "CD8_T"),
                                                    labels = c("Monocyte", "NK","B", "MAIT/NK T","CD4 T", "CD8 T"))

# levels = c("B", "CD4_T", "CD8_T", "MAIT", "NK",
#            "Monocytes"),
# labels = c("B", "CD4 T", "CD8 T", "MAIT/NK T", "NK",
#            "Monocyte")

cervix.integrated$condition <- factor(cervix.integrated$condition,
                                      levels=c("uninfected", "infected"),
                                      labels=c("CT-","CT+"))

Idents(cervix.integrated) <- "major_celltypes_updated"

### UMAP with cell type colored

# Get UMAP ranges to place short axis stubs
um <- Embeddings(cervix.integrated, "umap")
xr <- range(um[, 1], na.rm = TRUE)
yr <- range(um[, 2], na.rm = TRUE)
stubx <- diff(xr) * 0.06
stuby <- diff(yr) * 0.06

UMAP_2A <- DimPlot(
  cervix.integrated, reduction = "umap",
  group.by = "major_celltypes_updated",
  cols = color_palette,
  label = FALSE,
  pt.size = 2,
  raster = TRUE,
  raster.dpi = c(600, 600),
  shuffle = TRUE
) +
  theme_jci +
  coord_fixed(1, clip = "off") +
  theme(
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.line  = element_blank(),

    panel.grid = element_blank(),
    plot.margin = margin(2, 2, 2, 2),
    legend.position = "right",
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key.size = unit(0.22, "cm"),
    legend.text = element_text(size = 6)
  ) +
  guides(colour = guide_legend(override.aes = list(size = 2), ncol = 1)) +
  annotate("segment", x = xr[1], xend = xr[1] + stubx, y = yr[1], yend = yr[1], linewidth = 0.3) +
  annotate("segment", x = xr[1], xend = xr[1], y = yr[1], yend = yr[1] + stuby, linewidth = 0.3) +
  annotate("text", x = xr[1] + stubx * 1.1, y = yr[1] - stuby * 0.25, label = "UMAP_1", size = 7/ggplot2::.pt) +
  annotate("text", x = xr[1] - stubx * 0.25, y = yr[1] + stuby * 1.1, label = "UMAP_2", angle = 90, vjust = 0, size = 7/ggplot2::.pt)

save_svg(
  UMAP_2A,
  file.path(output_dir, "Fig2A_UMAP_celltype.svg"),
  width_in = 3, height_in = 3
)

ggsave(paste0("/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Figures/Fig2/Fig2A_UMAP_celltype.tiff"),
       UMAP_2A,
       width = 2, height = 2, units = "in", dpi = 600)

UMAP_2C <- DimPlot(
  cervix.integrated,
  reduction = "umap",
  group.by = "major_celltypes_updated",
  split.by = "condition",
  label = FALSE,
  cols = color_palette,
  combine = TRUE,
  pt.size = 2.5,
  raster = TRUE,
  raster.dpi = c(600, 600),
  shuffle = TRUE
) +
  theme_jci +
  # remove all axis stuff and legend
  labs(x = NULL, y = NULL, title = NULL) +
  scale_x_continuous(breaks = NULL, labels = NULL) +
  scale_y_continuous(breaks = NULL, labels = NULL) +
  coord_fixed(1, clip = "off") +
  NoLegend() +
  theme(
    axis.title   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    axis.line    = element_blank(),
    panel.grid   = element_blank(),
    legend.position = "none",

    # keep only CT+/CT- facet titles (make the strip minimal)
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_text(size = 7)
  )

save_svg(
  UMAP_2C,
  file.path(output_dir, "Fig2C_UMAP_celltype_condition.svg"),
  width_in = 3, height_in = 1.5
)

ggsave(paste0("/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Figures/Fig2/Fig2C_UMAP_celltype_condition.tiff"),
       UMAP_2C,
       width = 6, height = 3, units = "in", dpi = 600)

### Violin cell type

violin.marker <- c("CD14", "FCGR3B", "LYZ", "S100A8", "KLRB1","GNLY", "NKG7",
                   "CD79A", "MS4A1", "CD3D",  "TRAC", "SLC4A10","CD4",
                   "CD8A",  "CD8B")

violin_2B <- VlnPlot(
  cervix.integrated,
  features = violin.marker,
  stack = TRUE,            # stack all genes in one panel
  flip  = TRUE,            # genes on Y (more readable in small height)
  pt.size = 0,             # no points to avoid clutter
  combine = TRUE
) +
  theme_jci +
  scale_y_continuous(breaks = NULL, labels = NULL) +
  coord_cartesian(clip = "off") +
  theme(
    legend.position = "none",
    axis.title   = element_blank(),
    axis.text.y  = element_text(angle = 45, size = 6),  # gene labels (6 pt ≈ JCI minimum)
    axis.text.x  = element_text(size = 6, angle = 45, hjust = 1, vjust = 1),
    axis.ticks   = element_blank(),
    axis.line    = element_blank(),
    panel.grid   = element_blank(),
    strip.text.y.right = element_text(angle = 0, size = 6),
    panel.spacing.y = unit(0, "pt"),
    plot.margin  = margin(2, 2, 2, 2)
  )

save_svg(
  violin_2B,
  file.path(output_dir, "Fig2B_violin_celltype.svg"),
  width_in =2, height_in = 3
)

### Enrichment plot

Chi2_12celltypes <- read.csv(file = "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Cervix/Cervix_enrichment_analysis/Cervix_major_celltypes/Cervix_INFvsUNINF_12celltypes_Chi2Enrich.csv")

inf_byCelltype_plot <- Chi2_12celltypes %>%
  mutate(infected = InCluster_infected/inf_total,
         uninfected = InCluster_uninfected / uninf_total) %>%
  dplyr::select(cell_type, infected, uninfected)

inf_byCelltype_plot$cell_type <- factor(inf_byCelltype_plot$cell_type,
                                        levels = c("Monocytes" , "NK", "B", "MAIT","CD4_T", "CD8_T"),
                                        labels = c("Monocyte", "NK","B", "MAIT/NK T","CD4 T", "CD8 T"))

inf_byCelltype_plot <- melt(inf_byCelltype_plot, id.vars = "cell_type")
inf_byCelltype_plot$value <- as.numeric(inf_byCelltype_plot$value)
inf_byCelltype_plot$posting <- as.numeric(inf_byCelltype_plot$value)*100

inf_byCelltype_plot$status <- factor(ifelse(inf_byCelltype_plot$variable == "infected", "CT+", "CT-"),
                                     levels= c("CT-", "CT+"))

heatmap_2C <- ggplot(inf_byCelltype_plot, aes(y = status, x = cell_type)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = "white", high = "red", limits = c(0, 1), labels = function(x) x * 100) +
  geom_text(aes(label = round(posting, digits = 1)), size = 3.5) +
  labs(fill = "Proportion")  +
  coord_cartesian(xlim = c(1, 6), ylim = c(1, 2), clip = "off") + theme_jci +
  theme(axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

bar_df <- inf_byCelltype_plot %>%
  mutate(status = factor(status, levels = c("CT-", "CT+"))) %>%
  tidyr::complete(status, cell_type, fill = list(value = 0))

bar_2C <- ggplot(bar_df, aes(x = status, y = value)) +
  # bars: do not create legend
  geom_col(aes(fill = cell_type), position = "fill", width = 0.7, show.legend = FALSE) +
  # legend-only points (drive a point-style legend; hidden on plot)
  geom_point(
    aes(fill = cell_type),
    alpha = 0,      # invisible on panel
    shape = 21,     # circle uses 'fill'
    size  = 6,      # adjust circle size here
    stroke = 0,     # no outline
    show.legend = TRUE
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = color_palette, drop = FALSE, name = NULL) +
  labs(x = NULL, y = "% of Cervical Leukocytes ") +
  coord_cartesian(ylim = c(0, 1), expand = FALSE) +
  theme_jci +
  theme(
    axis.text.x      = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
    legend.position  = "right",
    legend.key       = element_blank(),   # no key box behind the circle
    legend.background= element_blank(),
    legend.key.size   = unit(4, "mm"),   # size of each legend key
    legend.key.height = unit(4, "mm"),
    legend.key.width  = unit(4, "mm"),
    legend.spacing.y  = unit(0.5, "mm"), # vertical gap between keys
    legend.box.margin = margin(0, 0, 0, 0) # trim outer legend padding
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(
        alpha  = 1,   # visible in legend
        shape  = 21,
        size   = 4,   # legend circle size (match 'geom_point' size)
        color  = NA,  # <-- no border
        stroke = 0
      )
    )
  )

ggsave(
  filename = file.path(output_dir, "Fig2C_bar_enrichment.svg"),
  plot = bar_2C,
  device = svglite::svglite,
  width = 2.5, height = 2.5, units = "in",
  bg = "white", limitsize = FALSE
)

ggsave("/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Figures/Fig2/Fig2C_heatmap_enrichment.tiff",
       plot = heatmap_2C, dpi = 600,
       width = 4, height = 1.8, unit = "in")

ggsave("/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Figures/Fig2/Fig2C_bar_enrichment.tiff",
       plot = bar_2C, dpi = 600,
       width = 3, height =3, unit = "in")

###############################################
# blood
##############################################

blood.integrated <- readRDS("/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/PBMC/SCTransformV2_Oct26/blood_all_filtered_normalized_integrated_clustered_prepMarker_majorCT.RDS")
table(blood.integrated$ID, blood.integrated$condition)

blood.integrated$condition <- factor(blood.integrated$condition,
                                      levels=c("uninfected", "infected"),
                                      labels=c("CT-","CT+"))

blood.integrated$major_celltypes <- factor(blood.integrated$major_celltypes,
                                           levels = c("NK_T","CD4_T", "CD8_T", "gd_T", "Undecided"),
                                           labels = c("MAIT/NK T", "CD4 T", "CD8 T",  "gd T", "Undecided"))

blood.integrated <- blood.integrated[, blood.integrated$major_celltypes != "Undecided"]

cell_type_colors <- c("CD4 T" = "#FF7F00",  # Customize colors as needed
                      "CD8 T" = "#A65628",
                      "MAIT/NK T" = "#984EA3",
                      "gd T" = "#00CED1")

# UMAP

um_blood <- Embeddings(blood.integrated, "umap")
xr <- range(um_blood[, 1], na.rm = TRUE)
yr <- range(um_blood[, 2], na.rm = TRUE)
stubx <- diff(xr) * 0.06
stuby <- diff(yr) * 0.06

UMAP_2E <- DimPlot(
  blood.integrated, reduction = "umap",
  group.by = "major_celltypes",
  cols = cell_type_colors,
  label = FALSE,
  pt.size = 2,
  raster = TRUE,
  raster.dpi = c(600, 600),
  shuffle = TRUE
) +
  theme_jci +
  coord_fixed(1, clip = "off") +
  theme(
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.line  = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(2, 2, 2, 2),
    legend.position = "right",
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key.size = unit(0.22, "cm"),
    legend.text = element_text(size = 6)
  ) +
  guides(colour = guide_legend(override.aes = list(size = 2), ncol = 1)) +
  annotate("segment", x = xr[1], xend = xr[1] + stubx, y = yr[1], yend = yr[1], linewidth = 0.3) +
  annotate("segment", x = xr[1], xend = xr[1], y = yr[1], yend = yr[1] + stuby, linewidth = 0.3) +
  annotate("text", x = xr[1] + stubx * 1.1, y = yr[1] - stuby * 0.25,
           label = "UMAP_1", size = 7/ggplot2::.pt) +
  annotate("text", x = xr[1] - stubx * 0.25, y = yr[1] + stuby * 1.1,
           label = "UMAP_2", angle = 90, vjust = 0, size = 7/ggplot2::.pt)

save_svg(
  UMAP_2E,
  file.path(output_dir, "Fig2E_UMAP_celltype.svg"),
  width_in = 3, height_in = 3
)

# +
#   theme(
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     axis.line = element_line(color = "black"),
#     axis.title = element_text(size = 8),
#     axis.text.x = element_text(size = 6),
#     axis.text.y = element_text(size = 6),
#     text = element_text(size = 8),
#     plot.margin = margin(2, 2, 2, 2),
#     # plot.title = element_blank(),
#     legend.position = "inside",
#     legend.position.inside = c(0.35, 0.8),
#     legend.justification = c("right", "bottom"),
#     legend.background = element_rect(fill = "white", color = NA),
#     legend.key.size = unit(0.3, "cm"),
#     legend.text = element_text(size = 8)
#   )

UMAP_2F <- DimPlot(
  blood.integrated,
  reduction = "umap",
  group.by ="major_celltypes",
  split.by = "condition",
  label = FALSE,
  cols = cell_type_colors,
  combine = TRUE,
  pt.size = 2,
  raster = TRUE,
  raster.dpi = c(600, 600),
  shuffle = TRUE
) +
  theme_jci +
  # remove all axis stuff and legend
  labs(x = NULL, y = NULL, title = NULL) +
  scale_x_continuous(breaks = NULL, labels = NULL) +
  scale_y_continuous(breaks = NULL, labels = NULL) +
  coord_fixed(1, clip = "off") +
  NoLegend() +
  theme(
    axis.title   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    axis.line    = element_blank(),
    panel.grid   = element_blank(),
    legend.position = "none",

    # keep only CT+/CT- facet titles (make the strip minimal)
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_text(size = 7)
  )

save_svg(
  UMAP_2F,
  file.path(output_dir, "Fig2F_UMAP_celltype_condition.svg"),
  width_in = 3, height_in = 1.5
)
#
# +
#   theme(
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     axis.line = element_line(color = "black"),
#     axis.title = element_text(size = 8),
#     axis.text.x = element_text(size = 6),
#     axis.text.y = element_text(size = 6),
#     text = element_text(size = 8),
#     plot.margin = margin(2, 2, 2, 2),
#     plot.title = element_blank(),
#     legend.position = "inside",
#     legend.position.inside = c(0.15, 0.8),
#     legend.justification = c("right", "bottom"),
#     legend.background = element_rect(fill = "white", color = NA),
#     legend.key.size = unit(0.3, "cm"),
#     legend.text = element_text(size = 8)
#   )

# Violin plot

violin.marker2 <- c("CD3D", "CD3E", "TRAC", "SLC4A10", "KLRB1","GNLY", "NKG7", "CD4", "CD8A", "CD8B",  "TRDC","TRGC1")

Idents(blood.integrated) <- "major_celltypes"

violin_2E <- VlnPlot(
  blood.integrated,
  features = violin.marker2,
  stack = TRUE,            # stack all genes in one panel
  flip  = TRUE,            # genes on Y (more readable in small height)
  pt.size = 0,             # no points to avoid clutter
  combine = TRUE
) +
  theme_jci +
  scale_y_continuous(breaks = NULL, labels = NULL) +
  coord_cartesian(clip = "off") +
  theme(
    legend.position = "none",
    axis.title   = element_blank(),
    axis.text.y  = element_text(angle = 45, size = 6),  # gene labels (6 pt ≈ JCI minimum)
    axis.text.x  = element_text(size = 6, angle = 45, hjust = 1, vjust = 1),
    axis.ticks   = element_blank(),
    axis.line    = element_blank(),
    panel.grid   = element_blank(),
    strip.text.y.right = element_text(angle = 0, size = 6),
    panel.spacing.y = unit(0, "pt"),
    plot.margin  = margin(2, 2, 2, 2)
  )

save_svg(
  violin_2E,
  file.path(output_dir, "Fig2E_violin_celltype.svg"),
  width_in =2, height_in = 2.5
)

  # theme(
  #   panel.background = element_blank(),
  #   panel.grid.major = element_blank(),
  #   panel.grid.minor = element_blank(),
  #   axis.line = element_line(color = "black"),
  #   axis.title = element_text(size = 2),      # Smaller axis title
  #   plot.title = element_blank(),
  #   text = element_text(size = 2),            # Smaller general text
  #   axis.text.x = element_text(size = 2, angle = 45, hjust = 1),
  #   axis.text.y = element_text(size = 2),
  #   plot.margin = margin(2, 2, 2, 2)
  # )

ggsave(paste0("/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Figures/Fig2/Fig2E_violin_celltype.tiff"),
       violin_2E,
       width = 4, height = 5, units = "in", dpi = 600)

# Enrichment plot

Chi2_12celltypes_blood <- read.csv(file = "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Blood/Blood_majorCT_enrichment_analysis/Blood_INFvsUNINF_4celltypes_Chi2Enrich.csv")

inf_byCelltype_plot_blood <- Chi2_12celltypes_blood %>%
  mutate(infected = InCluster_infected/inf_total,
         uninfected = InCluster_uninfected / uninf_total) %>%
  dplyr::select(cell_type, infected, uninfected)

inf_byCelltype_plot_blood$cell_type <- factor(inf_byCelltype_plot_blood$cell_type,
                                              levels = c("NK_T","CD4_T", "CD8_T", "gd_T", "Undecided"),
                                              labels = c("MAIT/NK T", "CD4 T", "CD8 T",  "gd T", "Undecided"))

inf_byCelltype_plot_blood <- melt(inf_byCelltype_plot_blood, id.vars = "cell_type")
inf_byCelltype_plot_blood$value <- as.numeric(inf_byCelltype_plot_blood$value)
inf_byCelltype_plot_blood$posting <- as.numeric(inf_byCelltype_plot_blood$value)*100

inf_byCelltype_plot_blood$status <- factor(ifelse(inf_byCelltype_plot_blood$variable == "infected", "CT+", "CT-"),
                                     levels= c("CT-", "CT+"))

heatmap_2F <- ggplot(inf_byCelltype_plot_blood, aes(y = status, x = cell_type)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = "white", high = "red", limits = c(0, 1)) +
  geom_text(aes(label = round(posting, digits = 1)), size = 3.5) +
  labs(fill = "Proportion") + theme_jci +
  theme(axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

bar_df_blood <- inf_byCelltype_plot_blood %>%
  mutate(status = factor(status, levels = c("CT-", "CT+"))) %>%
  tidyr::complete(status, cell_type, fill = list(value = 0))

bar_2F <- ggplot(bar_df_blood, aes(x = status, y = value)) +
  # bars: do NOT create legend (otherwise you'll get boxes)
  geom_col(aes(fill = cell_type),
           position = "fill", width = 0.7, show.legend = FALSE) +

  # legend-only points (hidden on plot, but define circle keys)
  geom_point(aes(fill = cell_type),
             alpha = 0,       # invisible on the panel
             shape = 21,      # filled circle
             size  = 4,       # legend circle size
             stroke = 0,      # no outline
             show.legend = TRUE) +

  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = cell_type_colors, drop = FALSE, name = NULL) +
  labs(x = NULL, y = "% of PBMC-derived T Cells") +
  coord_cartesian(ylim = c(0, 1), expand = FALSE) +
  theme_jci +
  theme(
    axis.text.x       = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
    legend.position   = "right",
    legend.key        = element_blank(),
    legend.background = element_blank(),
    legend.key.size   = unit(4, "mm"),
    legend.key.height = unit(4, "mm"),
    legend.key.width  = unit(4, "mm"),
    legend.spacing.y  = unit(0.5, "mm"),
    legend.box.margin = margin(0, 0, 0, 0)
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(
        alpha  = 1,   # visible in legend
        shape  = 21,
        size   = 4,
        color  = NA,  # no border
        stroke = 0
      )
    )
  )

ggsave(
  filename = file.path(output_dir, "Fig2F_bar_enrichment.svg"),
  plot = bar_2F,
  device = svglite::svglite,
  width = 2.5, height = 2.5, units = "in",
  bg = "white", limitsize = FALSE
)

  # theme(
  #   text = element_text(size = 6),
  #   axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1, face = "bold"),
  #   axis.text.y = element_text(size = 6, face = "bold"),
  #   axis.title = element_blank(), # Remove axis titles
  #   axis.line = element_blank(),  # Remove axis lines
  #   axis.ticks = element_blank(), # Remove axis ticks
  #   panel.background = element_blank(), # Remove panel background
  #   plot.title = element_blank() # Remove plot title
  # ) +
  # coord_cartesian(xlim = c(1, 4), ylim = c(1, 2), clip = "off")

ggsave("/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Figures/Fig2/Fig2F_heatmap_enrichment.tiff",
       plot = heatmap_2F, dpi = 600, unit = "in",
       width = 4, height = 1.8)

ggsave("/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Figures/Fig2/Fig2F_bar_enrichment.tiff",
       plot = bar_2F, dpi = 600,
       width = 3, height =3, unit = "in")

###############################################
# combine figures
##############################################

