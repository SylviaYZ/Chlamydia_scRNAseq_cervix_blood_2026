# CD8 subtypes

library(dplyr)
library(ggplot2)
library(Seurat)
library(table1)
library(ggrepel)
library(EnhancedVolcano)
library(ReactomePA)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(scales)
library(RColorBrewer)
library(patchwork)
library(cowplot)
library(slingshot)
library(tradeSeq)
library(ComplexHeatmap)
library(tidyr)
library(ggpubr)
library(reshape2)
packageVersion("Seurat")

output_dir <- "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Figures/Fig8"

default_cols <- hue_pal()(4)
color_palette  <- c(default_cols[1], "#39B600", "#D89000", default_cols[4])

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

cervix.CD8T.integrated <- readRDS("/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_CD8T_normalized_integrated_clustered_prepMarker_SingleR.RDS")

cervix.CD8T.integrated$condition <- factor(cervix.CD8T.integrated$condition,
                                           levels=c("uninfected", "infected"),
                                           labels=c("CT-","CT+"))

cervix.CD8T.integrated$SingleR.pruned <- factor(cervix.CD8T.integrated$SingleR.pruned,
                                                levels = c("Naive CD8 T cells",
                                                           "Central memory CD8 T cells",
                                                           "Effector memory CD8 T cells",
                                                           "Terminal effector CD8 T cells"),
                                                labels = c("TN",
                                                           "TCM",
                                                           "TEM",
                                                           "TTE"))

DefaultAssay(cervix.CD8T.integrated) <- "SCT"

um <- Embeddings(cervix.CD8T.integrated, "umap")
xr <- range(um[, 1], na.rm = TRUE)
yr <- range(um[, 2], na.rm = TRUE)
stubx <- diff(xr) * 0.06
stuby <- diff(yr) * 0.06

UMAP_8A <- DimPlot(
  cervix.CD8T.integrated, reduction = "umap", group.by = "SingleR.pruned",
  label = FALSE,
  pt.size = 3,
  raster = TRUE,
  raster.dpi = c(600, 600),
  shuffle = TRUE,
  cols = color_palette
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
    legend.background = element_rect(fill = "white", color = NA),
    legend.key.size = unit(0.22, "cm"),
    legend.text = element_text(size = 6)
  ) +
  annotate("segment", x = xr[1], xend = xr[1] + stubx, y = yr[1], yend = yr[1], linewidth = 0.3) +
  annotate("segment", x = xr[1], xend = xr[1], y = yr[1], yend = yr[1] + stuby, linewidth = 0.3) +
  annotate("text", x = xr[1] + stubx * 1.1, y = yr[1] - stuby * 0.25, label = "UMAP_1", size = 7/ggplot2::.pt) +
  annotate("text", x = xr[1] - stubx * 0.25, y = yr[1] + stuby * 1.1, label = "UMAP_2", angle = 90, vjust = 0, size = 7/ggplot2::.pt)

save_svg(
  UMAP_8A,
  file.path(output_dir, "Fig8A_UMAP_celltype.svg"),
  width_in = 3, height_in = 3
)

UMAP_8C <- DimPlot(
  cervix.CD8T.integrated, reduction = "umap", group.by = "SingleR.pruned",
  split.by = "condition",
  label = FALSE,
  combine = TRUE,
  pt.size = 3,
  raster = TRUE,
  raster.dpi = c(600, 600),
  shuffle = TRUE,
  cols = color_palette
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
  UMAP_8C,
  file.path(output_dir, "Fig8C_UMAP_celltype_condition.svg"),
  width_in = 3, height_in = 1.5
)

Idents(cervix.CD8T.integrated) <- "SingleR.pruned"

combined_averages <- AverageExpression(cervix.CD8T.integrated, return.seurat = TRUE)

heatmap <- DoHeatmap(combined_averages, features = c("CCR7", "TCF7", "LEF1", "KLF2", "MAL",
                                                     "MYC", "GZMK", "GZMM","GZMH",
                                                     "NKG7", "CD74", "CD38", "GZMA", "IFNG",
                                                     "ZNF683", "ITGAE", "ITGA1", "XCL1",
                                                     "XCL2", "JAML","PRF1", "GNLY",
                                                     "GZMB", "HAVCR2", "LAG3", "CTLA4", "TIGIT",
                                                     "ENTPD1", "TOX2", "CX3CR1", "FGFBP2",
                                                     "FCGR3A", "SPON2", "KLRG1"
                                                     ), label = TRUE , group.colors = color_palette,
                     draw.lines = FALSE)  +
  scale_fill_gradientn(name = "Standardized\nExpression",
                       colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu"))) +
  guides(color = "none")+
  theme(
    #axis.text.x = element_text(size = 18),  # Cell type labels (X-axis)
    axis.text.y = element_text(size = 8),  # Gene labels (Y-axis)
    legend.title = element_text(size = 6),
    text = element_text(size = 6),
    legend.position = "bottom",
    legend.justification = "center",
    legend.direction = "horizontal"
  )

ggsave(file.path(output_dir, "heatmap_8B.tiff"),
       heatmap,
       width = 3, height = 7.5, units = "in", dpi = 600)

save_svg(
  heatmap,
  file.path(output_dir, "Fig8_heatmap.svg"),
  width_in = 2.5, height_in = 7
)

# Heatmap for enrichment
Chi2_12celltypes <- read.csv(file = "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Cervix/Cervix_CD8T_Analysis/Cervix_CD8T_Enrichment/Cervix_INFvsUNINF_CD8T_SingleR_Chi2Enrich.csv")

inf_byCelltype_plot <- Chi2_12celltypes %>%
  mutate(infected = InCluster_infected/inf_total,
         uninfected = InCluster_uninfected / uninf_total) %>%
  dplyr::select(cell_type, infected, uninfected)

inf_byCelltype_plot$cell_type <- factor(inf_byCelltype_plot$cell_type,
                                        levels = c("Naive CD8 T cells",
                                                   "Central memory CD8 T cells",
                                                   "Effector memory CD8 T cells",
                                                   "Terminal effector CD8 T cells"),
                                        # labels = c("CD8 TN",
                                        #            "CD8 TCM",
                                        #            "CD8 TEM",
                                        #            "CD8 TTE"),
                                        labels = c("TN",
                                                   "TCM",
                                                   "TEM",
                                                   "TTE"))

inf_byCelltype_plot <- melt(inf_byCelltype_plot, id.vars = "cell_type")
inf_byCelltype_plot$value <- as.numeric(inf_byCelltype_plot$value)
inf_byCelltype_plot$posting <- as.numeric(inf_byCelltype_plot$value)*100

inf_byCelltype_plot$status <- factor(ifelse(inf_byCelltype_plot$variable == "infected", "CT+", "CT-"),
                                     levels= c("CT-", "CT+"))

bar_df <- inf_byCelltype_plot %>%
  mutate(status = factor(status, levels = c("CT-", "CT+"))) %>%
  tidyr::complete(status, cell_type, fill = list(value = 0))

heat1 <- ggplot(inf_byCelltype_plot, aes(y = status, x = cell_type)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = "white", high = "red", limits = c(0, 1)) +
  geom_text(aes(label = round(posting, digits = 1)), size = 5) +
  labs(fill = "Proportion") +
  theme(
    text = element_text(size = 8),
    #axis.text.x = element_text(size = 18, angle = 45, hjust = 1, vjust = 1,face = "bold"),
    axis.text.x = element_text(size = 8, face = "bold"),
    axis.text.y = element_text(size = 8, face = "bold"),
    axis.title = element_blank(), # Remove axis titles
    axis.line = element_blank(),  # Remove axis lines
    axis.ticks = element_blank(), # Remove axis ticks
    panel.background = element_blank(), # Remove panel background
    plot.title = element_blank() # Remove plot title
  ) +
  coord_cartesian(xlim = c(1, 4), ylim = c(1, 2), clip = "off")

heat1

ggsave(file.path(output_dir, "heatmap_8D.tiff"),
       heat1,
       width = 5, height = 2, units = "in", dpi = 600)

bar_8D <- ggplot(bar_df, aes(x = status, y = value, fill = cell_type)) +
  geom_col(aes(fill = cell_type),
           position = "fill", width = 0.7, show.legend = FALSE) +
  # legend-only points (hidden on plot, but define circle keys)
  geom_point(aes(fill = cell_type),
             alpha = 0,       # invisible on the panel
             shape = 21,      # filled circle
             size  = 4,       # legend circle size
             stroke = 0,      # no outline
             show.legend = TRUE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(
    values = color_palette,
    drop = FALSE,
    name = NULL
  ) +
  labs(x = NULL, y = "% of Total CD8") +
  coord_cartesian(ylim = c(0, 1), expand = FALSE) +
  theme_jci +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
    legend.position = "right"
  )+
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
  filename = file.path(output_dir, "Fig8D_bar_enrichment.svg"),
  plot = bar_8D,
  device = svglite::svglite,
  width = 2.5, height = 2.5, units = "in",
  bg = "white", limitsize = FALSE
)

