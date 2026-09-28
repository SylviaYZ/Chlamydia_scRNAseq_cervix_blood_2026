# CD4 DE trajectory

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
library(grid)
library(svglite)

packageVersion("Seurat")

# ============================================================
# Output directory
# ============================================================

output_dir <- "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Figures/Fig5/"

# ============================================================
# Final CD4 T-cell color palette
# ============================================================

# Finalized colors used across figures

celltype_colors <- c(
  "TN"     = "#8DD3C7",
  "Tfh"    = "#FFFFB3",
  "Th1"    = "#FB8072",
  "Th1/17" = "#80B1D3",
  "Th17"   = "#FDB462",
  "Th2"    = "#B3DE69",
  "TTE"    = "#FCCDE5",
  "Treg"   = "#BEBADA"
)

# Aggregated palette used when Th1/Th1-17/Th17/Th2 are
# combined into the meta-Th population.
#
# Meta-Th remains gray.

color_palette_aggregated <- c(
  "TN"   = celltype_colors["TN"],
  "Tfh"  = celltype_colors["Tfh"],
  "Th"   = "grey70",
  "TTE"  = celltype_colors["TTE"],
  "Treg" = celltype_colors["Treg"]
)

# ============================================================
# Differential expression volcano + Reactome network
# ============================================================

DE.volcano.plot <- function(
    ct,
    title,
    gene.bold,
    net.plot = FALSE,
    net.prop.gene.show = 0.5
) {

  # ----------------------------------------------------------
  # Constants
  # ----------------------------------------------------------

  base_pt  <- 6
  title_pt <- 8
  label_pt <- 6
  pt2gg    <- 1 / ggplot2::.pt

  # ----------------------------------------------------------
  # Load DE results
  # ----------------------------------------------------------

  DE <- read.csv(
    paste0(
      "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Cervix/Cervix_CD4T_analysis/Cervix_CD4T_DE_INFvsUNINF/ALL_INFvsUNINF/",
      ct,
      "_INFvsUNINF_All_Wilcox.csv"
    )
  )

  DE <- DE %>%
    dplyr::filter(
      p_val_adj < 0.05
    ) %>%
    dplyr::mutate(
      fc_label_keep = (
        avg_log2FC >= 2 |
          avg_log2FC <= -2
      ),
      label.print = dplyr::case_when(
        fc_label_keep &
          avg_log2FC >
          stats::quantile(
            avg_log2FC,
            0.99,
            na.rm = TRUE
          ) ~ X,

        fc_label_keep &
          avg_log2FC <
          stats::quantile(
            avg_log2FC,
            0.01,
            na.rm = TRUE
          ) ~ X,

        fc_label_keep &
          p_val_adj <
          stats::quantile(
            p_val_adj,
            0.02,
            na.rm = TRUE
          ) ~ X,

        TRUE ~ ""
      )
    )

  # ----------------------------------------------------------
  # Volcano colors
  # ----------------------------------------------------------

  keyvals <- ifelse(
    DE$avg_log2FC < 0 &
      DE$p_val_adj < 0.05,
    "royalblue",
    ifelse(
      DE$avg_log2FC > 0 &
        DE$p_val_adj < 0.05,
      "red",
      "black"
    )
  )

  names(keyvals)[keyvals == "red"] <-
    "Up-regulated"

  names(keyvals)[keyvals == "royalblue"] <-
    "Down-regulated"

  # ----------------------------------------------------------
  # Volcano plot
  # ----------------------------------------------------------

  DEplot <- EnhancedVolcano(
    DE,
    lab = DE$X,
    selectLab = DE$label.print,
    x = "avg_log2FC",
    y = "p_val_adj",
    pCutoff = 0.05,
    FCcutoff = 1,
    colCustom = keyvals,
    title = title,
    subtitle = "",
    caption = "",
    legendLabels = c("", ""),
    drawConnectors = TRUE,
    widthConnectors = 0.2,
    maxoverlapsConnectors = 30,
    pointSize = 0.45,
    labSize = 2,
    titleLabSize = title_pt,
    subtitleLabSize = base_pt,
    captionLabSize = base_pt,
    axisLabSize = base_pt,
    colAlpha = 0.7,
    labFace = "plain"
  ) +
    theme_minimal(
      base_size = base_pt
    ) +
    theme(
      legend.position = "none",

      plot.title = element_text(
        hjust = 0.5,
        size = title_pt
      ),

      axis.title = element_text(
        size = base_pt
      ),

      axis.text = element_text(
        size = base_pt
      ),

      plot.margin = unit(
        c(1, 1, 1, 1),
        "mm"
      ),

      axis.line = element_line(
        color = "black",
        linewidth = 0.25
      ),

      axis.ticks = element_line(
        color = "black",
        linewidth = 0.25
      ),

      axis.ticks.length = unit(
        0.8,
        "mm"
      ),

      panel.grid.major = element_line(
        color = "grey85",
        linewidth = 0.2
      ),

      panel.grid.minor = element_blank()
    ) +
    labs(
      y = "-Log10 adjusted P",
      x = "Log2 fold change"
    )

  ggsave(
    paste0(
      output_dir,
      ct,
      "_DE.svg"
    ),
    DEplot,
    device = svglite::svglite,
    width = 2,
    height = 2,
    units = "in",
    bg = "white",
    limitsize = FALSE
  )

  # ==========================================================
  # Reactome network
  # ==========================================================

  if (net.plot) {

    # --------------------------------------------------------
    # Reactome enrichment
    # --------------------------------------------------------

    DE.Entrez.SYMBOL <- bitr(
      DE$X,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    ) %>%
      dplyr::mutate(
        Mapping = "SYMBOL",
        USE = 1 * !duplicated(SYMBOL)
      )

    enrichedP <- enrichPathway(
      gene = DE.Entrez.SYMBOL$ENTREZID,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      readable = TRUE
    )

    FC <- setNames(
      DE$avg_log2FC,
      DE$X
    )

    show_cat <- enrich_category[[ct]]

    res <- enrichedP@result

    if (is.numeric(show_cat)) {

      res_sel <- dplyr::arrange(
        res,
        p.adjust
      ) |>
        dplyr::slice_head(
          n = show_cat
        )

    } else if (is.character(show_cat)) {

      res_sel <- dplyr::filter(
        res,
        Description %in% show_cat
      )

    } else {

      res_sel <- res
    }

    cat_names <- res_sel$Description

    genes_by_term <- setNames(
      lapply(
        strsplit(
          res_sel$geneID,
          "/"
        ),
        function(x) {
          unique(
            trimws(x)
          )
        }
      ),
      res_sel$Description
    )

    # --------------------------------------------------------
    # Top proportion of |FC| genes per pathway
    # --------------------------------------------------------

    label_genes <- unique(
      unlist(
        lapply(
          names(genes_by_term),
          function(cat) {

            genes <- genes_by_term[[cat]]

            vec <- FC[
              intersect(
                names(FC),
                genes
              )
            ]

            vec <- vec[
              !is.na(vec)
            ]

            if (!length(vec)) {
              return(
                character(0)
              )
            }

            abs_vec <- abs(vec)

            thr <- as.numeric(
              stats::quantile(
                abs_vec,
                net.prop.gene.show,
                na.rm = TRUE
              )
            )

            sel <- names(abs_vec)[
              abs_vec >= thr
            ]

            if (!length(sel)) {

              sel <- names(
                sort(
                  abs_vec,
                  decreasing = TRUE
                )
              )[1]

            }

            sel
          }
        )
      )
    )

    # --------------------------------------------------------
    # CNET
    # --------------------------------------------------------

    cnet <- cnetplot(
      enrichedP,
      foldChange = FC,
      node_label = "none",
      showCategory = show_cat,
      circular = FALSE,
      layout = "kk"
    ) +

      scale_color_gradient2(
        name = "log2FC",
        low = "blue",
        mid = "white",
        high = "red",
        midpoint = 0
      ) +

      ggraph::scale_edge_width(
        range = c(
          0.02,
          0.06
        )
      ) +

      ggraph::scale_edge_alpha(
        range = c(
          0.18,
          0.35
        ),
        guide = "none"
      ) +

      theme_void(
        base_size = 5
      ) +

      theme(
        legend.position = "bottom",

        legend.title = element_text(
          size = 5
        ),

        legend.text = element_text(
          size = 5
        ),

        legend.key.height = unit(
          1.6,
          "mm"
        ),

        legend.key.width = unit(
          3.6,
          "mm"
        ),

        plot.margin = unit(
          c(1, 1, 1, 1),
          "mm"
        ),

        plot.title = element_text(
          hjust = 0.5,
          size = 7
        )
      ) +

      guides(
        color = guide_colorbar(
          title.position = "top",
          barheight = unit(
            2,
            "mm"
          ),
          barwidth = unit(
            12,
            "mm"
          )
        )
      ) +

      coord_cartesian(
        clip = "off"
      )

    # --------------------------------------------------------
    # Remove existing size/radius/color scales
    # --------------------------------------------------------

    rm_idx_size <- which(
      vapply(
        cnet$scales$scales,
        function(s) {
          any(
            c(
              "size",
              "radius"
            ) %in% s$aesthetics
          )
        },
        logical(1)
      )
    )

    rm_idx_col <- which(
      vapply(
        cnet$scales$scales,
        function(s) {
          any(
            c(
              "colour",
              "color"
            ) %in% s$aesthetics
          )
        },
        logical(1)
      )
    )

    rm_idx <- unique(
      c(
        rm_idx_size,
        rm_idx_col
      )
    )

    if (length(rm_idx)) {

      cnet$scales$scales <-
        cnet$scales$scales[
          -rm_idx
        ]
    }

    # --------------------------------------------------------
    # Hide built-in node layers
    # --------------------------------------------------------

    node_layers <- which(
      sapply(
        cnet$layers,
        function(l) {
          inherits(
            l$geom,
            "GeomPoint"
          ) ||
            inherits(
              l$geom,
              "GeomNodePoint"
            )
        }
      )
    )

    if (length(node_layers)) {

      for (i in node_layers) {

        cnet$layers[[i]]$aes_params$alpha <- 0
        cnet$layers[[i]]$aes_params$size <- 0
      }
    }

    # --------------------------------------------------------
    # Node coordinates
    # --------------------------------------------------------

    node_df <- cnet$data

    cat_df <- subset(
      node_df,
      name %in% cat_names
    )

    gene_all_df <- subset(
      node_df,
      !(name %in% cat_names)
    )

    gene_all_df$fc <- FC[
      gene_all_df$name
    ]

    gene_lab_df <- subset(
      gene_all_df,
      name %in% label_genes
    )

    # Absolute node sizes

    cat_df$..sz <- 2
    gene_all_df$..sz <- 1

    # --------------------------------------------------------
    # Redraw nodes
    # --------------------------------------------------------

    cnet <- cnet +

      ggplot2::geom_point(
        data = gene_all_df,
        aes(
          x = x,
          y = y,
          color = fc,
          size = ..sz
        ),
        inherit.aes = FALSE,
        alpha = 0.95,
        show.legend = FALSE
      ) +

      ggplot2::geom_point(
        data = cat_df,
        aes(
          x = x,
          y = y,
          size = ..sz
        ),
        inherit.aes = FALSE,
        shape = 21,
        fill = "grey95",
        color = "grey40",
        stroke = 0.45,
        show.legend = FALSE
      ) +

      scale_size_identity() +

      scale_color_gradient2(
        name = "log2FC",
        low = "blue",
        mid = "white",
        high = "red",
        midpoint = 0,
        na.value = "grey80"
      ) +

      guides(
        color = guide_colorbar(
          title.position = "top",
          barheight = unit(
            2,
            "mm"
          ),
          barwidth = unit(
            12,
            "mm"
          )
        )
      ) +

      theme(
        legend.position = "bottom",

        legend.title = element_text(
          size = 5
        ),

        legend.text = element_text(
          size = 5
        )
      )

    # --------------------------------------------------------
    # Labels
    # --------------------------------------------------------

    cat_size_pt <- 4.2
    gene_size_pt <- 4.8

    cat_size_gg <- cat_size_pt * pt2gg
    gene_size_gg <- gene_size_pt * pt2gg

    cnet <- cnet +

      ggrepel::geom_text_repel(
        data = cat_df,
        aes(
          x = x,
          y = y,
          label = name
        ),
        inherit.aes = FALSE,
        size = cat_size_gg,
        segment.size = 0.03,
        max.overlaps = Inf,
        box.padding = 0.015,
        point.padding = 0.04,
        min.segment.length = 0,
        force = 0.02
      ) +

      ggrepel::geom_text_repel(
        data = gene_lab_df,
        aes(
          x = x,
          y = y,
          label = name
        ),
        inherit.aes = FALSE,
        size = gene_size_gg,
        segment.size = 0.03,
        max.overlaps = Inf,
        box.padding = 0.015,
        point.padding = 0.035,
        min.segment.length = 0,
        force = 0.015
      )

    ggsave(
      paste0(
        output_dir,
        ct,
        "_DE_net.svg"
      ),
      cnet,
      device = svglite::svglite,
      width = 1.5,
      height = 2,
      units = "in",
      bg = "white",
      limitsize = FALSE
    )
  }
}

# ============================================================
# Load CD4 object for DE analysis
# ============================================================

cervix.CD4T.integrated <- readRDS(
  "/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_CD4T_normalized_integrated_clustered_prepMarker_SingleR_manual.RDS"
)

Idents(cervix.CD4T.integrated) <-
  "SingleR.pruned.aggregated"

# ============================================================
# Reactome categories to display
# ============================================================

enrich_category <- list(

  "Naive CD4 T cells" = c(
    "Signaling by Interleukins",
    "Signaling by TGF-beta Receptor Complex",
    "FCERI mediated MAPK activation"
  ),

  "Follicular helper T cells" = c(
    "Signaling by Interleukins",
    "TNFR2 non-canonical NF-kB pathway",
    "Costimulation by the CD28 family"
  ),

  "T-helper cells" = c(
    "Signaling by Interleukins",
    "Interferon Signaling",
    "Signaling by TGF-beta Receptor Complex"
  ),

  "T regulatory cells" = c(
    "Signaling by Interleukins",
    "Interferon Signaling",
    "RUNX1 and FOXP3 control the development of regulatory T lymphocytes (Tregs)",
    "Signaling by TGF-beta Receptor Complex"
  )
)

# ============================================================
# Naive CD4 T cells
# ============================================================

# ============================================================
# Follicular helper T cells
# ============================================================

DE.volcano.plot(
  ct = "Follicular helper T cells",
  title = "CD4+ Tfh cells",
  gene.bold = NA,
  net.plot = TRUE
)

# ============================================================
# T-helper cells
# ============================================================

DE.volcano.plot(
  ct = "T-helper cells",
  title = "CD4+ Th cells",
  gene.bold = NA,
  net.plot = TRUE
)

# ============================================================
# T regulatory cells
# ============================================================

DE.volcano.plot(
  ct = "T regulatory cells",
  title = "CD4+ Treg cells",
  gene.bold = NA,
  net.plot = TRUE
)

# ============================================================
# Trajectory analysis heatmap
# ============================================================

sls.model.CD4 <- readRDS(
  "/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_CD4_slingshot_tradeSeq_model.RDS"
)

condRes <- readRDS(
  "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Cervix/Cervix_CD4T_analysis/Cervix_CD4T_trajectory/Cervix_CD4T_PCA_THelper/cervix_CD4lineage_condTest.RDS"
)

condRes <- condRes %>%
  filter(
    padj < 0.05
  ) %>%
  arrange(
    padj
  )

sum(
  condRes$padj <= 0.05,
  na.rm = TRUE
)

# ============================================================
# Recalculate smooth expression for significant genes
# ============================================================

yhatSmooth2 <- predictSmooth(
  sls.model.CD4,
  gene = rownames(condRes),
  nPoints = 50,
  tidy = FALSE
)

yhatSmooth2Scaled <- t(
  scale(
    t(yhatSmooth2)
  )
)

heatSmooth2_infected <- heatmap(
  yhatSmooth2Scaled[, 1:50],
  scale = "none",
  keep.dendro = TRUE,
  Rowv = TRUE,
  Colv = FALSE
)

# ============================================================
# Cluster trajectory genes
# ============================================================

heatSmooth2_infected_row.clusters <-
  as.hclust(
    heatSmooth2_infected$Rowv
  )

infected.DE.Traj.cluster <- cutree(
  heatSmooth2_infected_row.clusters,
  k = 5
)

infected.DE.Traj.cluster.df <- data.frame(
  gene = names(
    infected.DE.Traj.cluster
  ),
  cluster = infected.DE.Traj.cluster
)

infected.DE.Traj.cluster.df.output <- left_join(
  infected.DE.Traj.cluster.df,
  condRes %>%
    mutate(
      gene = rownames(condRes)
    ),
  by = "gene"
)

all.equal(
  names(
    infected.DE.Traj.cluster
  ),
  rownames(
    yhatSmooth2Scaled
  )
)

# ============================================================
# Infected trajectory heatmap
# ============================================================

infected_yhatSmooth2Scaled_long <-
  as.data.frame(
    yhatSmooth2Scaled[, 1:50]
  )

colnames(
  infected_yhatSmooth2Scaled_long
) <- 1:50

infected_yhatSmooth2Scaled_long$gene <-
  rownames(
    infected_yhatSmooth2Scaled_long
  )

infected_yhatSmooth2Scaled_long <-
  gather(
    infected_yhatSmooth2Scaled_long,
    key = "pseudotime",
    value = "scaled_predicted_expression",
    -gene
  ) %>%
  mutate(
    pseudotime = as.numeric(
      pseudotime
    )
  )

infected_yhatSmooth2Scaled_long <- left_join(
  infected_yhatSmooth2Scaled_long,
  infected.DE.Traj.cluster.df,
  by = "gene"
)

infected_yhatSmooth2Scaled_long$cluster <- factor(
  infected_yhatSmooth2Scaled_long$cluster,
  levels = c(
    2,
    5,
    3,
    1,
    4
  )
)

infected_heatmap_plot <- ggplot(
  infected_yhatSmooth2Scaled_long,
  aes(
    pseudotime,
    y = gene,
    fill = scaled_predicted_expression
  )
) +
  geom_tile() +

  scale_fill_gradient2(
    name = "Predicted\nExpression",
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(
      -3,
      7
    )
  ) +

  facet_grid(
    rows = vars(cluster),
    space = "free",
    scale = "free"
  ) +

  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),

    axis.title.x = element_text(
      size = 8
    ),

    axis.title.y = element_text(
      size = 8
    ),

    plot.title = element_text(
      size = 8
    )
  )

print(
  infected_heatmap_plot
)

# ============================================================
# Uninfected trajectory heatmap
# ============================================================

uninfected_yhatSmooth2Scaled_long <-
  as.data.frame(
    yhatSmooth2Scaled[, 51:100]
  )

colnames(
  uninfected_yhatSmooth2Scaled_long
) <- 1:50

uninfected_yhatSmooth2Scaled_long$gene <-
  rownames(
    uninfected_yhatSmooth2Scaled_long
  )

uninfected_yhatSmooth2Scaled_long <-
  gather(
    uninfected_yhatSmooth2Scaled_long,
    key = "pseudotime",
    value = "scaled_predicted_expression",
    -gene
  ) %>%
  mutate(
    pseudotime = as.numeric(
      pseudotime
    )
  )

uninfected_yhatSmooth2Scaled_long <- left_join(
  uninfected_yhatSmooth2Scaled_long,
  infected.DE.Traj.cluster.df,
  by = "gene"
)

uninfected_yhatSmooth2Scaled_long$cluster <- factor(
  uninfected_yhatSmooth2Scaled_long$cluster,
  levels = c(
    2,
    5,
    3,
    1,
    4
  )
)

uninfected_heatmap_plot <- ggplot(
  uninfected_yhatSmooth2Scaled_long,
  aes(
    pseudotime,
    y = gene,
    fill = scaled_predicted_expression
  )
) +
  geom_tile() +

  scale_fill_gradient2(
    name = "Predicted\nExpression",
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(
      -3,
      7
    )
  ) +

  facet_grid(
    rows = vars(cluster),
    space = "free",
    scale = "free"
  ) +

  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),

    axis.title.x = element_text(
      size = 8
    ),

    axis.title.y = element_text(
      size = 8
    ),

    plot.title = element_text(
      size = 8
    )
  )

print(
  uninfected_heatmap_plot
)

# ============================================================
# Heatmap legend formatting
# ============================================================

small_legend <- list(

  guides(
    fill = guide_colorbar(
      direction = "horizontal",
      title.position = "left",
      barheight = unit(
        3,
        "mm"
      ),
      barwidth = unit(
        17,
        "mm"
      ),
      ticks = FALSE
    )
  ),

  theme(
    legend.position = "bottom",

    legend.title = element_text(
      size = 7
    ),

    legend.text = element_text(
      size = 6
    ),

    legend.box.margin = margin(
      2,
      0,
      2,
      0
    ),

    legend.margin = margin(
      0,
      0,
      0,
      0
    ),

    legend.spacing.x = unit(
      2,
      "mm"
    ),

    legend.spacing.y = unit(
      2,
      "mm"
    )
  )
)

infected_heatmap_plot2 <-
  infected_heatmap_plot +
  small_legend

uninfected_heatmap_plot2 <-
  uninfected_heatmap_plot +
  small_legend

p_save <- ggarrange(

  uninfected_heatmap_plot2 +
    ggtitle(
      "  CT-"
    ),

  infected_heatmap_plot2 +
    ggtitle(
      "  CT+"
    ),

  ncol = 2,
  common.legend = TRUE,
  legend = "bottom"
)

p_save

ggsave(
  filename = file.path(
    output_dir,
    "trajectory_heatmap.tiff"
  ),
  plot = p_save,
  width = 8,
  height = 5,
  dpi = 600
)

ggsave(
  file.path(
    output_dir,
    "Fig5_trajectory_heatmap.svg"
  ),
  p_save,
  device = svglite::svglite,
  width = 3,
  height = 3,
  units = "in",
  bg = "white",
  limitsize = FALSE
)

# ============================================================
# Trajectory UMAP
# ============================================================

# IMPORTANT:
# Do NOT regenerate colors with ggplotColours().
#
# The trajectory plot uses the same finalized palette as the
# other CD4 figures.
#
# TN  = #8DD3C7
# Tfh = #FFFFB3
# Th  = grey70 (aggregated/meta-Th)
# TTE = #FCCDE5

color_palette1 <- c(
  "TN"  = unname(
    celltype_colors["TN"]
  ),

  "Tfh" = unname(
    celltype_colors["Tfh"]
  ),

  "Th"  = "grey70",

  "TTE" = unname(
    celltype_colors["TTE"]
  )
)

# ============================================================
# JCI-style theme
# ============================================================

font_family <- "Arial"

theme_jci <- theme_minimal(
  base_size = 8,
  base_family = font_family
) +
  theme(
    panel.grid = element_blank(),

    axis.line = element_line(
      color = "black"
    ),

    axis.ticks = element_line(
      color = "black"
    ),

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

    legend.text = element_text(
      size = 8
    ),

    legend.title = element_text(
      size = 8,
      face = "plain"
    )
  )

save_svg <- function(
    plot,
    file,
    width_in,
    height_in,
    bg = "white"
) {

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
# Reload Seurat object for trajectory
# ============================================================

cervix.CD4T.integrated <- readRDS(
  "/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_CD4T_normalized_integrated_clustered_prepMarker_SingleR_manual.RDS"
)

# Remove Tregs from trajectory analysis
cervix.CD4T.integrated <-
  cervix.CD4T.integrated[
    ,
    cervix.CD4T.integrated$SingleR.pruned.aggregated !=
      "T regulatory cells"
  ]

# ============================================================
# Convert to SingleCellExperiment
# ============================================================

cervix.CD4T.integrated.SCE <-
  as.SingleCellExperiment(
    cervix.CD4T.integrated,
    assay = "SCT"
  )

reducedDim(
  cervix.CD4T.integrated.SCE,
  "umap"
) <- Embeddings(
  cervix.CD4T.integrated,
  "umap"
)

reducedDim(
  cervix.CD4T.integrated.SCE,
  "pca"
) <- Embeddings(
  cervix.CD4T.integrated,
  "pca"
)

# ============================================================
# Slingshot trajectory
# ============================================================

cervix.CD4T.integrated.SCE <- slingshot(
  cervix.CD4T.integrated.SCE,
  reducedDim = "pca",
  clusterLabels =
    colData(
      cervix.CD4T.integrated.SCE
    )$SingleR.pruned.aggregated,
  start.clus = "Naive CD4 T cells"
)

# ============================================================
# Embed trajectory curves into UMAP
# ============================================================

cervix.CD4T.embed <- embedCurves(
  cervix.CD4T.integrated.SCE,
  reducedDim(
    cervix.CD4T.integrated.SCE,
    "umap"
  )
)

cervix.CD4T.curve <- slingCurves(
  cervix.CD4T.embed,
  as.df = TRUE
)

cervix.CD4T.mst <- slingMST(
  cervix.CD4T.embed,
  as.df = TRUE
) %>%
  mutate(
    shape = ifelse(
      Cluster == "Naive CD4 T cells",
      17,
      19
    )
  )

cervix.CD4T.mst %>%
  dplyr::select(
    Lineage,
    Order,
    Cluster
  )

# ============================================================
# Create UMAP dataframe
# ============================================================

df.CD4 <- cbind(
  reducedDims(
    cervix.CD4T.integrated.SCE
  )$umap,
  cervix.CD4T.integrated.SCE$SingleR.pruned.aggregated
) %>%
  as.data.frame()

colnames(
  df.CD4
)[3] <- "SingleR.pruned.aggregated"

df.CD4$umap_1 <- as.numeric(
  df.CD4$umap_1
)

df.CD4$umap_2 <- as.numeric(
  df.CD4$umap_2
)

# ============================================================
# Aggregate cell-type labels
# ============================================================

df.CD4$SingleR.pruned.aggregated <- factor(
  df.CD4$SingleR.pruned.aggregated,

  levels = c(
    "Naive CD4 T cells",
    "Follicular helper T cells",
    "T-helper cells",
    "Terminal effector CD4 T cells"
  ),

  labels = c(
    "TN",
    "Tfh",
    "Th",
    "TTE"
  )
)

# ============================================================
# UMAP axis stubs
# ============================================================

xr <- range(
  df.CD4$umap_1,
  na.rm = TRUE
)

yr <- range(
  df.CD4$umap_2,
  na.rm = TRUE
)

stubx <- diff(xr) * 0.06
stuby <- diff(yr) * 0.06

# ============================================================
# Base trajectory UMAP
# ============================================================

p.CD4 <- ggplot(
  df.CD4,
  aes(
    x = umap_1,
    y = umap_2
  )
) +

  geom_point(
    aes(
      color = SingleR.pruned.aggregated
    ),
    size = 0.001
  ) +

  theme_classic() +

  labs(
    color = "Cell type"
  ) +

  # ----------------------------------------------------------
# UPDATED FINAL CELL-TYPE COLORS
# ----------------------------------------------------------

scale_color_manual(
  values = color_palette1,
  breaks = c(
    "TN",
    "Tfh",
    "Th",
    "TTE"
  ),
  limits = c(
    "TN",
    "Tfh",
    "Th",
    "TTE"
  ),
  drop = FALSE
) +

  theme_jci +

  coord_fixed(
    1,
    clip = "off"
  ) +

  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),

    panel.grid = element_blank(),

    plot.margin = margin(
      2,
      2,
      2,
      2
    ),

    legend.position = "right",

    legend.background = element_rect(
      fill = "white",
      color = NA
    ),

    legend.key.size = unit(
      0.22,
      "cm"
    ),

    legend.text = element_text(
      size = 6
    )
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

# ============================================================
# Add trajectory curves
# ============================================================

plot_traj <- p.CD4 +

  geom_path(
    data = cervix.CD4T.curve %>%
      arrange(
        Order
      ),
    aes(
      group = Lineage
    )
  ) +

  guides(
    color = guide_legend(
      title = NULL,
      override.aes = list(
        size = 3
      )
    )
  ) +

  theme(
    text = element_text(
      size = 8
    ),

    legend.text = element_text(
      size = 8
    )
  )

plot_traj

# ============================================================
# Save trajectory UMAP
# ============================================================

save_svg(
  plot_traj,
  file.path(
    output_dir,
    "Fig6G_UMAP_trajectory.svg"
  ),
  width_in = 3,
  height_in = 3
)
