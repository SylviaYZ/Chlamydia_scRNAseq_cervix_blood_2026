# CD8 DE, pathways and trajectory

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

output_dir <- "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Figures/Fig9/"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

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

DE.volcano.plot <- function(ct, title, gene.bold,
                            net.plot = FALSE,
                            net.prop.gene.show = 0.5){

  # --- constants ---
  base_pt  <- 6
  title_pt <- 8
  label_pt <- 6
  pt2gg    <- 1/ggplot2::.pt

  # Run cervix/CD8/05_subtype_DE.R first; it generates this input file.
  DE <- read.csv(paste0("/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Cervix/Cervix_CD8T_Analysis/Cervix_CD8T_DE_INFvsUNINF/ALL_INFvsUNINF/",ct,"_INFvsUNINF_All_Wilcox.csv"))

  DE <- DE %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    dplyr::mutate(
      fc_label_keep = avg_log2FC >= 2 | avg_log2FC <= -2,
      label.print = dplyr::case_when(
        fc_label_keep & avg_log2FC > stats::quantile(avg_log2FC, 0.99, na.rm = TRUE) ~ X,
        fc_label_keep & avg_log2FC < stats::quantile(avg_log2FC, 0.01, na.rm = TRUE) ~ X,
        fc_label_keep & p_val_adj < stats::quantile(p_val_adj,   0.02, na.rm = TRUE) ~ X,
        TRUE ~ ""
      )
    )

  keyvals <- ifelse(
    DE$avg_log2FC < 0 & DE$p_val_adj < 0.05, 'royalblue',
    ifelse(DE$avg_log2FC > 0 & DE$p_val_adj < 0.05, 'red', 'black')
  )
  names(keyvals)[keyvals == 'red']       <- 'Up-regulated'
  names(keyvals)[keyvals == 'royalblue'] <- 'Down-regulated'

  # --- VOLCANO (2×2") ---
  DEplot <- EnhancedVolcano(
    DE,
    lab = DE$X,
    selectLab = DE$label.print,
    x = 'avg_log2FC', y = 'p_val_adj',
    pCutoff = 0.05, FCcutoff = 1,
    colCustom = keyvals,
    title = title, subtitle = "", caption = "",
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
    theme_minimal(base_size = base_pt) +
    theme(
      legend.position = "none",
      plot.title  = element_text(hjust = 0.5, size = title_pt),
      axis.title  = element_text(size = base_pt),
      axis.text   = element_text(size = base_pt),
      plot.margin = unit(c(1, 1, 1, 1), "mm"),
      axis.line  = element_line(color = "black", linewidth = 0.25),
      axis.ticks = element_line(color = "black", linewidth = 0.25),
      axis.ticks.length = unit(0.8, "mm"),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.2),
      panel.grid.minor = element_blank()
    ) +
    labs(y = "-Log10 adjusted P", x= "Log2 fold change")

  ggsave(
    paste0(output_dir, ct, "_DE.svg"),
    DEplot, device = svglite::svglite,
    width = 2, height = 2, units = "in",
    bg = "white", limitsize = FALSE
  )

  if (net.plot) {
    # --- Reactome enrichment ---
    DE.Entrez.SYMBOL <- bitr(
      DE$X, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db
    ) %>% dplyr::mutate(Mapping = "SYMBOL", USE = 1*!duplicated(SYMBOL))

    enrichedP <- enrichPathway(
      gene = DE.Entrez.SYMBOL$ENTREZID,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      readable = TRUE
    )

    FC <- setNames(DE$avg_log2FC, DE$X)
    show_cat <- enrich_category[[ct]]

    res <- enrichedP@result
    if (is.numeric(show_cat)) {
      res_sel <- dplyr::arrange(res, p.adjust) |> dplyr::slice_head(n = show_cat)
    } else if (is.character(show_cat)) {
      res_sel <- dplyr::filter(res, Description %in% show_cat)
    } else {
      res_sel <- res
    }
    cat_names <- res_sel$Description

    genes_by_term <- setNames(
      lapply(strsplit(res_sel$geneID, "/"), function(x) unique(trimws(x))),
      res_sel$Description
    )

    # top 25% |FC| per shown pathway (no cap)
    label_genes <- unique(unlist(lapply(names(genes_by_term), function(cat) {
      genes <- genes_by_term[[cat]]
      vec <- FC[intersect(names(FC), genes)]
      vec <- vec[!is.na(vec)]
      if (!length(vec)) return(character(0))
      abs_vec <- abs(vec)
      thr <- as.numeric(stats::quantile(abs_vec, net.prop.gene.show, na.rm = TRUE))
      sel <- names(abs_vec)[abs_vec >= thr]
      if (!length(sel)) sel <- names(sort(abs_vec, decreasing = TRUE))[1]
      sel
    })))

    # --- CNET (1.5" × 2"): tiny gene nodes, thin edges, stacked legend ---
    cnet <- cnetplot(
      enrichedP,
      foldChange   = FC,
      node_label   = "none",          # we'll label manually
      showCategory = show_cat,
      circular     = FALSE, layout = "kk"
    ) +
      # keep color legend for log2FC
      scale_color_gradient2(
        name = "log2FC",
        low = "blue", mid = "white", high = "red", midpoint = 0
      ) +
      # edge aesthetics (thin)
      ggraph::scale_edge_width(range = c(0.02, 0.06)) +
      ggraph::scale_edge_alpha(range = c(0.18, 0.35), guide = "none") +
      theme_void(base_size = 5) +                 # smaller base just for cnet
      theme(
        legend.position   = "bottom",             # stacked legends at bottom
        legend.title      = element_text(size = 5),
        legend.text       = element_text(size = 5),
        legend.key.height = unit(1.6, "mm"),
        legend.key.width  = unit(3.6, "mm"),
        plot.margin       = unit(c(1, 1, 1, 1), "mm"),
        plot.title        = element_text(hjust = 0.5, size = 7)
      ) +
      guides(
        color = guide_colorbar(title.position = "top",
                               barheight = unit(2, "mm"), barwidth = unit(12, "mm"))
      ) +
      coord_cartesian(clip = "off")

    # ---- remove any size/radius/colour scales so our sizes/colors aren't rescaled ----
    rm_idx_size <- which(vapply(cnet$scales$scales,
                                function(s) any(c("size","radius") %in% s$aesthetics),
                                logical(1)))
    rm_idx_col  <- which(vapply(cnet$scales$scales,
                                function(s) any(c("colour","color") %in% s$aesthetics),
                                logical(1)))
    rm_idx <- unique(c(rm_idx_size, rm_idx_col))
    if (length(rm_idx)) cnet$scales$scales <- cnet$scales$scales[-rm_idx]

    # ---- hide built-in node layers (keeps edges) ----
    node_layers <- which(sapply(cnet$layers, function(l)
      inherits(l$geom, "GeomPoint") || inherits(l$geom, "GeomNodePoint")))
    if (length(node_layers)) {
      for (i in node_layers) {
        cnet$layers[[i]]$aes_params$alpha <- 0
        cnet$layers[[i]]$aes_params$size  <- 0
      }
    }

    # ---- coordinates ----
    node_df      <- cnet$data
    cat_df       <- subset(node_df, name %in% cat_names)

    # ALL gene nodes as dots (even if not labeled)
    gene_all_df  <- subset(node_df, !(name %in% cat_names))
    gene_all_df$fc <- FC[gene_all_df$name]
    # only top 25% get labels
    gene_lab_df  <- subset(gene_all_df, name %in% label_genes)

    # absolute sizes (identity scale)
    cat_df$..sz      <- 2    # hubs
    gene_all_df$..sz <- 1     # all gene dots

    # ---- redraw dots (ALL genes), then hubs on top ----
    cnet <- cnet +
      ggplot2::geom_point(
        data = gene_all_df,
        aes(x = x, y = y, color = fc, size = ..sz),
        inherit.aes = FALSE, alpha = 0.95, show.legend = FALSE
      ) +
      ggplot2::geom_point(
        data = cat_df,
        aes(x = x, y = y, size = ..sz),
        inherit.aes = FALSE, shape = 21, fill = "grey95",
        color = "grey40", stroke = 0.45, show.legend = FALSE
      ) +
      scale_size_identity() +
      # re-add a clean color scale (replaces any prior one we stripped)
      scale_color_gradient2(
        name = "log2FC", low = "blue", mid = "white", high = "red",
        midpoint = 0, na.value = "grey80"
      ) +
      guides(
        color = guide_colorbar(title.position = "top",
                               barheight = unit(2, "mm"),
                               barwidth  = unit(12, "mm"))
      ) +
      theme(
        legend.position = "bottom",
        legend.title = element_text(size = 5),
        legend.text  = element_text(size = 5)
      )

    # ---- labels (sub-6 pt; categories smaller than genes) ----
    cat_size_pt  <- 4.2; gene_size_pt <- 4.8
    cat_size_gg  <- cat_size_pt  * pt2gg
    gene_size_gg <- gene_size_pt * pt2gg

    cnet <- cnet +
      ggrepel::geom_text_repel(
        data = cat_df,
        aes(x = x, y = y, label = name),
        inherit.aes = FALSE, size = cat_size_gg,
        segment.size = 0.03, max.overlaps = Inf,
        box.padding = 0.015, point.padding = 0.04,
        min.segment.length = 0, force = 0.02
      ) +
      ggrepel::geom_text_repel(
        data = gene_lab_df,
        aes(x = x, y = y, label = name),
        inherit.aes = FALSE, size = gene_size_gg,
        segment.size = 0.03, max.overlaps = Inf,
        box.padding = 0.015, point.padding = 0.035,
        min.segment.length = 0, force = 0.015
      )

    ggsave(
      paste0(output_dir, ct, "_DE_net.svg"),
      cnet, device = svglite::svglite,
      width = 1.5, height = 2, units = "in",
      bg = "white", limitsize = FALSE
    )
  }
}

# Run cervix/CD8/03_SingleR_annotation.R first; it generates this input file.
cervix.CD8T.integrated <- readRDS("/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_CD8T_normalized_integrated_clustered_prepMarker_SingleR.RDS")

enrich_category <- list("Naive CD8 T cells" = c("The NLRP3 inflammasome", "RIP-mediated NFkB activation via ZBP1",
                                                "TAK1-dependent IKK and NF-kappa-B activation"),
                        "Central memory CD8 T cells" = c("Signaling by Interleukins",
                                                         "Interferon alpha/beta signaling",
                                                         "Downregulation of TGF-beta receptor signaling"),
                        "Effector memory CD8 T cells" = c("Autophagy",
                                                         "Signaling by Interleukins",
                                                         "Downregulation of TGF-beta receptor signaling"),
                        "Terminal effector CD8 T cells" = c("Signaling by Interleukins",
                                                            "Downregulation of TGF-beta receptor signaling",
                                                            "Interferon Signaling"))

DE.volcano.plot(ct = "Central memory CD8 T cells",
                title = "CD8+ TCM cells",
                gene.bold =  NA,
                net.plot = TRUE)

DE.volcano.plot(ct = "Effector memory CD8 T cells",
                title = "CD8+ TEM cells",
                gene.bold =  NA,
                net.plot = TRUE)

DE.volcano.plot(ct = "Terminal effector CD8 T cells",
                title = "CD8+ TTE cells",
                gene.bold =  NA,
                net.plot = TRUE)

# CD8 trajectory UMAP and heatmap
library(SingleCellExperiment)

# Run cervix/CD8/03_SingleR_annotation.R first; it generates this annotated CD8 object.
cervix.CD8T.integrated <- readRDS("/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_CD8T_normalized_integrated_clustered_prepMarker_SingleR.RDS")

# Run cervix/CD8/07_trajectory_condition_DE.R first; it saves the model, condition test and gene clusters below.
sls.model.CD8 <- readRDS("/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_CD8_slingshot_tradeSeq_model.RDS")
trajectory_dir <- "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Cervix/Cervix_CD8T_Analysis/Cervix_CD8T_trajectory/Cervix_CD8T_trajectory_SingleR"
condRes <- readRDS(file.path(trajectory_dir, "cervix_CD8lineage_condTest.RDS"))
infected.DE.Traj.cluster.df <- readRDS(file.path(trajectory_dir, "infected_DE_Traj_cluster.RDS"))

# Recreate the original PCA-based Slingshot curve and project it onto the saved UMAP.
DefaultAssay(cervix.CD8T.integrated) <- "SCT"
cervix.CD8T.integrated.SCE <- as.SingleCellExperiment(cervix.CD8T.integrated, assay = "SCT")
reducedDim(cervix.CD8T.integrated.SCE, "umap") <- Embeddings(cervix.CD8T.integrated, "umap")
reducedDim(cervix.CD8T.integrated.SCE, "pca") <- Embeddings(cervix.CD8T.integrated, "pca")

cervix.CD8T.integrated.SCE <- slingshot(
  cervix.CD8T.integrated.SCE,
  reducedDim = "pca",
  clusterLabels = colData(cervix.CD8T.integrated.SCE)$SingleR.pruned,
  start.clus = "Naive CD8 T cells"
)
cervix.CD8T.embed <- embedCurves(
  cervix.CD8T.integrated.SCE,
  reducedDim(cervix.CD8T.integrated.SCE, "umap")
)
cervix.CD8T.curve <- slingCurves(cervix.CD8T.embed, as.df = TRUE)
colnames(cervix.CD8T.curve)[1:2] <- c("umap_1", "umap_2")

df.CD8 <- as.data.frame(Embeddings(cervix.CD8T.integrated, "umap")[, 1:2])
colnames(df.CD8) <- c("umap_1", "umap_2")
df.CD8$SingleR.pruned <- factor(
  cervix.CD8T.integrated$SingleR.pruned,
  levels = c("Naive CD8 T cells", "Central memory CD8 T cells",
             "Effector memory CD8 T cells", "Terminal effector CD8 T cells"),
  labels = c("TN", "TCM", "TEM", "TTE")
)
stopifnot(!anyNA(df.CD8$SingleR.pruned))

# Use the same subtype colors and labels as the CD8 figures.
trajectory_colors <- c(TN = "#F8766D", TCM = "#39B600",
                       TEM = "#D89000", TTE = "#C77CFF")
xr <- range(df.CD8$umap_1)
yr <- range(df.CD8$umap_2)
stubx <- diff(xr) * 0.06
stuby <- diff(yr) * 0.06

plot_traj <- ggplot(df.CD8, aes(umap_1, umap_2)) +
  geom_point(aes(color = SingleR.pruned), size = 0.001) +
  scale_color_manual(values = trajectory_colors, drop = FALSE) +
  theme_jci +
  coord_fixed(1, clip = "off") +
  theme(
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(2, 2, 2, 2),
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = NA),
    legend.key.size = grid::unit(0.22, "cm"),
    legend.text = element_text(size = 8)
  ) +
  annotate("segment", x = xr[1], xend = xr[1] + stubx,
           y = yr[1], yend = yr[1], linewidth = 0.3) +
  annotate("segment", x = xr[1], xend = xr[1],
           y = yr[1], yend = yr[1] + stuby, linewidth = 0.3) +
  annotate("text", x = xr[1] + stubx * 1.1, y = yr[1] - stuby * 0.25,
           label = "UMAP_1", size = 7 / ggplot2::.pt) +
  annotate("text", x = xr[1] - stubx * 0.25, y = yr[1] + stuby * 1.1,
           label = "UMAP_2", angle = 90, vjust = 0, size = 7 / ggplot2::.pt) +
  geom_path(
    data = cervix.CD8T.curve %>% arrange(Order),
    aes(x = umap_1, y = umap_2, group = Lineage),
    inherit.aes = FALSE
  ) +
  guides(color = guide_legend(title = NULL, override.aes = list(size = 3)))

ggsave(file.path(output_dir, "Fig4G_UMAP_trajectory.svg"), plot_traj,
       device = svglite::svglite, width = 3, height = 3, units = "in",
       bg = "white", limitsize = FALSE)

# Keep the original condition-test cutoff, 50 points and row-wise scaling.
condRes <- condRes %>% filter(padj < 0.05) %>% arrange(padj)
stopifnot(nrow(condRes) > 0,
          !anyDuplicated(infected.DE.Traj.cluster.df$gene),
          setequal(rownames(condRes), infected.DE.Traj.cluster.df$gene))

# The original model orders CT+ first and CT- second; check before splitting columns.
condition_levels <- levels(colData(sls.model.CD8)$tradeSeq$conditions)
if (!identical(condition_levels, c("infected", "uninfected"))) {
  stop("Expected model condition levels: infected, uninfected. Check the saved CD8 model.")
}
yhatSmooth2 <- predictSmooth(
  sls.model.CD8, gene = rownames(condRes), nPoints = 50, tidy = FALSE
)
if (ncol(yhatSmooth2) != 100L) {
  stop("Expected one CD8 lineage with 50 points for each of the two conditions.")
}
yhatSmooth2Scaled <- t(scale(t(yhatSmooth2)))
if (any(!is.finite(yhatSmooth2Scaled))) {
  stop("Non-finite scaled predictions; check the saved model and significant genes.")
}

# Reuse the five gene clusters saved by the analysis; match the original CD8 figure order.
make_trajectory_heatmap <- function(columns, title) {
  plot_data <- as.data.frame(yhatSmooth2Scaled[, columns, drop = FALSE])
  colnames(plot_data) <- as.character(1:50)
  plot_data$gene <- rownames(plot_data)
  plot_data <- plot_data %>%
    pivot_longer(-gene, names_to = "pseudotime", values_to = "scaled_predicted_expression") %>%
    mutate(pseudotime = as.numeric(pseudotime)) %>%
    left_join(infected.DE.Traj.cluster.df, by = "gene") %>%
    mutate(cluster = factor(cluster, levels = c(4, 1, 3, 5, 2)))
  stopifnot(!anyNA(plot_data$cluster))

  ggplot(plot_data, aes(pseudotime, gene, fill = scaled_predicted_expression)) +
    geom_tile() +
    scale_fill_gradient2(
      name = "Predicted\nExpression",
      low = "blue", mid = "white", high = "red",
      midpoint = 0, limits = c(-3, 7)
    ) +
    facet_grid(rows = vars(cluster), space = "free", scales = "free") +
    labs(title = title) +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_text(size = 8),
      axis.title.y = element_text(size = 8),
      plot.title = element_text(size = 8),
      legend.position = "bottom",
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 6),
      legend.box.margin = margin(2, 0, 2, 0),
      legend.margin = margin(0, 0, 0, 0),
      legend.spacing.x = grid::unit(2, "mm"),
      legend.spacing.y = grid::unit(2, "mm")
    ) +
    guides(fill = guide_colorbar(
      direction = "horizontal", title.position = "left",
      barheight = grid::unit(3, "mm"), barwidth = grid::unit(17, "mm"),
      ticks = FALSE
    ))
}

infected_heatmap_plot <- make_trajectory_heatmap(1:50, "   CT+")
uninfected_heatmap_plot <- make_trajectory_heatmap(51:100, "   CT-") +
  theme(plot.title = element_text(face = "bold"))
p_save <- ggarrange(
  uninfected_heatmap_plot, infected_heatmap_plot,
  ncol = 2, common.legend = TRUE, legend = "bottom"
)

ggsave(file.path(output_dir, "Fig4H_trajectory_heatmap.svg"), p_save,
       device = svglite::svglite, width = 3, height = 3, units = "in",
       bg = "white", limitsize = FALSE)
ggsave(file.path(output_dir, "Fig4H_trajectory_heatmap.tiff"), p_save,
       width = 8, height = 5, units = "in", dpi = 600)
