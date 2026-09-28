# innate B DE pathways

library(dplyr)
library(ggplot2)
library(Seurat)
library(scCustomize)
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
library(ggraph)

output_dir <- "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Figures/Fig3/"

# DE analysis in cervix

cervix.integrated <- readRDS("/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_all_filtered_normalized_integrated_clustered_prepMarker_majorCT.RDS")

DefaultAssay(object = cervix.integrated) <- "SCT"

cervix.integrated$major_celltypes[which(cervix.integrated$major_celltypes%in% c("CD14_Mono","CD16_Mono"))] <- "Monocyte"

Idents(cervix.integrated) <- "major_celltypes"

celltypes <- unique(cervix.integrated$major_celltypes)

enrich_category <- list("CD4_T" = c("Signaling by Interleukins", "Interferon alpha/beta signaling",
                                    "Interferon gamma signaling", "Toll-like Receptor Cascades"),

                        "CD8_T" = c("Interferon gamma signaling",  "Interferon alpha/beta signaling",
                                    "Signaling by Interleukins",
                                    "Costimulation by the CD28 family"),

                        "B" = c( "Interferon gamma signaling", "Interferon alpha/beta signaling",
                                 "Signaling by Interleukins"),

                        "Monocyte" = c("Interferon gamma signaling", "MHC class II antigen presentation",
                                      "Interferon alpha/beta signaling", "Antigen processing-Cross presentation"),

                        "MAIT" = c("Cell surface interactions at the vascular wall", "TNFR2 non-canonical NF-kB pathway",
                                    "Interferon alpha/beta signaling", "Signaling by Interleukins"),

                         "NK" = c("Interferon alpha/beta signaling", "Interferon gamma signaling", "Toll-like Receptor Cascades") )

DE.volcano.plot <- function(ct, title, gene.bold,
                            net.plot = FALSE,
                            net.prop.gene.show = 0.5){

  # --- constants ---
  base_pt  <- 6
  title_pt <- 8
  label_pt <- 6
  pt2gg    <- 1/ggplot2::.pt

  DE <- read.csv(paste0("/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Cervix/Cervix_DE_INFvsUNINF/All_INFvsUNINF/", ct, "_INFvsUNINF_All_Wilcox.csv"))

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

  # --- VOLCANO  ---
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
    DEplot, device = "svg",
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
      circular     = FALSE, layout = "kk",
      size_edge = 0.2
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

###  B cells
DE.volcano.plot(ct = "B",
                title = "B cells",
                gene.bold =  NULL,
                net.plot = TRUE)

### Monocytes
DE.volcano.plot(ct = "Monocyte",
                title = "Monocytes",
                gene.bold =  NULL,
                net.plot = "TRUE")
### CD4 T cells
# DE.volcano.plot(ct = "CD4_T",
#                 title = "CD4+ T cells",
#                 gene.bold =  NULL,
#                 net.plot = TRUE)

### CD8 T cells
# DE.volcano.plot(ct = "CD8_T",
#                 title = "CD8+ T cells",
#                 gene.bold =  NULL,
#                 net.plot = TRUE)

### MAIT cells
DE.volcano.plot(ct = "MAIT",
                title = "MAIT/NK T cells",
                gene.bold =  NULL,
                net.plot = TRUE)

### NK cells
DE.volcano.plot(ct = "NK",
                title = "NK cells",
                gene.bold =  NULL,
                net.plot = TRUE)
