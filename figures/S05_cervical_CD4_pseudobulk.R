# S05 cervical CD4 T cell pseudobulk analysis

library(ggplot2)
library(Seurat)
library(SeuratObject)
library(Matrix)
library(edgeR)
library(ggrepel)
library(statmod)
library(grid)
library(svglite)

output_dir <- "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Figures/SuppFig5"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

font_family <- "Arial"

theme_jci <- theme_bw(base_size = 7, base_family = font_family) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.25),
    axis.ticks.length = unit(1.5, "pt"),
    axis.text = element_text(size = 6.5, color = "black"),
    axis.title = element_text(size = 7, face = "plain"),
    strip.background = element_blank(),
    strip.text = element_text(size = 6.5, margin = margin(1, 0, 2, 0)),
    plot.title = element_blank(),
    plot.margin = margin(3, 3, 3, 3)
  )

save_svg <- function(plot, file, width_in, height_in, bg = "white") {
  ggsave(
    filename = file,
    plot = plot,
    device = svglite::svglite,
    width = width_in, height = height_in, units = "in",
    bg = bg, limitsize = FALSE
  )
}

condition_cols <- c(uninfected = "#2979A5", infected = "#CF6045")

plot_genes <- c(
  "IFNG", "CCL20", "CXCL13", "CXCR6", "CCL5", "CXCR3",
  "CD69", "BHLHE40", "KLF2", "LEF1", "TCF7"
)

format_p <- function(p) {
  formatC(p, format = "g", digits = 2)
}

# =========================
# CD4 T cells: participant QC
# =========================
# Read metadata before selecting CD4 T cells.
obj <- readRDS(
  "/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_all_filtered_normalized_integrated_clustered_prepMarker_majorCT.RDS"
)
md <- obj[[]]
md$participant <- trimws(as.character(md$ID))
md$group <- tolower(trimws(as.character(md$condition)))
ids <- sort(unique(md$participant))
sample_info <- data.frame(participant = ids, stringsAsFactors = FALSE)
sample_info$condition <- md$group[match(ids, md$participant)]

is_cd4 <- !is.na(md$major_celltypes) & md$major_celltypes == "CD4_T"
sample_info$n_CD4 <- as.integer(table(factor(md$participant[is_cd4], levels = ids)))
sample_info$condition <- factor(sample_info$condition, levels = c("uninfected", "infected"))
sample_info$retained <- sample_info$n_CD4 >= 10L
cd4 <- subset(obj, cells = rownames(md)[is_cd4])
cd4$pb_participant <- md[colnames(cd4), "participant"]
rm(obj)
invisible(gc())

# =========================
# Raw RNA counts and participant pseudobulk
# =========================
counts <- GetAssayData(cd4, assay = "RNA", slot = "counts")
counts <- counts[, colnames(cd4), drop = FALSE]
targets_present <- plot_genes[plot_genes %in% rownames(counts)]
ret <- sample_info[sample_info$retained, , drop = FALSE]
rownames(ret) <- ret$participant
keep_cells <- cd4$pb_participant %in% ret$participant
aggregation <- Matrix::sparseMatrix(i = seq_len(sum(keep_cells)),
                                    j = match(cd4$pb_participant[keep_cells], ret$participant), x = 1,
                                    dims = c(sum(keep_cells), nrow(ret)),
                                    dimnames = list(colnames(cd4)[keep_cells], ret$participant))
pb <- as.matrix(counts[, keep_cells, drop = FALSE] %*% aggregation)

# =========================
# Pseudobulk differential expression
# =========================
# Filter/normalize/model using all expressed genes, not only the target panel.
# This is a between-participant unpaired comparison; no ID fixed effect.
design <- model.matrix(~condition, data = ret)
y0 <- edgeR::DGEList(counts = pb, samples = ret)
gene_keep <- edgeR::filterByExpr(y0, design = design)
y <- y0[gene_keep, , keep.lib.sizes = FALSE]
y <- edgeR::calcNormFactors(y, method = "TMM")
# Keep the original edgeR workflow, including NB dispersion estimation.
# edgeR v3 requires this estimation step before quasi-likelihood fitting.
y <- edgeR::estimateDisp(y, design, robust = TRUE)
fit <- edgeR::glmQLFit(y, design, robust = TRUE)
qlf <- edgeR::glmQLFTest(fit, coef = "conditioninfected")
all_de <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
all_de$gene <- rownames(all_de)
names(all_de)[names(all_de) == "FDR"] <- "FDR_all_tested_genes"
ret$filtered_library_size <- y$samples$lib.size
ret$TMM_factor <- y$samples$norm.factors
ret$effective_library_size <- ret$filtered_library_size * ret$TMM_factor

panel <- data.frame(gene = plot_genes, stringsAsFactors = FALSE)
for (v in setdiff(names(all_de), "gene")) panel[[v]] <- all_de[[v]][match(panel$gene, all_de$gene)]
panel$BH_adjusted_p_value <- panel$FDR_all_tested_genes

# Descriptive logCPM for all present targets, including filtered ones. Use the
# same effective library sizes as the fitted model; do not re-normalize the displayed genes.
target_logcpm <- edgeR::cpm(pb[targets_present, , drop = FALSE],
                            lib.size = ret$effective_library_size, log = TRUE, prior.count = 2)
expression_long <- data.frame(
  gene = rep(targets_present, times = nrow(ret)),
  participant = rep(ret$participant, each = length(targets_present)),
  logCPM = as.vector(target_logcpm), stringsAsFactors = FALSE)
expression_long$condition <- ret$condition[match(expression_long$participant, ret$participant)]

# =========================
# Original single-cell comparison
# =========================
# Import original cell-level FindMarkers results exactly; match by gene symbol.
# The CSV's first unnamed column is the gene row-name column from write.csv.
sc <- read.csv(
  "/work/users/s/y/syzhang/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Cervix/Cervix_DE_INFvsUNINF/All_INFvsUNINF/CD4_T_INFvsUNINF_All_Wilcox.csv",
  row.names = 1, check.names = FALSE, stringsAsFactors = FALSE
)
sc$gene <- trimws(rownames(sc))
# Match pseudobulk-tested genes to the original single-cell results.
all_comparison <- all_de
sc_index <- match(all_comparison$gene, sc$gene)
all_comparison$singlecell_log2FC <- sc$avg_log2FC[sc_index]
all_comparison$singlecell_adjusted_p <- sc$p_val_adj[sc_index]
all_comparison$within_scatter_limits <-
  is.finite(all_comparison$singlecell_log2FC) & is.finite(all_comparison$logFC) &
  all_comparison$singlecell_log2FC >= -5 &
  all_comparison$singlecell_log2FC <= 15 &
  all_comparison$logFC >= -5 &
  all_comparison$logFC <= 15
all_scatter_data <- subset(
  all_comparison,
  singlecell_adjusted_p < 0.05 & is.finite(singlecell_log2FC) & is.finite(logFC)
)

# =========================
# Supplemental Figure 5A: single-cell versus pseudobulk log2FC
# =========================
d <- all_scatter_data
selected <- d[d$gene %in% plot_genes, , drop = FALSE]
labels <- selected[selected$within_scatter_limits, , drop = FALSE]

# Spearman correlation uses all SC-significant genes with both fold changes.
spearman_test <- cor.test(
  d$singlecell_log2FC, d$logFC,
  method = "spearman", exact = FALSE
)
correlation_label <- sprintf(
  "Spearman \u03c1 = %.2f\np = %s",
  unname(spearman_test$estimate),
  format.pval(spearman_test$p.value, digits = 2, eps = .Machine$double.eps)
)
correlation_label <- sub("p = <", "p < ", correlation_label, fixed = TRUE)

CD4_p1 <- ggplot(d, aes(singlecell_log2FC, logFC)) +
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.25) +
  geom_vline(xintercept = 0, color = "grey80", linewidth = 0.25) +
  geom_abline(
    slope = 1, intercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.25
  ) +
  geom_point(color = "grey70", shape = 16, size = 0.7) +
  geom_point(
    data = selected, color = "black", shape = 16, size = 0.7
  ) +
  ggrepel::geom_text_repel(
    data = labels, aes(label = gene), seed = 17,
    family = font_family, color = "black", size = 7 / ggplot2::.pt,
    max.overlaps = Inf, max.time = 2, max.iter = 20000,
    box.padding = 0.25, point.padding = 0.15, min.segment.length = 0,
    segment.size = 0.2, segment.color = "black",
    xlim = c(-5, 15), ylim = c(-5, 15)
  ) +
  annotate(
    "text", x = 14.5, y = -4.5, label = correlation_label,
    hjust = 1, vjust = 0, family = font_family,
    size = 7 / ggplot2::.pt, color = "black"
  ) +
  scale_x_continuous(breaks = seq(-5, 15, 5)) +
  scale_y_continuous(breaks = seq(-5, 15, 5)) +
  coord_fixed(
    ratio = 1, xlim = c(-5, 15), ylim = c(-5, 15), expand = FALSE
  ) +
  labs(
    x = "Single-cell log2FC (CT+ / CT-)",
    y = "Pseudobulk log2FC (CT+ / CT-)"
  ) +
  theme_jci +
  theme(legend.position = "none")

save_svg(
  CD4_p1,
  file.path(output_dir, "CD4_singlecell_vs_pseudobulk.svg"),
  width_in = 3.5, height_in = 3
)

# =========================
# Supplemental Figure 6B: participant expression grid
# =========================
boxplot_genes <- plot_genes
d <- expression_long[expression_long$gene %in% plot_genes, , drop = FALSE]
d$gene <- factor(d$gene, levels = boxplot_genes)
ann <- panel[match(boxplot_genes, panel$gene), , drop = FALSE]
q <- ann$BH_adjusted_p_value
ann$label <- paste0(ann$gene, "\nadj. pval. = ", format_p(q))
labels <- setNames(ann$label, ann$gene)

CD4_p2 <- ggplot(d, aes(condition, logCPM, color = condition)) +
  geom_boxplot(
    aes(fill = condition), width = 0.5,
    outlier.shape = NA, alpha = 0.2, linewidth = 0.3
  ) +
  geom_point(
    position = position_jitter(width = 0.12, height = 0, seed = 17),
    size = 0.95
  ) +
  facet_wrap(
    ~gene, ncol = 6, scales = "free", drop = FALSE,
    labeller = labeller(gene = labels)
  ) +
  scale_color_manual(values = condition_cols) +
  scale_fill_manual(values = condition_cols) +
  scale_x_discrete(labels = c(uninfected = "CT-", infected = "CT+")) +
  labs(x = NULL, y = "TMM-normalized log2 CPM") +
  theme_jci +
  theme(
    legend.position = "none",
    panel.spacing.x = unit(4, "pt"),
    panel.spacing.y = unit(5, "pt")
  )

save_svg(
  CD4_p2,
  file.path(output_dir, "CD4_pseudobulk_expression.svg"),
  width_in = 7, height_in = 3.8
)
