# helper gene heatmap

library(dplyr)
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(tibble)
library(grid)
library(svglite)

# ============================================================
# SETTINGS
# ============================================================

output_dir <- "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Figures/Fig6/"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cervix.CD4T.integrated <- readRDS(
  "/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_CD4T_normalized_integrated_clustered_prepMarker_SingleR_manual.RDS"
)

obj <- cervix.CD4T.integrated

celltype_col  <- "SingleR.pruned"
infection_col <- "condition"

celltype_levels_raw <- c(
  "Th17 cells",
  "Th1/Th17 cells",
  "Th1 cells"
)

infection_levels_raw <- c(
  "uninfected",
  "infected"
)

# Preserve the original group order
group_levels <- c(
  "Th17_CT-",
  "Th17_CT+",
  "Th1/17_CT-",
  "Th1/17_CT+",
  "Th1_CT-",
  "Th1_CT+"
)

# ============================================================
# SELECTED GENES — FIXED ORIGINAL ORDER
# ============================================================

genes_shared <- c(
  "CXCR6", "BHLHE40", "IRF4", "BATF", "ICOS",
  "STAT4", "IRF9", "RELB", "NFKB2", "IL2RB",
  "NAMPT", "PMEPA1", "SKIL",
  "KLF2", "LEF1", "TCF7", "JUNB"
)

genes_th17 <- c(
  "AHR", "CCL20", "IL1R1", "IL18R1",
  "MIR155HG", "IL4I1"
)

genes_th117 <- c(
  "IFNG", "IL12RB2", "MX1", "ISG15", "IRF7",
  "TNF", "LTA", "PDCD1", "TIGIT",
  "SMAD7", "SKI"
)

genes_th1 <- c(
  "CCR5", "OASL", "CTLA4", "LAG3", "CXCL13"
)

# Combine genes in their original order without visual splits
genes_heatmap <- unique(c(
  genes_shared,
  genes_th17,
  genes_th117,
  genes_th1
))

# ============================================================
# PREPARE CELL METADATA
# ============================================================

meta_sub <- obj@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::filter(
    .data[[celltype_col]] %in% celltype_levels_raw,
    .data[[infection_col]] %in% infection_levels_raw
  ) %>%
  dplyr::mutate(
    celltype = dplyr::recode(
      .data[[celltype_col]],
      "Th17 cells" = "Th17",
      "Th1/Th17 cells" = "Th1/17",
      "Th1 cells" = "Th1"
    ),
    infection = dplyr::recode(
      .data[[infection_col]],
      "uninfected" = "CT-",
      "infected" = "CT+"
    ),
    celltype = factor(
      celltype,
      levels = c("Th17", "Th1/17", "Th1")
    ),
    infection = factor(
      infection,
      levels = c("CT-", "CT+")
    ),
    group = factor(
      paste(celltype, infection, sep = "_"),
      levels = group_levels
    )
  ) %>%
  dplyr::arrange(group)

cells_use <- meta_sub$cell

if (length(cells_use) == 0L) {
  stop("No cells matched the selected cell types and conditions.")
}

cat("\nCells included in each group:\n")

print(
  meta_sub %>%
    dplyr::count(celltype, infection, .drop = FALSE)
)

cat("\nTotal cells:", length(cells_use), "\n")

missing_groups <- setdiff(
  group_levels,
  as.character(unique(meta_sub$group))
)

if (length(missing_groups) > 0L) {
  stop(
    "Cannot calculate average expression for groups with no cells: ",
    paste(missing_groups, collapse = ", ")
  )
}

# ============================================================
# GET SCT-NORMALIZED EXPRESSION
# ============================================================

DefaultAssay(obj) <- "SCT"

sct_data <- GetAssayData(
  obj,
  assay = "SCT",
  layer = "data"
)

# Preserve the original gene order
genes_present <- genes_heatmap[
  genes_heatmap %in% rownames(sct_data)
]

genes_missing <- setdiff(
  genes_heatmap,
  genes_present
)

if (length(genes_missing) > 0L) {
  message(
    "Genes not found in SCT data: ",
    paste(genes_missing, collapse = ", ")
  )
}

if (length(genes_present) == 0L) {
  stop("None of the selected genes were found in SCT data.")
}

expr <- sct_data[
  genes_present,
  cells_use,
  drop = FALSE
]

# ============================================================
# AVERAGE EXPRESSION PER GROUP
# ============================================================

# Average SCT-normalized expression across all cells within
# each cell-type/infection group, including zero-expression cells.
#
# Result:
#   rows    = genes
#   columns = the six groups in the fixed order above

expr_avg <- do.call(
  cbind,
  lapply(group_levels, function(group_name) {
    
    group_cells <- meta_sub$cell[
      as.character(meta_sub$group) == group_name
    ]
    
    group_expr <- expr[
      ,
      group_cells,
      drop = FALSE
    ]
    
    if (inherits(group_expr, "Matrix")) {
      Matrix::rowMeans(group_expr)
    } else {
      rowMeans(group_expr)
    }
  })
)

rownames(expr_avg) <- genes_present
colnames(expr_avg) <- group_levels

# ============================================================
# Z-SCORE GROUP AVERAGES BY GENE
# ============================================================

# Z-score each gene across the six group averages.
#
# IMPORTANT:
# The actual Z-scores remain uncapped.
# Only the heatmap COLOR SCALE below is capped at -1.5 to +1.5.

expr_scaled <- t(
  scale(
    t(expr_avg)
  )
)

# Remove genes with undefined Z-scores, such as genes with
# identical average expression across all six groups.

keep_genes <- apply(
  expr_scaled,
  1,
  function(x) all(is.finite(x))
)

if (any(!keep_genes)) {
  message(
    "Removed genes with no variance or non-finite values: ",
    paste(
      rownames(expr_scaled)[!keep_genes],
      collapse = ", "
    )
  )
}

expr_scaled <- expr_scaled[
  keep_genes,
  ,
  drop = FALSE
]

if (nrow(expr_scaled) == 0L) {
  stop("No genes remain after scaling group averages.")
}

# ============================================================
# GROUP ANNOTATIONS
# ============================================================

group_meta <- meta_sub %>%
  dplyr::distinct(
    group,
    celltype,
    infection
  ) %>%
  dplyr::arrange(group)

stopifnot(
  identical(
    as.character(group_meta$group),
    colnames(expr_scaled)
  )
)

column_celltype  <- group_meta$celltype
column_infection <- group_meta$infection

# ============================================================
# COLORS — COLOR SCALE CAPPED AT -1.5 TO +1.5
# ============================================================

# The underlying Z-scores remain unchanged.
#
# Any value <= -1.5 receives the darkest blue.
# Any value >= +1.5 receives the darkest red.

heat_col <- circlize::colorRamp2(
  seq(-1.5, 1.5, length.out = 4),
  rev(RColorBrewer::brewer.pal(n = 4, name = "RdBu"))
)

celltype_colors <- c(
  "Th17"   = "#FDB462",
  "Th1/17" = "#80B1D3",
  "Th1"    = "#FB8072"
)

# Purple distinguishes CT+ from the red heatmap colors
infection_colors <- c(
  "CT-" = "grey80",
  "CT+" = "#7B3294"
)

# ============================================================
# TOP ANNOTATION
# ============================================================

ha <- HeatmapAnnotation(
  `Cell type` = column_celltype,
  `Infection Status` = column_infection,
  
  col = list(
    `Cell type` = celltype_colors,
    `Infection Status` = infection_colors
  ),
  
  simple_anno_size = unit(2.5, "mm"),
  
  annotation_name_gp = gpar(
    fontsize = 7
  ),
  
  annotation_legend_param = list(
    
    `Cell type` = list(
      title = "Cell type",
      at = c(
        "Th17",
        "Th1/17",
        "Th1"
      ),
      ncol = 1,
      title_position = "topleft",
      title_gp = gpar(
        fontsize = 7
      ),
      labels_gp = gpar(
        fontsize = 7
      )
    ),
    
    `Infection Status` = list(
      title = "Infection Status",
      at = c(
        "CT-",
        "CT+"
      ),
      ncol = 1,
      title_position = "topleft",
      title_gp = gpar(
        fontsize = 7
      ),
      labels_gp = gpar(
        fontsize = 7
      )
    )
  )
)

# ============================================================
# CREATE HEATMAP
# ============================================================

ht <- Heatmap(
  expr_scaled,
  
  name = "Z-score",
  col = heat_col,
  
  # ----------------------------------------------------------
  # Columns
  # ----------------------------------------------------------
  
  # One column per group, preserving the specified group order
  cluster_columns = FALSE,
  
  column_order = seq_len(
    ncol(expr_scaled)
  ),
  
  show_column_dend = FALSE,
  show_column_names = TRUE,
  
  column_labels = gsub(
    "_",
    " ",
    colnames(expr_scaled),
    fixed = TRUE
  ),
  
  column_names_rot = 45,
  
  column_names_gp = gpar(
    fontsize = 7
  ),
  
  # ----------------------------------------------------------
  # Rows
  # ----------------------------------------------------------
  
  # Preserve the original selected-gene order
  cluster_rows = FALSE,
  
  row_order = seq_len(
    nrow(expr_scaled)
  ),
  
  show_row_dend = FALSE,
  show_row_names = TRUE,
  
  row_names_side = "left",
  
  row_names_gp = gpar(
    fontsize = 7
  ),
  
  # ----------------------------------------------------------
  # Annotation
  # ----------------------------------------------------------
  
  top_annotation = ha,
  
  # ----------------------------------------------------------
  # Heatmap legend
  # ----------------------------------------------------------
  
  heatmap_legend_param = list(
    title = "Mean expression\nZ-score",
    
    at = c(
      -1.5,
      -1,
      0,
      1,
      1.5
    ),
    
    labels = c(
      "-1.5",
      "-1",
      "0",
      "1",
      "1.5"
    ),
    
    direction = "vertical",
    
    title_position = "topleft",
    
    title_gp = gpar(
      fontsize = 7
    ),
    
    labels_gp = gpar(
      fontsize = 7
    ),
    
    legend_height = unit(
      25,
      "mm"
    )
  ),
  
  border = FALSE,
  use_raster = FALSE
)

# ============================================================
# SAVE HEATMAP
# ============================================================

# The matrix itself is NOT capped.
# The color representation saturates at Z = -1.5 and +1.5.

svglite::svglite(
  file.path(
    output_dir,
    paste0(
      "Th17_Th117_Th1_SelectedGenes_",
      "GroupAverage_FixedOrder_",
      "ColorCapped1.5_Heatmap.svg"
    )
  ),
  width = 5,
  height = 5.5
)

draw(
  ht,
  merge_legends = TRUE,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()