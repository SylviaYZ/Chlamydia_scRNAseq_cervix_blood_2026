# export trajectory gene table

library(dplyr)
library(openxlsx)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)

select <- dplyr::select

# -----------------------------
# Read differential expression output
# -----------------------------
gene_output <- readRDS(
  "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Cervix/Cervix_CD4T_analysis/Cervix_CD4T_trajectory/Cervix_CD4T_PCA_THelper/infected_DE_Traj_cluster.RDS"
)

head(gene_output)
colnames(gene_output)

# -----------------------------
# Reorder / relabel clusters
# -----------------------------
gene_output$cluster_update <- dplyr::case_when(
  gene_output$cluster == 2 ~ 1,
  gene_output$cluster == 5 ~ 2,
  gene_output$cluster == 3 ~ 3,
  gene_output$cluster == 1 ~ 4,
  gene_output$cluster == 4 ~ 5,
  TRUE ~ NA_real_
)

# -----------------------------
# Check required columns
# -----------------------------
required_cols <- c("gene", "waldStat", "df", "pvalue", "padj", "cluster_update")

missing_cols <- setdiff(required_cols, colnames(gene_output))

if (length(missing_cols) > 0) {
  stop(
    paste(
      "The following required columns are missing from gene_output:",
      paste(missing_cols, collapse = ", ")
    )
  )
}

# -----------------------------
# Helper functions
# -----------------------------

# Valid Excel sheet names
sanitize_sheet <- function(x) {
  x <- gsub("[:\\\\/\\?\\*\\[\\]]", "_", as.character(x))
  ifelse(nchar(x) > 25, paste0(substr(x, 1, 25), "..."), x)
}

# Format p-values for Excel display
# Important: this returns character values, so Excel will not display very small values as 0.
format_pval <- function(x, digits = 3) {
  ifelse(
    is.na(x),
    NA,
    ifelse(
      x == 0,
      "<1E-16",
      format(x, scientific = TRUE, digits = digits)
    )
  )
}

# -----------------------------
# Output paths
# -----------------------------
output_dir <- "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Cervix/Cervix_CD4T_analysis/Cervix_CD4T_trajectory/Cervix_CD4T_PCA_THelper"

de_output_file <- file.path(output_dir, "DE_by_cluster_results.xlsx")

reactome_dir <- file.path(output_dir, "Reactome_pathway")

if (!dir.exists(reactome_dir)) {
  dir.create(reactome_dir, recursive = TRUE)
}

# -----------------------------
# Create DE Excel workbook
# -----------------------------
wb <- createWorkbook()

# Keep only rows with non-NA cluster_update
df <- gene_output %>%
  dplyr::filter(!is.na(.data$cluster_update))

# Excel styles
header_style <- createStyle(
  textDecoration = "bold",
  halign = "center",
  border = "Bottom"
)

wald_style <- createStyle(numFmt = "0.0000")

integer_style <- createStyle(numFmt = "0")

text_style <- createStyle(numFmt = "@")

footnote_style <- createStyle(
  textDecoration = "italic",
  fontColour = "#666666",
  wrapText = TRUE
)

# -----------------------------
# Loop over clusters and write DE tables
# -----------------------------
for (cl in sort(unique(df$cluster_update))) {

  sheet_name <- paste0("cluster_", sanitize_sheet(cl))
  addWorksheet(wb, sheet_name)

  out <- df %>%
    dplyr::filter(.data$cluster_update == cl) %>%
    dplyr::arrange(.data$padj) %>%
    dplyr::select(dplyr::all_of(c("gene", "waldStat", "df", "pvalue", "padj"))) %>%
    dplyr::mutate(
      pvalue = format_pval(.data$pvalue),
      padj = format_pval(.data$padj)
    )

  colnames(out) <- c(
    "Gene",
    "Wald statistic",
    "Degrees of freedom",
    "Nominal p-value",
    "Adjusted p-value"
  )

  writeData(
    wb,
    sheet = sheet_name,
    x = out,
    startRow = 1,
    startCol = 1
  )

  # Header formatting
  addStyle(
    wb,
    sheet = sheet_name,
    style = header_style,
    rows = 1,
    cols = 1:ncol(out),
    gridExpand = TRUE
  )

  # Apply number/text formatting only if table has rows
  if (nrow(out) > 0) {

    # Wald statistic
    addStyle(
      wb,
      sheet = sheet_name,
      style = wald_style,
      rows = 2:(nrow(out) + 1),
      cols = which(colnames(out) %in% c("Wald statistic")),
      gridExpand = TRUE
    )

    # Degrees of freedom
    addStyle(
      wb,
      sheet = sheet_name,
      style = integer_style,
      rows = 2:(nrow(out) + 1),
      cols = which(colnames(out) %in% c("Degrees of freedom")),
      gridExpand = TRUE
    )

    # P-values are text because zero values are displayed as <1E-16
    addStyle(
      wb,
      sheet = sheet_name,
      style = text_style,
      rows = 2:(nrow(out) + 1),
      cols = which(colnames(out) %in% c("Nominal p-value", "Adjusted p-value")),
      gridExpand = TRUE
    )
  }

  # Add table footnote
  footnote_row <- nrow(out) + 3

  footnote_text <- paste(
    "Footnote:",
    "P-values reported as <1E-16 were stored as 0 due to numerical precision limits in the tradeSeq Wald test p-value calculation",
    "and should be interpreted as smaller than 1 × 10^-16, not as true zero.",
    "Wald statistic represents the Wald test statistic from tradeSeq conditionTest;",
    "Degrees of freedom represents the degrees of freedom used in the test;",
    "Nominal p-value is the unadjusted p-value;",
    "Adjusted p-value is the Benjamini-Hochberg adjusted p-value."
  )

  writeData(
    wb,
    sheet = sheet_name,
    x = footnote_text,
    startRow = footnote_row,
    startCol = 1
  )

  mergeCells(
    wb,
    sheet = sheet_name,
    cols = 1:ncol(out),
    rows = footnote_row
  )

  addStyle(
    wb,
    sheet = sheet_name,
    style = footnote_style,
    rows = footnote_row,
    cols = 1,
    gridExpand = TRUE
  )

  # Improve readability
  setColWidths(
    wb,
    sheet = sheet_name,
    cols = 1:ncol(out),
    widths = "auto"
  )

  freezePane(
    wb,
    sheet = sheet_name,
    firstActiveRow = 2
  )
}

# Save DE workbook
saveWorkbook(
  wb,
  file = de_output_file,
  overwrite = TRUE
)
