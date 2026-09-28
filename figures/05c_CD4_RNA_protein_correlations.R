#  CD4 RNA protein correlations

library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(Seurat)
library(ggplot2)
library(stringr)

# ============================ Paths ============================
BASE <- "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Cervix/Cervix_major_cytokine"
if (!dir.exists(BASE)) dir.create(BASE, recursive = TRUE, showWarnings = FALSE)

path_map_cytgene <- file.path(BASE, "Cytokines to Gene Symbols1.xlsx")
path_cyt         <- file.path(BASE, "all_samples_batch_correct_clean_cytokine.csv")

out_dir <- "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Figures/Fig3_CD4"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================ Read Seurat ============================
seu <- readRDS(
  "/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_all_filtered_normalized_integrated_clustered_prepMarker_majorCT.RDS"
)

DefaultAssay(seu) <- "SCT"
seu$major_celltypes[seu$major_celltypes %in% c("CD14_Mono", "CD16_Mono")] <- "Monocyte"
seu <- subset(seu, subset = condition == "infected" & major_celltypes == "CD4_T")

# ============================ Plot helper ============================
make_plot <- function(dat, xvar, yvar, xlab, ylab, title) {
  x <- as.numeric(dat[[xvar]])
  y <- as.numeric(dat[[yvar]])
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]

  if (length(x) < 3 || length(unique(x)) < 2 || length(unique(y)) < 2) {
    spearman_r <- NA_real_
    spearman_p <- NA_real_
  } else {
    ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
    spearman_r <- unname(ct$estimate)
    spearman_p <- ct$p.value
  }

  rho_lab <- if (is.na(spearman_r)) {
    "Spearman ρ = NA"
  } else {
    paste0(
      "Spearman ρ = ", sprintf("%.2f", spearman_r),
      "\np = ", sprintf("%.4f", spearman_p)
    )
  }

  ggplot(dat, aes(x = .data[[xvar]], y = .data[[yvar]])) +
    geom_point(size = 1, color = "black", alpha = 0.9) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.4) +
    annotate(
      "text",
      x = Inf, y = -Inf,
      label = rho_lab,
      hjust = 1.02, vjust = -0.2,
      size = 2.0
    ) +
    labs(
      title = title,
      x = xlab,
      y = ylab
    ) +
    theme_minimal(base_size = 8) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.3),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.text = element_text(color = "black", size = 7),
      axis.title = element_text(color = "black", size = 8),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
      plot.margin = margin(t = 5, r = 5, b = 2, l = 5)
    )
}

# ============================ Read mapping for clipping bounds ============================
map_cytgene <- readxl::read_xlsx(path_map_cytgene)
cn <- names(map_cytgene)

analyte_col <- cn[grepl("analyte label", cn, ignore.case = TRUE)][1]
lower_col   <- cn[grepl("lower", cn, ignore.case = TRUE)][1]
upper_col   <- cn[grepl("upper", cn, ignore.case = TRUE)][1]

if (is.na(lower_col)) lower_col <- cn[3]
if (is.na(upper_col)) upper_col <- cn[4]

clip_tbl <- map_cytgene %>%
  transmute(
    Cytokine  = as.character(.data[[analyte_col]]),
    lower_lim = suppressWarnings(as.numeric(.data[[lower_col]])),
    upper_lim = suppressWarnings(as.numeric(.data[[upper_col]]))
  ) %>%
  filter(!is.na(Cytokine), Cytokine != "") %>%
  distinct()

# ============================ Read cytokine data ============================
cyt_raw <- read_csv(path_cyt, show_col_types = FALSE)
id_col  <- names(cyt_raw)[1]

cyt_long <- cyt_raw %>%
  rename(ID = all_of(id_col)) %>%
  mutate(ID = gsub("-", "", as.character(ID))) %>%
  pivot_longer(-ID, names_to = "Cytokine", values_to = "CytokineExpr_log2_raw") %>%
  filter(!is.na(CytokineExpr_log2_raw)) %>%
  left_join(clip_tbl, by = "Cytokine") %>%
  mutate(
    CytokineExpr = 2^CytokineExpr_log2_raw,
    CytokineExpr = if_else(!is.na(lower_lim) & CytokineExpr < lower_lim, lower_lim, CytokineExpr),
    CytokineExpr = if_else(!is.na(upper_lim) & CytokineExpr > upper_lim, upper_lim, CytokineExpr),
    CytokineExpr_log2_clip = log2(CytokineExpr)
  )

bad_pairs <- tibble::tribble(
  ~Cytokine,                  ~ID,
  "YKL40.CHI3L1..77.",        "3159",
  "CXCL11.I.TAC..19.",        "3119",
  "CXCL6.GCP.2..15.",         "3119",
  "CXCL9.MIG..47.",           "3119"
)

cyt_long <- cyt_long %>%
  anti_join(bad_pairs, by = c("Cytokine", "ID"))

# ============================ Cytokine summaries ============================
bca1_dat <- cyt_long %>%
  filter(Cytokine == "BCA.1..15.") %>%
  group_by(ID) %>%
  summarise(
    Cytokine_log2 = mean(CytokineExpr_log2_clip, na.rm = TRUE),
    .groups = "drop"
  )

cxcl_avg <- cyt_long %>%
  filter(Cytokine %in% c("CXCL9.MIG..47.", "CXCL10.IP.10..48.", "CXCL11.I.TAC..19.")) %>%
  group_by(ID) %>%
  summarise(
    Cytokine_linear_mean = mean(CytokineExpr, na.rm = TRUE),
    Cytokine_log2 = log2(Cytokine_linear_mean),
    .groups = "drop"
  )

rantes_dat <- cyt_long %>%
  filter(Cytokine == "RANTES..74.") %>%
  group_by(ID) %>%
  summarise(
    Cytokine_log2 = mean(CytokineExpr_log2_clip, na.rm = TRUE),
    .groups = "drop"
  )

# ============================ RNA mean expression in CD4_T ============================
genes_needed <- c("CXCL13", "CXCR3", "CCL5")
genes_present <- intersect(genes_needed, rownames(GetAssayData(seu, slot = "data", assay = "SCT")))

if (length(genes_present) < length(genes_needed)) {
  stop("Not all required genes (CXCL13, CXCR3, CCL5) are present in SCT assay.")
}

sc_df <- FetchData(
  seu,
  vars = c("ID", genes_needed),
  slot = "data"
)

sc_df <- sc_df %>%
  mutate(ID = gsub("-", "", as.character(ID)))

rna_mean <- sc_df %>%
  group_by(ID) %>%
  summarise(
    CXCL13_mean = mean(CXCL13, na.rm = TRUE),
    CXCR3_mean  = mean(CXCR3,  na.rm = TRUE),
    CCL5_mean   = mean(CCL5,   na.rm = TRUE),
    .groups = "drop"
  )

# ============================ Zero-expression table ============================
rna_zero_table <- sc_df %>%
  pivot_longer(
    cols = all_of(genes_needed),
    names_to = "RNA",
    values_to = "expr"
  ) %>%
  group_by(ID, RNA) %>%
  summarise(
    total_cells = sum(!is.na(expr)),
    zero_cells = sum(expr == 0, na.rm = TRUE),
    frac_zero = zero_cells / total_cells,
    .groups = "drop"
  ) %>%
  arrange(RNA, ID)

write_csv(
  rna_zero_table,
  file.path(out_dir, "CD4T_RNA_zero_fraction_by_participant.csv")
)

# ============================ Keep only participants in both datasets for each plot ============================
plot0_ids <- intersect(bca1_dat$ID, rna_mean$ID)
plot1_ids <- intersect(cxcl_avg$ID, rna_mean$ID)
plot2_ids <- intersect(rantes_dat$ID, rna_mean$ID)

plot0_dat <- bca1_dat %>%
  filter(ID %in% plot0_ids) %>%
  inner_join(
    rna_mean %>%
      filter(ID %in% plot0_ids) %>%
      select(ID, CXCL13_mean),
    by = "ID"
  )

plot1_dat <- cxcl_avg %>%
  filter(ID %in% plot1_ids) %>%
  inner_join(
    rna_mean %>%
      filter(ID %in% plot1_ids) %>%
      select(ID, CXCR3_mean),
    by = "ID"
  )

plot2_dat <- rantes_dat %>%
  filter(ID %in% plot2_ids) %>%
  inner_join(
    rna_mean %>%
      filter(ID %in% plot2_ids) %>%
      select(ID, CCL5_mean),
    by = "ID"
  )

# ============================ Make plots ============================
p0 <- make_plot(
  dat   = plot0_dat,
  xvar  = "Cytokine_log2",
  yvar  = "CXCL13_mean",
  xlab  = "Log2 CXCL13 (BCA1) pg/mL",
  ylab  = "Mean CXCL13 expression",
  title = NULL
)

p1 <- make_plot(
  dat   = plot1_dat,
  xvar  = "Cytokine_log2",
  yvar  = "CXCR3_mean",
  xlab  = "Log2 mean CXCL9/CXCL10/CXCL11 pg/mL",
  ylab  = "Mean CXCR3 expression",
  title = NULL
)

p2 <- make_plot(
  dat   = plot2_dat,
  xvar  = "Cytokine_log2",
  yvar  = "CCL5_mean",
  xlab  = "Log2 CCL5 (RANTES) pg/mL",
  ylab  = "Mean CCL5 expression",
  title = NULL
)

# ============================ Save plots as SVG ============================
ggsave(
  filename = file.path(out_dir, "CD4T_BCA1_vs_CXCL13.svg"),
  plot = p0,
  width = 2,
  height = 2
)

ggsave(
  filename = file.path(out_dir, "CD4T_CXCL9_10_11_vs_CXCR3.svg"),
  plot = p1,
  width = 2.7,
  height = 2
)

ggsave(
  filename = file.path(out_dir, "CD4T_RANTES_vs_CCL5.svg"),
  plot = p2,
  width = 2,
  height = 2
)
