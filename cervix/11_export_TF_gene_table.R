# Export broad cervical CD4/CD8 TF results for Supplementary Table 5

suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(writexl)
})

# File paths

# Run cervix/09_collectri_viper_inference.R first; it saves this VIPER activity RDS.

viper_rds <- paste0(
  "/Chlamydia-SingleCell/",
  "SCTransformV2 Oct26 2023/Cervix/Cervix_major_TF_VIPER/",
  "Cervix_major_TF_VIPER.rds"
)

# Run cervix/05_annotate_major_cell_types.R first; it saves this annotated cervical RDS.

seurat_rds <- paste0(
  "/Chlamydia-SingleCell/",
  "CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/",
  "cervix_all_filtered_normalized_integrated_clustered_",
  "prepMarker_majorCT.RDS"
)

output_file <- paste0(
  "/Chlamydia-SingleCell/",
  "SCTransformV2 Oct26 2023/Cervix/Cervix_major_TF_VIPER/",
  "VIPER_CD4_CD8_Wilcoxon_results.xlsx"
)

# Load data

# Run cervix/09_collectri_viper_inference.R first; it generates this input file.
acts <- readRDS(viper_rds)

# Run cervix/05_annotate_major_cell_types.R first; it generates this input file.
cervix.integrated <- readRDS(seurat_rds)

# Match the cell-type naming used in the original VIPER analysis
cervix.integrated$major_celltypes[
  cervix.integrated$major_celltypes %in%
    c("CD14_Mono", "CD16_Mono")
] <- "Monocyte"

# Create cell metadata

cell_info <- data.frame(
  condition = colnames(cervix.integrated),
  major_celltypes = cervix.integrated$major_celltypes,
  CT_infection = cervix.integrated$condition,
  stringsAsFactors = FALSE
)

# Add cell metadata to VIPER scores
viper_scores <- acts %>%
  dplyr::left_join(
    cell_info,
    by = "condition"
  )

# Check that cells were matched correctly
cat("Number of VIPER rows:", nrow(viper_scores), "\n")

cat(
  "Rows missing cell-type metadata:",
  sum(is.na(viper_scores$major_celltypes)),
  "\n"
)

# Keep CD4 and CD8 T cells

viper_Tcells <- viper_scores %>%
  dplyr::filter(
    major_celltypes %in% c("CD4_T", "CD8_T"),
    CT_infection %in% c("infected", "uninfected")
  )

# Wilcoxon analysis

wilcox_results <- viper_Tcells %>%
  dplyr::group_by(major_celltypes, source) %>%
  dplyr::summarise(
    n_infected = sum(
      CT_infection == "infected" & !is.na(score)
    ),
    
    n_uninfected = sum(
      CT_infection == "uninfected" & !is.na(score)
    ),
    
    mean_infected = mean(
      score[CT_infection == "infected"],
      na.rm = TRUE
    ),
    
    mean_uninfected = mean(
      score[CT_infection == "uninfected"],
      na.rm = TRUE
    ),
    
    median_infected = median(
      score[CT_infection == "infected"],
      na.rm = TRUE
    ),
    
    median_uninfected = median(
      score[CT_infection == "uninfected"],
      na.rm = TRUE
    ),
    
    diff_mean = mean_infected - mean_uninfected,
    
    diff_median = median_infected - median_uninfected,
    
    wilcox_p = wilcox.test(
      score ~ CT_infection,
      exact = FALSE
    )$p.value,
    
    .groups = "drop"
  )

# BH adjustment

wilcox_results <- wilcox_results %>%
  dplyr::group_by(major_celltypes) %>%
  dplyr::mutate(
    adj_p = p.adjust(
      wilcox_p,
      method = "BH"
    )
  ) %>%
  dplyr::ungroup()

# Keep only positive mean differences

wilcox_results <- wilcox_results %>%
  dplyr::filter(diff_mean > 0)

# Prepare individual sheets

CD4_results <- wilcox_results %>%
  dplyr::filter(major_celltypes == "CD4_T") %>%
  dplyr::arrange(desc(diff_mean)) %>%
  dplyr::rename(
    regulator = source
  ) %>%
  dplyr::select(
    regulator,
    n_infected,
    n_uninfected,
    mean_infected,
    mean_uninfected,
    median_infected,
    median_uninfected,
    diff_mean,
    diff_median,
    wilcox_p,
    adj_p
  )

CD8_results <- wilcox_results %>%
  dplyr::filter(major_celltypes == "CD8_T") %>%
  dplyr::arrange(desc(diff_mean)) %>%
  dplyr::rename(
    regulator = source
  ) %>%
  dplyr::select(
    regulator,
    n_infected,
    n_uninfected,
    mean_infected,
    mean_uninfected,
    median_infected,
    median_uninfected,
    diff_mean,
    diff_median,
    wilcox_p,
    adj_p
  )

# Print retained regulators for verification

cat(
  "\nCD4 T-cell regulators with diff_mean > 0:",
  nrow(CD4_results),
  "\n"
)

print(
  CD4_results %>%
    dplyr::select(
      regulator, diff_mean, diff_median, wilcox_p, adj_p
    )
)

cat(
  "\nCD8 T-cell regulators with diff_mean > 0:",
  nrow(CD8_results),
  "\n"
)

print(
  CD8_results %>%
    dplyr::select(
      regulator, diff_mean, diff_median, wilcox_p, adj_p
    )
)

# Export Excel workbook

dir.create(
  dirname(output_file),
  recursive = TRUE,
  showWarnings = FALSE
)

writexl::write_xlsx(
  list(
    CD4_T = CD4_results,
    CD8_T = CD8_results
  ),
  path = output_file
)

cat("\nExcel workbook written to:\n")
cat(output_file, "\n")
