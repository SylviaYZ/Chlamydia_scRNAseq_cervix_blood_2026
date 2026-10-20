# Prepare the public range-adjusted cytokine table for Figure 5

library(dplyr)
library(readr)

# Input supplied in data/; no preceding script is needed. Values are already batch-corrected and range-adjusted log2 concentrations.
goodcytokines <- read_csv("data/cytokines_range_adjusted_public.csv", show_col_types = FALSE) %>%
  filter(CT %in% c("CT-", "CT+")) %>%
  select(PTID, CT, `CXCL13 (BCA1) - Good`, `CCL5 (RANTES) - Good`,
         `CXCL9 - Good`, `CXCL10 (IP-10) - Good`, `CXCL11 - Good`)

out_dir <- "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Figures/Fig3_CD4/Cervix_Cytokines/cleandata"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(goodcytokines, file.path(out_dir, "goodcytokines_rangeapplied_KY.csv"))
