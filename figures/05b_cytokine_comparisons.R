#  cytokine comparisons

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)
library(stringr)

goodcytokines <- read.csv("/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Figures/Fig3_CD4/Cervix_Cytokines/cleandata/goodcytokines_rangeapplied_KY.csv")

# analytes to keep
analytes_keep <- c("CXCL13", "CXCL9", "CXCL10", "CXCL11", "CCL5")

# all cytokine columns
cytokine_cols <- setdiff(names(goodcytokines), c("PTID", "CT"))

# reshape to long and keep selected analytes only
goodcytokines_long <- goodcytokines %>%
  mutate(CT = factor(CT, levels = c("CT-", "CT+"))) %>%
  pivot_longer(
    cols = all_of(cytokine_cols),
    names_to = "Analyte",
    values_to = "log2_value"
  ) %>%
  filter(!is.na(CT)) %>%   # filter out 2 PTIDs with equivocal CT diagnostic
  mutate(
    Analyte_clean = str_remove(Analyte, "\\.\\.\\..*$"),
    Analyte_clean = case_when(
      str_detect(Analyte, regex("BCA1|CXCL13", ignore_case = TRUE)) ~ "CXCL13 (BCA1)",
      str_detect(Analyte, regex("RANTES|CCL5", ignore_case = TRUE)) ~ "CCL5 (RANTES)",
      str_detect(Analyte, regex("CXCL10|IP[- ]?10", ignore_case = TRUE)) ~ "CXCL10 (IP-10)",
      str_detect(Analyte, regex("CXCL9|MIG", ignore_case = TRUE)) ~ "CXCL9",
      str_detect(Analyte, regex("CXCL11|I[- ]?TAC", ignore_case = TRUE)) ~ "CXCL11",
      TRUE ~ Analyte_clean
    ),
    Analyte_base = case_when(
      str_detect(Analyte_clean, "^CXCL13") ~ "CXCL13",
      str_detect(Analyte_clean, "^CXCL9$") ~ "CXCL9",
      str_detect(Analyte_clean, "^CXCL10") ~ "CXCL10",
      str_detect(Analyte_clean, "^CXCL11") ~ "CXCL11",
      str_detect(Analyte_clean, "^CCL5") ~ "CCL5",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Analyte_base)) %>%
  distinct(PTID, CT, Analyte, log2_value, Analyte_clean, Analyte_base)

# make one figure per analyte
for (a in analytes_keep) {

  plot_dat <- goodcytokines_long %>%
    filter(Analyte_base == a)

  p <- ggplot(plot_dat, aes(x = CT, y = log2_value, color = CT)) +
    geom_boxplot(
      outlier.shape = NA,
      na.rm = TRUE,
      color = "black",
      width = 0.55
    ) +
    geom_beeswarm(
      na.rm = TRUE,
      cex = 2,
      size = 0.5
    ) +
    stat_compare_means(
      na.rm = TRUE,
      method = "wilcox.test",
      label = "p.signif",
      hide.ns = TRUE,
      label.x.npc = 0.5,
      size = 3
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.15))
    ) +
    scale_color_brewer(palette = "Dark2") +
    labs(
      title = unique(plot_dat$Analyte_clean)[1],
      x = NULL,
      y = paste0("Log2 ", unique(plot_dat$Analyte_clean)[1], " pg/mL")
    ) +
    theme_minimal(base_size = 8) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 8),
      axis.text.x = element_text(size = 7, color = "black"),
      axis.text.y = element_text(size = 7, color = "black"),
      axis.title.y = element_text(size = 8),
      axis.line = element_line(color = "black", linewidth = 0.3),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length = unit(1.5, "mm"),
      panel.grid = element_blank()
    )

  ggsave(
    filename = file.path(
      "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Figures/Fig3_CD4",
      paste0(a, "_CT-vsCT+.svg")
    ),
    plot = p,
    width = 1.5,
    height = 2
  )
}
