#  prepare cytokine table

#### Prep ####

library(tidyverse)
library(readxl)

## Import data
cytokines <- read_csv("data/all_samples_batch_correct_clean_cytokine.csv")
key <- read_excel("data/analyte_key.xlsx")
outcomes <- read_excel("data/TRAC1 Summary of Findings.xlsx",
                       na = "NA", skip = 7, n_max = 247)

## Data prep
key <- key %>%
  mutate(Label = paste(Protein, IsGood, sep = " - ")) %>%
  mutate(LowLim = log2(LowLim),
         UpLim = log2(UpLim))

rename_map <- setNames(key$Analyte, key$Label)

low <- setNames(key$LowLim, key$Label)
up <- setNames(key$UpLim,  key$Label)

outcomes <- outcomes %>%
  rename(PTID = 'Patient ID',
         CT = 'E_Cx_CT') %>%
  select(PTID, CT) %>%
  mutate(CT = if_else( # recode false-negatives as CT+ (pos by PCR)
    PTID %in% c("3-101", "3-144", "3-145", "3-148",
                "3-201", "3-209", "3-226"),
    -1,
    CT
  )) %>%
  mutate(CT = factor(
    case_when(
      CT == 0 ~ "CT-",
      CT == -1 ~ "CT+",
      CT == -2 ~ NA_character_
    )
  ))

cytokines <- cytokines %>%
  rename(PTID = id)

## Keep Good and Maybe cytokines, rename using label from key
goodcytokines <- cytokines %>%
  left_join(outcomes, by = 'PTID') %>%
  select(PTID, CT, all_of(key$Analyte)) %>%
  rename(any_of(rename_map))

## How many values will be replaced outside LowLim and UpLim
cap_counts_detail <- sapply(
  intersect(key$Label, names(goodcytokines)),
  \(nm) {
    x <- goodcytokines[[nm]]
    c(
      below = sum(x < low[nm], na.rm = TRUE),
      above = sum(x > up[nm],  na.rm = TRUE)
    )
  }
)

t(cap_counts_detail)
write_csv(
  as_tibble(t(cap_counts_detail), rownames = "Analyte"),
  "cap_counts_detail_KY.csv"
)

## Run replacement
goodcytokines <- goodcytokines %>%
  mutate(
    across(
      any_of(key$Label),
      ~ case_when(
        .x < low[cur_column()] ~ low[cur_column()],
        .x > up[cur_column()]  ~ up[cur_column()],
        TRUE                   ~ .x
      )
    )
  )

dir.create("cleandata", recursive = TRUE, showWarnings = FALSE)
write_csv(goodcytokines, "cleandata/goodcytokines_rangeapplied_KY.csv")
