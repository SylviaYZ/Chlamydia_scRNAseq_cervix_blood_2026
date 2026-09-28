# quality control

library(dplyr)
library(ggplot2)
library(Seurat)
library(scCustomize)
library(table1)
packageVersion("Seurat")
# Load data
blood.all <- readRDS(file = "/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/PBMC/SCTransformV2_Oct26/blood_all.RDS")

blood.all
# Examine sample quality
table1(~ nCount_RNA + nFeature_RNA | ID, data = blood.all@meta.data,
       overall = FALSE)

VlnPlot(blood.all, features = c("nFeature_RNA"),
        group.by = "ID") + ggtitle("nFeature_RNA by ID")

VlnPlot(blood.all, features = c( "nCount_RNA"),
        group.by = "ID") + ggtitle("nCount_RNA by ID")

table1(~ nCount_RNA + nFeature_RNA | orig.ident, data = blood.all@meta.data)

sample.info.prior.filter <- blood.all@meta.data %>% as.data.frame() %>%
  group_by(ID) %>%
  summarise(Chlamydia_infection = ifelse(unique(condition) == "infected",
                                         "Positive", "Negative"),
            Run = unique(orig.ident),
            Num_cells = n(),
            Median_RNA = median(nCount_RNA),
            Mean_RNA = mean(nCount_RNA),
            Median_GeneExpressed = median(nFeature_RNA),
            Mean_GeneExpressed = mean(nFeature_RNA))

xlsx::write.xlsx(sample.info.prior.filter, file = "sample_info_prior_filter.xlsx", sheetName = "Sheet1",
  col.names = TRUE,  append = FALSE)

# Seurat processing
## QC
### Check MT genes
mtgenes <- grepl("^MT-", rownames(blood.all))
rownames(blood.all)[mtgenes]

blood.all[["percent.mt"]] <- PercentageFeatureSet(blood.all, pattern = "^MT-")

### Check Ribosomal genes
rbgenes <- grepl("^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", rownames(blood.all))
rownames(blood.all)[rbgenes]

blood.all[["percent.rb"]] <- PercentageFeatureSet(blood.all, pattern ="^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")

### Filter cells
VlnPlot(blood.all, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

VlnPlot(blood.all, features = c("percent.mt", "percent.rb"), ncol = 2)

FeatureScatter(blood.all, feature1 = "nCount_RNA", feature2 = "percent.mt")

FeatureScatter(blood.all, feature1 = "percent.rb", feature2 = "percent.mt")

FeatureScatter(blood.all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

hist(blood.all@meta.data$nFeature_RNA, breaks = 100)

hist(blood.all@meta.data$nCount_RNA, breaks = 100)

hist(blood.all@meta.data$percent.mt, breaks = 100)

blood.all <- subset(blood.all,
                         subset = nFeature_RNA > 700 & nFeature_RNA < 2000 &
                           percent.mt < 10 & nCount_RNA > 1000)

blood.all

table1(~ nCount_RNA + nFeature_RNA + percent.mt| ID, data = blood.all@meta.data,
       overall = FALSE)

table1(~ nCount_RNA + nFeature_RNA + percent.mt| orig.ident, data = blood.all@meta.data)

VlnPlot(blood.all, features = c("nFeature_RNA"),
        group.by = "ID") + ggtitle("nFeature_RNA by ID after filtering")

VlnPlot(blood.all, features = c( "nCount_RNA"),
        group.by = "ID") + ggtitle("nCount_RNA by ID after filtering")

sample.info.after.filter <- blood.all@meta.data %>% as.data.frame() %>%
  group_by(ID) %>%
  summarise(Chlamydia_infection = ifelse(unique(condition) == "infected",
                                         "Positive", "Negative"),
            Run = unique(orig.ident),
            Num_cells = n(),
            Median_RNA = median(nCount_RNA),
            Mean_RNA = mean(nCount_RNA),
            Median_GeneExpressed = median(nFeature_RNA),
            Mean_GeneExpressed = mean(nFeature_RNA))

xlsx::write.xlsx(sample.info.after.filter, file = "sample_info_after_filter.xlsx", sheetName = "Sheet1",
  col.names = TRUE,  append = FALSE)

### Save filtered data

saveRDS(blood.all, file = "/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/PBMC/SCTransformV2_Oct26/blood_all_filtered.RDS")

