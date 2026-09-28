# quality control

library(dplyr)
library(ggplot2)
library(Seurat)
library(scCustomize)
library(table1)
packageVersion("Seurat")
# Load data
cervix.all <- readRDS(file = "/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_all.RDS")

cervix.all
# Remove samples from participants co-infected with TV
cervix.all <- cervix.all[, !cervix.all$ID %in% c(3121, 3157)]
# Examine sample quality
table1(~ nCount_RNA + nFeature_RNA | ID, data = cervix.all@meta.data,
       overall = FALSE)

VlnPlot(cervix.all, features = c("nFeature_RNA"),
        group.by = "ID") + ggtitle("nFeature_RNA by ID")

VlnPlot(cervix.all, features = c( "nCount_RNA"),
        group.by = "ID") + ggtitle("nCount_RNA by ID")

table1(~ nCount_RNA + nFeature_RNA | orig.ident, data = cervix.all@meta.data)

sample.info.prior.filter <- cervix.all@meta.data %>% as.data.frame() %>%
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
mtgenes <- grepl("^MT-", rownames(cervix.all))
rownames(cervix.all)[mtgenes]

cervix.all[["percent.mt"]] <- PercentageFeatureSet(cervix.all, pattern = "^MT-")

### Check Ribosomal genes
rbgenes <- grepl("^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", rownames(cervix.all))
rownames(cervix.all)[rbgenes]

cervix.all[["percent.rb"]] <- PercentageFeatureSet(cervix.all, pattern ="^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")

### Filter cells
VlnPlot(cervix.all, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

VlnPlot(cervix.all, features = c("percent.mt", "percent.rb"), ncol = 2)

FeatureScatter(cervix.all, feature1 = "nCount_RNA", feature2 = "percent.mt")

FeatureScatter(cervix.all, feature1 = "percent.rb", feature2 = "percent.mt")

FeatureScatter(cervix.all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

cervix.all <- subset(cervix.all,
                         subset = nFeature_RNA > 200 & nFeature_RNA < 3000 &
                           percent.mt < 15 & nCount_RNA > 250)

cervix.all

table1(~ nCount_RNA + nFeature_RNA + percent.mt| ID, data = cervix.all@meta.data,
       overall = FALSE)

table1(~ nCount_RNA + nFeature_RNA + percent.mt| orig.ident, data = cervix.all@meta.data)

VlnPlot(cervix.all, features = c("nFeature_RNA"),
        group.by = "ID") + ggtitle("nFeature_RNA by ID after filtering")

VlnPlot(cervix.all, features = c( "nCount_RNA"),
        group.by = "ID") + ggtitle("nCount_RNA by ID after filtering")

sample.info.after.filter <- cervix.all@meta.data %>% as.data.frame() %>%
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
### Filter samples
cervix.all <- cervix.all[, !cervix.all$ID %in% c("3027", "3038", 
                                                 "3040", "3044",
                                                 "3194", "3226",
                                                 "3013")]

table(cervix.all$ID)
### Save filtered data

saveRDS(cervix.all, file = "/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_all_filtered.RDS")

