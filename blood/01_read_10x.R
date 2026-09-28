# read 10x

library(dplyr)
library(ggplot2)
library(Seurat)
library(scCustomize)
library(table1)
packageVersion("Seurat")
# Load original data from 10X
data.dir <- "/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/PBMC/data"

run <- list.files(data.dir)
run

blood <- list()

# Function to convert 10X object to Seurat object
get10X <- function(folder.name){

  exprs.list <- Read10X(data.dir = paste0(data.dir, "/",folder.name))

  Antibody <- exprs.list$`Antibody Capture`
  RNA <- exprs.list$`Gene Expression`

  id <- apply(Antibody, 2, function(x){
    rownames(Antibody)[which.max(x)]
  })

  print(paste0("Number of cells from each participant:"))
  print(table(id))

  ID <- data.frame(ID = id)

  Seurat.obj <- CreateSeuratObject(counts = RNA,
                             project = folder.name,
                             assay= "RNA",
                             meta.data = ID,
                             min.cells = 0,
                             min.features = 0)

  Seurat.obj <- RenameCells(Seurat.obj,
                            # prefix to add cell names
                            add.cell.id = folder.name
                            )

  print(Seurat.obj)

  return(Seurat.obj)
}

# Run: aggr_p1
blood[["aggr_p1"]] <- get10X("aggr_p1")
# Run: p2
blood2.1 <- get10X("p2_1")
blood2.2 <- get10X("p2_2")

# One participant id was incorrect.
blood2.2$ID[which(blood2.2$ID == 3170)] <- 3015

blood[["p2"]] <- merge(x=blood2.1, y=blood2.2,
                       project = "p2")
# Merge all runs together
blood.all <- Merge_Seurat_List(blood)
blood.all

# 11 participants
length(unique(blood.all$ID))

table1(~ ID, data = blood.all@meta.data)
# Define infection status
blood.all$condition <- "infected"
blood.all$condition[blood.all$ID %in% c("3044", "3150", "3070")] <- "uninfected"

table1(~ condition, data = blood.all@meta.data)

# Save data

saveRDS(blood.all, file = "/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/PBMC/SCTransformV2_Oct26/blood_all.RDS")

