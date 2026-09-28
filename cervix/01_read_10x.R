# read 10x

library(dplyr)
library(ggplot2)
library(Seurat)
library(scCustomize)
library(table1)
packageVersion("Seurat")
# Load original data from 10X
data.dir <- "/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/data"

run <- list.files(data.dir)
run

cervix <- list()

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
                             min.features = 0,)

  Seurat.obj <- RenameCells(Seurat.obj,
                            # prefix to add cell names
                            add.cell.id = folder.name
                            )

  print(Seurat.obj)

  return(Seurat.obj)
}

# Run: aggr_cx1
cervix[["aggr_cx1"]] <- get10X("aggr_cx1")
# Run: cx2
cervix[["cx2"]] <- get10X("cx2")
# Run: cx3
cx3.1 <- get10X("cx3_1")
cx3.2 <- get10X("cx3_2")

cervix[["cx3"]] <- merge(x=cx3.1, y=cx3.2,
                         project = "cx3")
# Run: aggr_cx4
cervix[["aggr_cx4"]] <- get10X("aggr_cx4")
# Run: aggr_cx5
cervix[["aggr_cx5"]] <- get10X("aggr_cx5")
# Run: cx6_1
cervix[["cx6_1"]] <- get10X("cx6_1")
# Merge all runs together
cervix.all <- Merge_Seurat_List(cervix)
cervix.all

# 22 participants
length(unique(cervix.all$ID))

table1(~ ID, data = cervix.all@meta.data)
# Define infection status
cervix.all$condition <- "infected"
cervix.all$condition[cervix.all$ID %in% c("3044", "3038", "3159", "3213")] <- "uninfected"

table1(~ condition, data = cervix.all@meta.data)

# Save data

saveRDS(cervix.all, file = "/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/Cervix/SCTransformV2_Oct26/cervix_all.RDS")

