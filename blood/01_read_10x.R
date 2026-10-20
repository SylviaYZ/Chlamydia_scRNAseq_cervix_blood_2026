# Read GEO H5 matrices and restore cohort IDs

# Create output folders after updating the paths for your computer.
dir.create("/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/PBMC/SCTransformV2_Oct26", recursive = TRUE, showWarnings = FALSE)

library(dplyr)
library(ggplot2)
library(Seurat)
library(readxl)
library(scCustomize)
library(table1)
packageVersion("Seurat")
# Download the blood H5 matrices from GEO GSE334218 into data.dir.
# The filename prefix matches the GEO library name used in get10X() below.
# Keep the filenames shown below, or adjust the paths if GEO adds accession prefixes.
data.dir <- "/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/PBMC/data_masked_h5"

# Download Supplemental Table 1 from the manuscript supplementary materials and save as data/Supplemental Table 1.xlsx.
# Row 3 contains headers: Assigned ID = GEO masked_ID; Participant ID = cohort ID.
supp_table_file <- "data/Supplemental Table 1.xlsx"
id_map <- readxl::read_excel(supp_table_file, sheet = "Supp. Table 1",
                            skip = 2, col_types = "text") %>%
  filter(Tissue == "Blood") %>%
  transmute(masked_ID = trimws(`Assigned ID`), ID = trimws(`Participant ID`)) %>%
  distinct()

# Each masked ID must have one cohort ID.
stopifnot(nrow(id_map) > 0, !anyNA(id_map),
          all(nzchar(id_map$masked_ID)), all(nzchar(id_map$ID)),
          !anyDuplicated(id_map$masked_ID))

run <- list.files(data.dir, pattern = "_masked_feature_bc_matrix\\.h5$")
run

blood <- list()

# Read one GEO library and assign its masked participant IDs
get10X <- function(folder.name){

  # Input: the GEO H5 file for this library; RNA counts and cell barcodes are unchanged.
  h5_file <- file.path(data.dir, paste0(folder.name, "_masked_feature_bc_matrix.h5"))
  exprs.list <- Read10X_h5(filename = h5_file)
  stopifnot(all(c("Gene Expression", "Antibody Capture") %in% names(exprs.list)))

  Antibody <- exprs.list$`Antibody Capture`
  RNA <- exprs.list$`Gene Expression`

  # Some Seurat versions return a vector when the H5 has only one antibody feature.
  if (is.null(dim(Antibody))) {
    h5 <- hdf5r::H5File$new(h5_file, mode = "r")
    feature_names <- h5[["matrix/features/name"]][]
    feature_types <- h5[["matrix/features/feature_type"]][]
    h5$close_all()
    Antibody <- matrix(Antibody, nrow = 1,
                       dimnames = list(feature_names[feature_types == "Antibody Capture"],
                                       colnames(RNA)))
  }

  masked_id <- apply(Antibody, 2, function(x){
    rownames(Antibody)[which.max(x)]
  })

  print(paste0("Number of cells from each participant:"))
  print(table(masked_id))

  # Keep the GEO masked_ID and add its cohort ID from Supplemental Table 1.
  cell_metadata <- data.frame(masked_ID = masked_id,
                              ID = id_map$ID[match(masked_id, id_map$masked_ID)],
                              row.names = colnames(Antibody))
  if (anyNA(cell_metadata$ID)) {
    stop("Masked IDs missing from Supplemental Table 1: ",
         paste(unique(masked_id[is.na(cell_metadata$ID)]), collapse = ", "))
  }

  Seurat.obj <- CreateSeuratObject(counts = RNA,
                             project = folder.name,
                             assay= "RNA",
                             meta.data = cell_metadata,
                             min.cells = 0,
                             min.features = 0)

  Seurat.obj <- RenameCells(Seurat.obj,
                            # prefix to add cell names
                            add.cell.id = folder.name
                            )

  print(Seurat.obj)

  return(Seurat.obj)
}

# GEO library aggr_p1: aggr_p1_masked_feature_bc_matrix.h5.
blood[["aggr_p1"]] <- get10X("aggr_p1")
# Merge GEO libraries p2_1 and p2_2 into the p2 analysis group.
# GEO library p2_1: p2_1_masked_feature_bc_matrix.h5.
blood2.1 <- get10X("p2_1")
# GEO library p2_2: p2_2_masked_feature_bc_matrix.h5.
blood2.2 <- get10X("p2_2")

# Correct the original p2_2 barcode error after mapping to cohort ID; keep masked_ID unchanged.
# If the H5 already maps these cells to 3015, this correction makes no change.
blood2.2$ID[which(blood2.2$ID == "3170")] <- "3015"

blood[["p2"]] <- merge(x=blood2.1, y=blood2.2,
                       project = "p2")
# Merge all runs together
blood.all <- Merge_Seurat_List(blood)
blood.all

# Check participant counts
length(unique(blood.all$ID))

table1(~ ID, data = blood.all@meta.data)
# Define infection status
blood.all$condition <- "infected"
blood.all$condition[blood.all$ID %in% c("3044", "3150", "3070")] <- "uninfected"

table1(~ condition, data = blood.all@meta.data)

# Save data

saveRDS(blood.all, file = "/Chlamydia-SingleCell/CR 6.1.2 Reanalysis/PBMC/SCTransformV2_Oct26/blood_all.RDS")

