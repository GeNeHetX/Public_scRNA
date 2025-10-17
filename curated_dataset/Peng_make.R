library(Matrix)
library(dplyr)
library(readr)
library(stringr)
library(R.utils)
library(readxl)


source("/mnt/d/scpublicdataset/FonctionSingleCell.r")

project <- "Peng"
path <- "/mnt/d/scpublicdataset"
CreateDataset(project,path)
# download count-matrix.txt and all_celltypes.txt files from :
# https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA001063
# ---> https://ngdc.cncb.ac.cn/gsa/browse/CRA001160
#      ---> https://download.cncb.ac.cn/gsa/CRA001160/
# Download files : count-matrix.txt and all_celltype.txt

raw_dir <- "/mnt/d/scpublicdataset/Peng/01RawData"
out_dir <- "/mnt/d/scpublicdataset/Peng/03VerifiedDataSet"
dir.create(out_dir, showWarnings = FALSE)

mtx_file <- file.path(raw_dir, "count-matrix.txt")
celltype_file <- file.path(raw_dir, "all_celltype.txt")


expr <- read.table(mtx_file,
                   header = TRUE,
                   row.names = 1,
                   sep = " ",
                   quote = "\"",
                   check.names = FALSE)
celltypes <- read_tsv(celltype_file, col_names = TRUE)

celltypes <- celltypes %>%
  rename(cell_name = cell.name)

celltypes <- celltypes %>%
  mutate(sample = sub("_.*$", "", cell_name))

barcode_list <- list()
mtx_list <- list()
samplesAnnot_list <- list()

unique_samples <- unique(celltypes$sample)
for (samp in unique_samples) {
  cells_samp <- celltypes %>% filter(sample == samp)
  cell_indices <- which(celltypes$sample == samp)
  mtx_samp <- expr[, cell_indices, drop = FALSE]

  sampleID_std <- paste0(project, "_S", samp)  
  patientID_std <- paste0(project, "_P", samp) 
  oldSampleID <- samp

  raw_barcodes <- sub("^[^_]+_", "", cells_samp$cell_name)
  barcodes <- data.frame(
    barcodes = paste0(sampleID_std, "_", raw_barcodes),
    sampleID = sampleID_std,
    oldSampleID = oldSampleID,
    publishedClassL1 = cells_samp$cluster
  )

  barcode_list[[samp]] <- barcodes
  mtx_list[[samp]] <- mtx_samp

  hadTreatment <- FALSE
  treatmentInfo <- 'naive'
  samplePathologicalState <- ifelse(grepl("N", samp), "normal", "tumor")

  clinic_samp <- data.frame(
    sampleID = sampleID_std,
    oldSampleID = oldSampleID,
    patientID = patientID_std,
    patientSampling = "surgery",
    specimenConservation = "fresh",
    specimenSampling = NA,
    samplePathologicalState = samplePathologicalState,
    hadTreatment = hadTreatment,
    treatmentInfo = treatmentInfo,
    specimenOrgan = "pancreas",
    tissueUnitExtraction = "cell",
    technology = "Single_cell 10x 3' V2",
    disease = "PDAC",
    organism = "homo_sapiens"
  )

  samplesAnnot_list[[samp]] <- clinic_samp
}

cat("Processed samples:", length(samplesAnnot_list), "\n")

all_barcodes <- do.call(rbind, barcode_list)
all_samplesAnnot <- do.call(rbind, samplesAnnot_list)
all_mtx <- do.call(cbind, mtx_list)
all_mtx_sparse <- as(as.matrix(all_mtx), "dgCMatrix")

check_sample_consistency(df_samples = all_samplesAnnot, df_barcodes = all_barcodes,
                         sample_col = "sampleID", barcode_col = "barcodes")

export_data(count = all_mtx_sparse,
            CellsAnnot = all_barcodes,
            Clinic = all_samplesAnnot,
            file_path = out_dir)