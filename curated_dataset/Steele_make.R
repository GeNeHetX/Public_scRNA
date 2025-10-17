library(Matrix)
library(dplyr)
library(readr)
library(stringr)
library(R.utils)
library(readxl)


source("/mnt/d/scpublicdataset/FonctionSingleCell.r")

project <- "Steele"
path <- "/mnt/d/scpublicdataset"
CreateDataset(project, path)

# download CCCA files from :
# https://www.weizmann.ac.il/sites/3CA/pancreas

raw_dir <- "/mnt/d/scpublicdataset/Steele/01RawData/Data_Steele2020_Pancreas"
out_dir <- "/mnt/d/scpublicdataset/Steele/03VerifiedDataSet"
dir.create(out_dir, showWarnings = FALSE)

cell_file   <- file.path(raw_dir, "Cells.csv")
mtx_file    <- file.path(raw_dir, "Exp_data_UMIcounts.mtx")
genes_file  <- file.path(raw_dir, "Genes.txt")
samples_file <- file.path(raw_dir, "Samples.csv")

cells <- read_csv(cell_file) %>%
  mutate(raw_barcode = sub(".*_", "", cell_name))
mtx   <- readMM(mtx_file)
genes <- read_tsv(genes_file, col_names = FALSE)[[1]]
rownames(mtx) <- genes
samples_meta <- read_csv(samples_file)

barcode_list <- list()
mtx_list <- list()
samplesAnnot_list <- list()

unique_samples <- unique(cells$sample)

for (samp in unique_samples) {
  cells_samp <- cells %>% filter(sample == samp)
  cell_indices <- which(cells$sample == samp)
  mtx_samp <- mtx[, cell_indices, drop = FALSE]

  short_id <- gsub("^PDAC_TISSUE_", "", samp)
  if (grepl("^[0-9]+$", short_id)) {
    short_id <- sprintf("%02d", as.integer(short_id))
  }
  sample_short <- paste0("S", short_id)
  sampleID_std <- paste0(project, "_", sample_short)
  patient_short <- paste0("P", short_id)
  patientID_std <- paste0(project, "_", patient_short)
  oldSampleID <- samp

  barcodes <- data.frame(
    barcodes = paste0("Steele_", sample_short, "_", cells_samp$raw_barcode),
    sampleID = sampleID_std,
    oldSampleID = oldSampleID,
    publishedClassL1 = cells_samp$cell_type
  )

  barcode_list[[samp]] <- barcodes
  mtx_list[[samp]] <- as(mtx_samp, "dgCMatrix")

  hadTreatment <- FALSE
  treatmentInfo <- 'naive'
  samplePathologicalState <-"tumor"

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
    technology = "SingleCell 10x",
    disease = "PDAC",
    organism = "homo_sapiens"
  )

  samplesAnnot_list[[samp]] <- clinic_samp
}

cat("Processed samples:", length(samplesAnnot_list), "\n")

all_barcodes <- do.call(rbind, barcode_list)
all_samplesAnnot <- do.call(rbind, samplesAnnot_list)
all_mtx <- do.call(cbind, mtx_list)

check_sample_consistency(df_samples = all_samplesAnnot, df_barcodes = all_barcodes,
                         sample_col = "sampleID", barcode_col = "barcodes")

export_data(count = all_mtx,
            CellsAnnot = all_barcodes,
            Clinic = all_samplesAnnot,
            file_path = out_dir)