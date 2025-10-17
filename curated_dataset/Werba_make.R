library(Matrix)
library(dplyr)
library(readr)
library(stringr)
library(R.utils)
library(readxl)
library(GEOquery)

source("/mnt/d/scpublicdataset/FonctionSingleCell.r")

project <- "Werba"
path <- "/mnt/d/scpublicdataset"
CreateDataset(project, path)


# Get clinical data
Clinic=getGEO('GSE205013',GSEMatrix=TRUE)[["GSE205013_series_matrix.txt.gz"]]@phenoData@data
rownames(Clinic) = Clinic$title
Clinic$PatientID = unlist(lapply(str_split(Clinic$title, "_"),"[[",1))


# download data with :
system("wget -c https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205013/suppl/GSE205013_RAW.tar -P /mnt/d/scpublicdataset/Werba/01RawData")
destfile <- file.path(path, project, "01RawData/GSE205013_RAW.tar")
untar(tarfile = destfile,
      exdir = file.path(path, project, "01RawData"), list = FALSE)


raw_dir <- "/mnt/d/scpublicdataset/Werba/01RawData"
out_dir <- "/mnt/d/scpublicdataset/Werba/03VerifiedDataSet"
dir.create(out_dir, showWarnings = FALSE)

barcode_list <- list()
mtx_list <- list()
samplesAnnot_list <- list()


for (pid in Clinic$PatientID) {
  cat("Processing sample:", pid, "\n")
  barcodes_file <- list.files(raw_dir,
                              pattern = paste0(pid, "_barcodes\\.tsv\\.gz$"),
                              full.names = TRUE)
  features_file <- list.files(raw_dir,
                              pattern = paste0(pid, "_features\\.tsv\\.gz$"),
                              full.names = TRUE)
  mtx_file      <- list.files(raw_dir,
                              pattern = paste0(pid, "_matrix\\.mtx\\.gz$"),
                              full.names = TRUE)

  clinic_row <- Clinic[Clinic$PatientID == pid, ]
  barcodes <- read_tsv(barcodes_file, col_names = FALSE,
                       show_col_types = FALSE)[[1]]
  features <- read_tsv(features_file, col_names = FALSE,
                       show_col_types = FALSE)[[2]]
  mtx <- readMM(mtx_file)
  rownames(mtx) <- features
  colnames(mtx) <- barcodes

  short_id <- gsub("^P", "", pid)
  sample_short <- paste0("S", short_id)
  sampleID_std <- paste0("Werba_", sample_short)
  patientID_std <- paste0("Werba_P", short_id)
  oldSampleID <- pid

  barcodes_df <- data.frame(
    barcodes = paste0(sampleID_std, "_", barcodes),
    sampleID = sampleID_std,
    oldSampleID = oldSampleID,
    publishedClassL1 = NA
  )
  barcode_list[[pid]] <- barcodes_df

  mtx_list[[pid]] <- as(mtx, "dgCMatrix")

  patientSampling <- case_when(
    grepl("biopsy", clinic_row$`procedure:ch1`,
          ignore.case = TRUE) ~ "biopsy",
    grepl("resection", clinic_row$`procedure:ch1`,
          ignore.case = TRUE) ~ "surgery",
    TRUE ~ "unknown"
  )

  hadTreatment <- ifelse(grepl("Treated", clinic_row$`treatment:ch1`),
                         TRUE, FALSE)
  treatmentInfo <- ifelse(hadTreatment, "treated", "naive")
  specimenOrgan <- ifelse(clinic_row$source_name_ch1 == "Primary PDAC",
                          "pancreas",
                          "liver")
  clinic_samp <- data.frame(
    sampleID = sampleID_std,
    oldSampleID = oldSampleID,
    patientID = patientID_std,
    patientSampling = patientSampling,
    specimenConservation = "fresh",
    specimenSampling = NA,
    samplePathologicalState = "tumor",
    hadTreatment = hadTreatment,
    treatmentInfo = treatmentInfo,
    specimenOrgan = specimenOrgan,
    tissueUnitExtraction = "cell",
    technology = "SingleCell 10x 3' v3",
    disease = "PDAC",
    organism = "homo_sapiens"
  )
  samplesAnnot_list[[pid]] <- clinic_samp
}
 
all_barcodes <- do.call(rbind, barcode_list)
all_samplesAnnot <- do.call(rbind, samplesAnnot_list)
all_mtx <- do.call(cbind, mtx_list)

check_sample_consistency(df_samples = all_samplesAnnot,
                         df_barcodes = all_barcodes,
                         sample_col = "sampleID", barcode_col = "barcodes")

export_data(count = all_mtx,
            CellsAnnot = all_barcodes,
            Clinic = all_samplesAnnot,
            file_path = out_dir)