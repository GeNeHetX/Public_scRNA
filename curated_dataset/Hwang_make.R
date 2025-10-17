library(Matrix)
library(dplyr)
library(readr)
library(stringr)
library(R.utils)

source("/mnt/d/scpublicdataset/FonctionSingleCell.r")

project <- "Hwang"
path <- "/mnt/d/scpublicdataset"
CreateDataset(project,path)
# download CCCA files from : 
# https://www.weizmann.ac.il/sites/3CA/pancreas

raw_dir <- "/mnt/d/scpublicdataset/Hwang/01RawData/Data_Hwang2022_Pancreas"
out_dir <- "/mnt/d/scpublicdataset/Hwang/03VerifiedDataSet"
dir.create(out_dir, showWarnings = FALSE)

genes <- read_tsv(file.path(raw_dir, "genes.txt"), col_names = FALSE)
gene_names <- genes$X1

barcode_list <- list()
mtx_list <- list()
samplesAnnot_list <- list()

groups <- list.dirs(raw_dir, recursive = FALSE, full.names = TRUE)
for (grp in groups) {
  grp_name <- basename(grp)
  cell_file <- list.files(grp, pattern = "^Cells.*\\.csv$", full.names = TRUE)
  mtx_file <- list.files(grp, pattern = "^Exp_data_TP10K.*\\.mtx$", full.names = TRUE)
  cells <- read_csv(cell_file)
  mtx <- readMM(mtx_file)
  rownames(mtx) <- gene_names
  unique_samples <- unique(cells$sample)
  
  for (samp in unique_samples) {
    cells_samp <- cells %>% filter(sample == samp)
    cell_indices <- which(cells$sample == samp)
    mtx_samp <- mtx[, cell_indices, drop = FALSE]
    sampleID_std <- paste0("S", samp)   # new sampleID
    sampleID_std <- paste0(project, "_", sampleID_std)
    patientID_std <- paste0("P", samp)
    patientID_std <- paste0(project, "_", patientID_std)
    oldSampleID <- samp                  # original sample name
    barcodes <- data.frame(
      barcodes = paste0(sampleID_std, "_", cells_samp$cell_name),
      sampleID = sampleID_std,
      oldSampleID = oldSampleID,
      publishedClassL1 = cells_samp$cell_type
    )
    if ("cell_subtype" %in% colnames(cells_samp)) {
      barcodes$publishedClassL2 <- cells_samp$cell_subtype
    }
    barcode_list[[paste0(grp_name, "_", samp)]] <- barcodes
    mtx_list[[paste0(grp_name, "_", samp)]] <- mtx_samp
    hadTreatment <- ifelse(grepl("T", samp), TRUE, FALSE)
    treatmentInfo <- ifelse(hadTreatment, "neoadjuvant", "naive")
    clinic_samp <- cells_samp %>%
      select(sampleID = sample,
            oldSampleID = sample) %>%
      mutate(
        sampleID = sampleID_std,
        oldSampleID = oldSampleID,
        patientID = patientID_std,
        patientSampling = "surgery",
        specimenConservation = "frozen",
        specimenSampling = NA,
        samplePathologicalState = "tumor",
        hadTreatment = hadTreatment,
        treatmentInfo = treatmentInfo,
        specimenOrgan = "pancreas",
        tissueUnitExtraction = "nucleus",
        technology = "Single_cell 10x 3' V2 or V3",
        disease = "PDAC",
        organism = "homo_sapiens"
      ) %>%
      distinct()
    samplesAnnot_list[[paste0(grp_name, "_", samp)]] <- clinic_samp
  }

}

matrix_all <- do.call(cbind, mtx_list)
barcodes_all <- do.call(rbind, barcode_list)
Clinic <- do.call(rbind, samplesAnnot_list)

check_sample_consistency(df_samples = Clinic, df_barcodes = barcodes_all,
                         sample_col = "sampleID", barcode_col = "barcodes")

export_data(count = matrix_all,
            CellsAnnot = barcodes_all,
            Clinic = Clinic,
            file_path = out_dir)


# treatmentInfo == "None" ~ "naive",
# treatmentInfo == "CRT" ~ "FOLFIRINOX + radiotherapy with concurrent capecitabine or 5-FU",
# treatmentInfo == "CRTL**" ~ "FOLFIRINOX + losartan + radiotherapy with concurrent capecitabine or 5-FU/NCT01821729",
# treatmentInfo == "CRTN*" ~ "FOLFIRINOX + stereotactic body radiotherapy + nivolumab/NCT03563248",
# treatmentInfo == "CRTLN*" ~ "FOLFIRINOX + stereotactic body radiotherapy + losartan + nivolumab/NCT03563248",
# treatmentInfo == "Other" ~ "treatment regimen consisting of chemotherapy and/or radiotherapy combination not otherwise specified",

