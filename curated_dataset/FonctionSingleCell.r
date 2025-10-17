# Liste des fonctions pour valider l'exportation des données

library(Matrix)
library(dplyr)
library(readr)
library(utils)


# Fonction pour créer un dataset
CreateDataset <- function(project, path) {
  project_path <- file.path(path, project)
  
  dir.create(project_path, recursive = TRUE, showWarnings = FALSE)
  
  subfolders <- c("01RawData", "02Processes", "03VerifiedDataSet")
  sapply(file.path(project_path, subfolders), dir.create, showWarnings = FALSE)
  
  cat("Le dossier", project, "et ses sous-dossiers ont été créés avec succès dans :", path, "\n")
}



# Fonction pour vérifier si les samples sont uniques et correspondent aux barcodes
check_sample_consistency <- function(df_samples, df_barcodes,
                                     sample_col = "sampleID",
                                     barcode_col = "barcodes") {
  if (!all(c(sample_col, barcode_col) %in% c(names(df_samples), names(df_barcodes)))) {
    stop("Les colonnes spécifiées n'existent pas dans les dataframes.")
  }

  prefixes <- sub("^(([^_]+)_([^_]+)).*", "\\1", df_barcodes[[barcode_col]])

  valid_samples <- df_samples[[sample_col]]
  invalid <- !prefixes %in% valid_samples

  if (any(invalid)) {
    invalid_prefixes <- unique(prefixes[invalid])
    stop(paste0("Les barcodes suivants ne correspondent à aucun sample : ",
                paste(invalid_prefixes, collapse = ", ")))
  }

  return(TRUE)
}

# Fonction pour exporter les données dans plusieurs formats 

export_data <- function(count, CellsAnnot, Clinic, file_path) {
  # Création du dossier de sortie
  dir.create(file_path, showWarnings = FALSE, recursive = TRUE)

  # === 1. matrix.mtx.gz ===
  raw_mtx_path <- file.path(file_path, "matrix.mtx")
  Matrix::writeMM(count, file = raw_mtx_path)

  gz_mtx_path <- file.path(file_path, "matrix.mtx.gz")
  R.utils::gzip(raw_mtx_path, destname = gz_mtx_path, overwrite = TRUE)

  # === 2. features.tsv.gz ===
  features <- data.frame(geneID = rownames(count))
  write_tsv(features, file = gzfile(file.path(file_path, "features.tsv.gz")))

  # === 3. barcodes.tsv.gz ===
  if (!all(c("barcodes", "sampleID", "oldSampleID") %in% colnames(CellsAnnot))) {
    stop("Les colonnes 'barcodes' et/ou 'sampleID' sont manquantes dans CellsAnnot.")
}

  cols_to_write <- c("barcodes", "sampleID", "oldSampleID")

  if ("publishedClassL1" %in% colnames(CellsAnnot)) {
      cols_to_write <- c(cols_to_write, "publishedClassL1")
  }
  if ("publishedClassL2" %in% colnames(CellsAnnot)) {
      cols_to_write <- c(cols_to_write, "publishedClassL2")
  }
  options(check.names = FALSE)
  write_tsv(CellsAnnot[, cols_to_write, drop = FALSE],
            file = gzfile(file.path(file_path, "barcodes.tsv.gz")))


  # === 4. clinical.Annot.tsv ===
  required_clinic_cols <- c("patientID", "disease", "organism")
  missing_clinic_cols <- setdiff(required_clinic_cols, colnames(Clinic))
  if (length(missing_clinic_cols) > 0) {
    stop("Colonnes manquantes dans Clinic : ", paste(missing_clinic_cols, collapse = ", "))
  }

  clinical <- Clinic %>%
    select(all_of(required_clinic_cols)) %>%
    distinct()
  write_tsv(clinical, file = file.path(file_path, "clinicalAnnot.tsv"))

  # === 5. sample.Annot.tsv ===
  required_sample_cols <- c("sampleID","oldSampleID","patientID","patientSampling","specimenConservation", "samplePathologicalState", "hadTreatment","treatmentInfo", "specimenOrgan","tissueUnitExtraction","technology")
  missing_sample_cols <- setdiff(required_sample_cols, colnames(Clinic))
  if (length(missing_sample_cols) > 0) {
    stop("Colonnes manquantes dans Clinic pour samplesAnnot.tsv : ",
         paste(missing_sample_cols, collapse = ", "))
  }

  sample <- Clinic %>%
    select(all_of(required_sample_cols)) %>%
    distinct()
  write_tsv(sample, file = file.path(file_path, "samplesAnnot.tsv"))

  message("Export terminé avec succès dans le dossier : ", file_path)

}


