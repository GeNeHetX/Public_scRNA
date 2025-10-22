library(stringr)
library(purrr)
library(dplyr)
library(Seurat)
library(patchwork)
library(SeuratDisk)
library(readxl)
library(openxlsx)
library(Matrix)

source("/mnt/d/06_SingleCellDataset/FonctionSingleCell.r")


project <- "Lin"
path <- "/mnt/d/06_SingleCellDataset"
rawpath <- file.path(path, project, "01RawData")

CreateDataset(project = project, path = path)

 # === 1. Download data ===
download.file("https://www.dropbox.com/sh/ggn759pr2y36iix/AABZD0nhQ68OM6SG8LDsZbS4a?dl=1",
              file.path(rawpath,"Data_Lin2020_Pancreas.zip"))

unzip(zipfile = file.path(rawpath,"Data_Lin2020_Pancreas.zip"),
      exdir = rawpath, list = F,unzip = getOption("unzip"))
      

# === 2. Clinical and sample annotations ===
Clinic = as.data.frame(data.table::fread(file = file.path(rawpath,"Meta-data.csv")))
names(Clinic)[names(Clinic) == "sample"] <- "oldSampleID"
Clinic$patientID <- sprintf("Lin_P%02d", match(Clinic$patient, unique(Clinic$patient)))
Clinic$sampleID <- sprintf("Lin_S%02d", match(Clinic$patientID, unique(Clinic$patientID)))
Clinic$organism = 'homo_sapiens'
names(Clinic)[names(Clinic) == "Diagnosis"] <- "disease"
names(Clinic)[names(Clinic) == "Primary or Metastasis"] <- "tumorType"
Clinic$disease = 'PDAC'

Clinic$patientSampling <- ifelse(
  startsWith(Clinic$patient, "M"), "biopsy",
  ifelse(startsWith(Clinic$patient, "P"), "surgery", NA)
)

Clinic$specimenConservation = 'fresh'
Clinic$samplePathologicalState = 'tumor'
Clinic$hadTreatment = "NA"
Clinic <- Clinic %>%
  mutate(specimenOrgan = case_when(
    site == "Pancreas"     ~ "pancreas",
    site == "Liver"   ~ "liver",
    site == "Omentum" ~ "omentum"
  )) 

Clinic$tissueUnitExtraction = 'cell'
Clinic$technology = "SingleCell 10x 3' v2"
Clinic$treatmentInfo = "NA"


# === 3. Cells annotations ===
CellsAnnot = as.data.frame(data.table::fread(file = file.path(rawpath,"Cells.csv")))
sample_map  <- setNames(Clinic$sampleID, Clinic$oldSampleID)
names(sample_map) <- sub("^M", "MET", names(sample_map))
CellsAnnot$sampleID  <- sample_map[CellsAnnot$sample]
names(CellsAnnot)[names(CellsAnnot) == "sample"] <- "oldSampleID"

prefix <- sub(":.*", "", CellsAnnot$cell_name)
new_prefix <- sample_map[prefix]

CellsAnnot$cell_name <- paste0(new_prefix, "_", sub(".*:", "", CellsAnnot$cell_name))
names(CellsAnnot)[names(CellsAnnot) == "cell_name"] <- "barcodes"

class_map <- c(
  EMT         = "transition",
  Fibroblast  = "fibroblast",
  Macrophage  = "myeloid",
  Endothelial = "endothelial",
  ETC         = "epithelial_malignant",
  TIL         = "lymphoid_malignant"
)

CellsAnnot$publishedClassL1 <- class_map[CellsAnnot$cell_type]
names(CellsAnnot)[names(CellsAnnot) == "cell_type"] <- "publishedClassL2"

# === 4. Count matrix ===
count <- readMM(file.path(rawpath,"Exp_data_UMIcounts.mtx"))
genes <- readLines(file.path(rawpath,"Genes.txt"))
length(genes) == nrow(count)  
rownames(count) <- genes

# === 5. Export data ===
export_data(
  count = count,
  CellsAnnot = CellsAnnot,
  Clinic = Clinic,
  file_path = file.path(path, project, "03VerifiedDataSet")
)
