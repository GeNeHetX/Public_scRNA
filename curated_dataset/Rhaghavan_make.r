library(stringr)
library(purrr)
library(dplyr)
library(Matrix)
library(GEOquery)
library(Matrix)
library(Seurat)
library(readxl)
library(vroom)


source("/mnt/d/06_SingleCellDataset/FonctionSingleCell.r")

# Data accessibility
    # download CCCA files from : 
        # https://www.weizmann.ac.il/sites/3CA/pancreas
    # download Biopsy files from :
        # https://singlecell.broadinstitute.org/single_cell/study/SCP1644/microenvironment-drives-cell-state-plasticity-and-drug-response-in-pancreatic-cancer#study-visualize


project <- "Raghavan"
path <- "/mnt/d/06_SingleCellDataset"
rawpath <- file.path(path, project, "01RawData")

CreateDataset(project = project, path = path)

list.files(rawpath)
CellsAnnot = as.data.frame(data.table::fread(file.path(rawpath,"complete_MetaData_70170cells_scp.csv")))
CellsAnnot = CellsAnnot[-1,]
CellsAnnot = subset(CellsAnnot, donor_ID!="PANFR0580") #PAN NET
CellsAnnot = subset(CellsAnnot, donor_ID!="PANFR0489R") #Ne correspond a aucun PatientID

CellsAnnot = subset(CellsAnnot, sample.type == "Biopsy")

Clinic = data.table::fread(file.path(rawpath,"Data.clinic.csv"))

count1 = as.data.frame(read.csv(file.path(rawpath, "Biopsy_RawDGE_23042cells.csv")))
rownames(count1) = count1$X
count2 = as.data.frame(read.csv(file.path(rawpath, "Biopsy473_RawDGE_1370cells.csv")))
rownames(count2) = count2$X

counts <- Reduce(function(x, y) full_join(x, y, by = "X"), list(count1, count2))
rownames(counts) <- counts$X
counts_filtered <- counts[, colnames(counts) %in% CellsAnnot$NAME]
counts_filtered[is.na(counts_filtered)] <- 0


Clinic$patientID <- sprintf("Ragh_P%02d", match(Clinic$PatientID, unique(Clinic$PatientID)))
Clinic$sampleID <- sprintf("Ragh_S%02d", match(Clinic$patientID, unique(Clinic$patientID)))
patient_map <- setNames(Clinic$patientID, Clinic$PatientID)
sample_map  <- setNames(Clinic$sampleID, Clinic$patientID)

CellsAnnot$patientID <- patient_map[CellsAnnot$donor_ID]
CellsAnnot$sampleID  <- sample_map[CellsAnnot$patientID]
names(CellsAnnot)[names(CellsAnnot) == "Coarse_Cell_Annotations"] <- "publishedClassL1"

patient_code <- sub("^[^_]+_([^_]+)_.*", "\\1", CellsAnnot$NAME)

id_map <- setNames(CellsAnnot$sampleID, CellsAnnot$donor_ID)
CellsAnnot$NAME <- paste0(id_map[patient_code], "_", sub(".*_", "", CellsAnnot$NAME))

old_names <- colnames(counts_filtered)
patient_codes <- sub("^[^_]+_([^_]+)_.*", "\\1", old_names)
suffix <- sub(".*_", "", old_names)
new_names <- paste0(id_map[patient_codes], "_", suffix)
colnames(counts_filtered) <- new_names
counts_filtered <- counts_filtered[, colnames(counts_filtered) %in% CellsAnnot$NAME]

CellsAnnot = data.frame("barcodes" = colnames(counts_filtered),
                        "sampleID" = CellsAnnot$sampleID,
                        "publishedClassL1" = CellsAnnot$publishedClassL1)

View(CellsAnnot)
Clinic$organism = 'homo_sapiens'
names(Clinic)[names(Clinic) == "Histology"] <- "disease"
Clinic$patientSampling = 'biopsy'
Clinic$specimenConservation = 'fresh'
Clinic$samplePathologicalState = 'tumor'
Clinic$hadTreatment <- ifelse(
  tolower(Clinic$Treatment.of.primary.disease) == "none",
  "FALSE",
  "TRUE"
)
Clinic$Site.of.biopsy <- tolower(Clinic$Site.of.biopsy)
names(Clinic)[names(Clinic) == "Site.of.biopsy"] <- "specimenOrgan"
Clinic$tissueUnitExtraction = 'cell'
Clinic$technology = "Seq-Well Library Nextera XT DNA tagmentation"
Clinic$Metastatic.treatments.prior.to.biopsy <- gsub(";", "/", Clinic$Metastatic.treatments.prior.to.biopsy)

Clinic$Metastatic.treatments.prior.to.biopsy <- ifelse(
  tolower(Clinic$Metastatic.treatments.prior.to.biopsy) == "none",
  "naive",
  Clinic$Metastatic.treatments.prior.to.biopsy
)
names(Clinic)[names(Clinic) == "Metastatic.treatments.prior.to.biopsy"] <- "treatmentInfo"
View(Clinic)

counts_mat <- as.matrix(counts_filtered)
storage.mode(counts_mat) <- "integer" 

mat_sparse <- Matrix(counts_mat, sparse = TRUE)

rownames(mat_sparse) <- rownames(counts_mat)
colnames(mat_sparse) <- colnames(counts_mat)

View(CellsAnnot)


export_data(
  count = mat_sparse,
  CellsAnnot = CellsAnnot,
  Clinic = Clinic,
  file_path = file.path(path, project, "03VerifiedDataSet")
)

