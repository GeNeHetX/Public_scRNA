library(stringr)
library(purrr)
library(dplyr)
library(Seurat)
library(patchwork)
library(SeuratDisk)
library(readxl)
library(Matrix)
library(GEOquery)
library(vroom)


source("/mnt/d/06_SingleCellDataset/FonctionSingleCell.r")


project <- "Lander-DeNardo"
path <- "/mnt/d/06_SingleCellDataset"
rawpath <- file.path(path, project, "01RawData")

CreateDataset(project = project, path = path)

Clinic=getGEO('GSE207536',GSEMatrix=TRUE)[["GSE207536_series_matrix.txt.gz"]]@phenoData@data

download.file(
  url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE207536&format=file",
  destfile = file.path(rawpath, "GSE207536_RAW.tar"),
  mode = "wb",
  method = "curl"
)

untar(file.path(rawpath,"GSE207536_RAW.tar"),
     exdir = rawpath)

Clinic$oldPatientID=gsub("Patient ","Pt", str_extract(Clinic$title, "Patient \\d+"))
Clinic$patientID <- setNames(
  sprintf("BelleDN_P%02d", seq_along(unique(Clinic$oldPatientID))),
  unique(Clinic$oldPatientID)
)[Clinic$oldPatientID]
Clinic$oldSampleID=gsub("Patient ","Pt", gsub("-treatment FNA PDAC tissue biopsy|-treatment FNA PDAC tissue biopsy ","",Clinic$title))
Clinic$oldSampleID=gsub(" ","_",gsub("PRE","Pre",gsub("POST","Post",Clinic$oldSampleID)))
Clinic$sampleID <- {
  ids <- sub("_.*", "", Clinic$oldSampleID)        
  suffix <- sub(".*_", "", Clinic$oldSampleID) 
  mapping <- setNames(sprintf("LanderDN_S%02d", seq_along(unique(ids))),unique(ids))
  paste0(mapping[ids], suffix)
}
Correspondance <- as.data.frame(unique(Clinic[, c("oldSampleID", "sampleID")]), row.names = NULL)


files <- list.files(rawpath, pattern = "\\.csv\\.gz$", full.names = TRUE)
Sys.setenv("VROOM_CONNECTION_SIZE" = 5000000)

file <- files[1]
name <- gsub("POST", "Post",gsub("PRE", "Pre", str_extract(basename(file), "Pt\\d+_(?i:pre|post)\\d+")))
name <- Correspondance$sampleID[Correspondance$oldSampleID == name]
data <- vroom(file, skip = 7, delim = ",",  col_names = TRUE)
data <- as.data.frame(data)
rownames(data) = paste(name, data$Cell_Index, sep = "_")
data$Cell_Index <- NULL
data_t <- t(data)
data_sparse <- as(data_t, "sparseMatrix")

fused <- CreateSeuratObject(counts = data_sparse, project = name , min.cells = 1, min.features = 1)
fused$orig.ident = name

for (file in files[-1]) {
    name=NULL
    name <- gsub("POST", "Post",gsub("PRE", "Pre", str_extract(basename(file), "Pt\\d+_(?i:pre|post)\\d+")))
    if (is.na(name)) {
      name <- gsub("POST", "Post",gsub("PRE", "Pre", str_extract(basename(file), "Pt\\d+_(?i:pre|post)")))
    }
    name <- Correspondance$sampleID[Correspondance$oldSampleID == name]
    print(name)
    data <- vroom(file, skip = 7, delim = ",",  col_names = TRUE)
    data <- as.data.frame(data)
    rownames(data) = paste(name, data$Cell_Index, sep = "_")
    data$Cell_Index <- NULL
    data_t <- t(data)
    data_sparse <- as(data_t, "sparseMatrix")

    dataSeurat <- CreateSeuratObject(counts = data_sparse, project = name , min.cells = 1, min.features = 1)
    dataSeurat$orig.ident = name
    fused <- merge(fused, y = dataSeurat , add.cell.ids = NULL)

  }

fused <- JoinLayers(fused)
count = fused@assays$RNA$counts
CellsAnnot = data.frame("barcodes" = colnames(fused),
                        "sampleID" = fused$orig.ident,
                        "oldSampleID" = Correspondance$oldSampleID[match(fused$orig.ident, Correspondance$sampleID)])


Clinic$organism_ch1 = 'homo_sapiens'
names(Clinic)[names(Clinic) == "organism_ch1"] <- "organism"
Clinic$disease = 'PDAC'
Clinic$patientSampling = 'biopsy'
Clinic$specimenConservation = 'fresh'
Clinic$samplePathologicalState = 'tumor'
Clinic$hadTreatment <- ifelse(grepl("Post", Clinic$sampleID, ignore.case = TRUE),
                              "TRUE",
                              ifelse(grepl("Pre", Clinic$sampleID, ignore.case = TRUE),
                                     "FALSE",
                                     NA_character_))

 
Clinic <- Clinic %>%
  mutate(treatmentInfo = case_when(
    grepl("Pre", sampleID, ignore.case = TRUE) & !grepl("Post", sampleID, ignore.case = TRUE) ~ "naive",
    !grepl("Pre", sampleID, ignore.case = TRUE) & grepl("Post", sampleID, ignore.case = TRUE) ~ "SBRT and FAKi (VS-6063)",
    TRUE ~ NA_character_
  ))

Clinic$specimenOrgan = 'pancreas'
Clinic$tissueUnitExtraction = 'cell'
Clinic$technology = "SingleCell 10x 3' v3.1"
Clinic$samplePathologicalState = 'tumor'

export_data(
  count = count,
  CellsAnnot = CellsAnnot,
  Clinic = Clinic,
  file_path = file.path(path, project, "03VerifiedDataSet")
)

