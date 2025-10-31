library(stringr)
library(purrr)
library(dplyr)
library(GEOquery)
library(Matrix)
library(R.utils)
library(Seurat)

source("/mnt/d/06_SingleCellDataset/FonctionSingleCell.r")


project <- "Hajkarim"
path <- "/mnt/d/06_SingleCellDataset"
rawpath <- file.path(path, project, "01RawData")

CreateDataset(project = project, path = path)

Clinic=getGEO('GSE253429',GSEMatrix=TRUE)[["GSE253429_series_matrix.txt.gz"]]@phenoData@data
Clinic$title <- sub(".*\\[(.*)\\].*", "\\1", Clinic$title)
names(Clinic)[names(Clinic) == "title"] <- "oldSampleID"
Clinic$patientID <- sprintf("Hajk_P%02d", match(Clinic$oldSampleID, unique(Clinic$oldSampleID)))
Clinic$sampleID <- sprintf("Hajk_S%02d", match(Clinic$patientID, unique(Clinic$patientID)))

download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE253429&format=file",
              file.path(rawpath,"GSE253429_RAW.tar"),
              method = "wget") 
untar(tarfile = file.path(rawpath,"GSE253429_RAW.tar"),
      exdir = rawpath, list = F)

mapping <- data.frame(
  oldSampleID = Clinic$oldSampleID, 
  sampleID = Clinic$sampleID   
)

files <- list.files(rawpath, pattern = "\\.h5$", full.names = TRUE)
file.rename(files, file.path(rawpath, sub("^GSM[0-9]+_", "", basename(files))))
file.rename(
  file.path(rawpath, paste0(mapping$oldSampleID, "_filtered_feature_bc_matrix.h5")),
  file.path(rawpath, paste0(mapping$sampleID, "_filtered_feature_bc_matrix.h5"))
)

data <- Read10X_h5(file.path(rawpath, paste0(Clinic$sampleID[1], "_filtered_feature_bc_matrix.h5")))
colnames(data)= paste0(Clinic$sampleID[1],"_",colnames(data))
fused = CreateSeuratObject(counts = data, project = Clinic$sampleID[1])

for (i in Clinic$sampleID[-1]) {
    message("Processing sample: ", i)
    data <- Read10X_h5(file.path(rawpath, paste0(i, "_filtered_feature_bc_matrix.h5")))
    seurat_obj <- CreateSeuratObject(counts = data, project = i)
    colnames(seurat_obj)= paste0(i,"_",colnames(seurat_obj))

    fused <- merge(fused, y = seurat_obj, add.cell.ids = NULL)
}
gc()
fused <- JoinLayers(fused)
count = fused@assays$RNA$counts

CellsAnnot <- transform(
  data.frame(
    barcodes = colnames(fused),
    sampleID = sapply(strsplit(colnames(fused), "_"), function(x) paste(x[1:2], collapse = "_"))),
    oldSampleID = Clinic$oldSampleID[match(sampleID, Clinic$sampleID)]
)

Clinic$organism = 'homo_sapiens'
names(Clinic)[names(Clinic) == "recurrence site:ch1"] <- "specimenOrgan"
Clinic$tissueUnitExtraction = 'nucleus'
Clinic$disease = 'PDAC'
Clinic$patientSampling = 'surgery'
Clinic$specimenConservation = 'frozen'
Clinic$samplePathologicalState = 'tumor'
Clinic <- Clinic %>%
  mutate(treatmentInfo = case_when(
    `neoadjuvant therapy:ch1` == "no" ~ "naive",
    TRUE                     ~ "neoadjuvant"
  ))

Clinic <- Clinic %>%
  mutate(hadTreatment = case_when(
    treatmentInfo == "naive" ~ "FALSE",
    TRUE                     ~ "TRUE"
  ))
Clinic$technology = "10x Genomics 5'"

export_data(
  count = count,
  CellsAnnot = CellsAnnot,
  Clinic = Clinic,
  file_path = file.path(path, project, "03VerifiedDataSet")
)
