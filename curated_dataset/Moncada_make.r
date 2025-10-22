library(stringr)
library(purrr)
library(dplyr)
library(Matrix)

source("/mnt/d/06_SingleCellDataset/FonctionSingleCell.r")


project <- "Moncada"
path <- "/mnt/d/06_SingleCellDataset"
rawpath <- file.path(path, project, "01RawData")

CreateDataset(project = project, path = path)

download.file("https://www.dropbox.com/s/t8fxmri99ji93zj/Meta-data_Moncada2020_Pancreas.csv?dl=1",
              file.path(rawpath,"Meta-data_Moncada2020_Pancreas.csv"))


download.file("https://www.dropbox.com/sh/8n2w3t8queawf9p/AADoUWm2e9VhEuEc_pvkyNdha?dl=1",
              file.path(rawpath,"Data_Moncada2020_Pancreas.zip"))
unzip(file.path(rawpath,"Data_Moncada2020_Pancreas.zip"), exdir = rawpath)


Clinic = as.data.frame(data.table::fread(file = file.path(rawpath,"Meta-data_Moncada2020_Pancreas.csv")))
names(Clinic)[names(Clinic) == "sample"] <- "oldSampleID"
Clinic$patientID <- sprintf("Monc_P%02d", match(Clinic$patient, unique(Clinic$patient)))
Clinic$sampleID <- sprintf("Monc_S%02d", match(Clinic$patientID, unique(Clinic$patientID)))
names(Clinic)[names(Clinic) == "site"] <- "specimenOrgan"
Clinic$hadTreatment = "FAUX"
names(Clinic)[names(Clinic) == "treated_naive"] <- "treatmentInfo"
Clinic$tissueUnitExtraction = 'cell'
Clinic$organism = 'homo_sapiens'
Clinic$patientSampling = 'surgery'
Clinic$specimenConservation = 'fresh'
Clinic$samplePathologicalState = 'tumor'
Clinic$disease = 'PDAC'


CellsAnnot = as.data.frame(data.table::fread(file = file.path(rawpath,"Cells.csv")))
sample_map  <- setNames(Clinic$sampleID, Clinic$oldSampleID)
CellsAnnot$sampleID  <- sample_map[CellsAnnot$sample]
prefix <- sub("(^[^_]+_[^_]).*", "\\1", CellsAnnot$cell_name)
new_prefix <- sample_map[prefix]

CellsAnnot$cell_name <- paste0(new_prefix, sub(".*(^[^_]+_[^_])", "\\2", CellsAnnot$cell_name))
names(CellsAnnot)[names(CellsAnnot) == "cell_name"] <- "barcodes"
names(CellsAnnot)[names(CellsAnnot) == "sample"] <- "oldSampleID"
names(CellsAnnot)[names(CellsAnnot) == "cell_type"] <- "publishedClassL1"
names(CellsAnnot)[names(CellsAnnot) == "cell_subtype"] <- "publishedClassL2"


count <- readMM(file.path(rawpath,"Exp_data_UMIcounts.mtx"))
genes <- readLines(file.path(rawpath,"Genes.txt"))
length(genes) == nrow(count)  
rownames(count) <- genes

export_data(
  count = count,
  CellsAnnot = CellsAnnot,
  Clinic = Clinic,
  file_path = file.path(path, project, "03VerifiedDataSet")
)
