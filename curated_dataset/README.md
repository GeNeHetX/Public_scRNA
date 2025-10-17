# scRNA data Curation
Cleaned, preprocessed, and curated public scRNA-seq datasets prepared for integrative analyses.

💾 List of dataset : 
[Database summary](https://insermfrance-my.sharepoint.com/:x:/g/personal/camille_pignolet_inserm_fr/Ebhhv4x2GspAifGHqlAOo0oBZtLiLijglbYzuOBcGtm_pw?rtime=3cTYMmUN3kg)

# Instruction
## 📘 Dataset Structure

At the end for one dataset you have:

- → **Raw scRNA matrix** (`matrix.mtx.gz`)  
- → **Gene annotation** (`features.tsv.gz`)  
- → **Barcode annotation** (`barcodes.tsv.gz`)  
  - Must contain at least:  
    `barcodes`, `originalBarcodes`, `sampleID`, `publishedClassL1`

---

### 🤓 Naming Convention

**Sample name example (for Carpenter):**  
`Carp_S01` (you can change the second part but it must start with **S**)  ⇒ **01**, **02**, **03**, **04**, **05**, …

**Patient name example (for Carpenter):**  
`Carp_P01`

---

### 📄 Annotation Files

- → **Sample annotation** (`samplesAnnot.tsv`)  
  Must contain at least:  
  `sampleID`, `oldSampleID`, `patientID`, `patientSampling`, `specimenOrgan`,  
  `specimenConservation`, `samplePathologicalState`, `treatmentInfo`,  
  `hadTreatment`, `tissueUnitExtraction`, `technology`

- → **Patient annotation** (`clinicalAnnot.tsv`)  
  Must contain at least:  
  `patientID`, `disease`, `organism`

---

### 🧩 Script to Check

- → Respect dictionary modalities and column names:  
  [Dictionary link](https://insermfrance-my.sharepoint.com/:x:/g/personal/usama_akhtar_inserm_fr/EStZGXBOnOZEqynK3WyjE6EBSlQqLy6Ie6TctiEsI8pr2A?e=wNK7fL)

- → Verify that `patientID` in **sample annotation** exists in **patient annotation**  
  (also check that `sampleID` in **barcode annotation** exists in **sample annotation**)

- → Ensure all barcodes in **matrix.mtx.gz** are present in **barcode annotation**

