# ==============================================================================
# Reference-based integration of single-cell samples
#
# Part of the analysis for:
#   Dapas et al. (2025) "Cellular and molecular dysregulation of the esophageal
#   epithelium in systemic sclerosis."
#
# Purpose
#   Normalize each QC-passing sample with SCTransform (v2) and integrate all
#   samples with Seurat's reference-based CCA anchoring, using three pairs of
#   proximal/distal healthy-control (HC) samples as the integration reference.
#   (Corresponds to the Methods "Sample integration and cell type annotation"
#   section.)
#
# Pipeline position
#   Downstream of 01_sample_qc.R, upstream of 03_annotate_celltypes.R.
#
# Sample-naming convention (used to assign condition/biopsy location)
#   Samples are read in under their (de-identified) IDs. Condition is the ID prefix
#   ('HC'/'GERD'/'SSc'); biopsy location is the '-P' (Proximal) / '-D' (Distal) suffix.
#
# Inputs
#   <data_dir>/<sample>/analysis/sample_QC/seuratObj_<sample>.RData
#       (per-sample QC-passing Seurat objects from 01_sample_qc.R)
#
# Outputs (written under results_dir)
#   eso_intRefFeatures.RData    selected integration features
#   eso_intRefAnchors.RData     integration anchors
#   eso_integratedRefObj.RData  integrated Seurat object (PCA computed)
# ==============================================================================


####  LOAD LIBRARIES  ##############################################################################

suppressPackageStartupMessages({
  library(methods)
  library(Seurat)
  library(sctransform)
  library(ggplot2)
  library(glmGamPoi)
  library(future)
  library(future.apply)
  library(cowplot)
  library(scales)
})


####  SET PARAMETERS  ##############################################################################

set.seed(0123456789)
plan("multisession", workers = 6)

# Memory options (sized for a high-memory node; adjust to your allocation)
Sys.setenv("R_MAX_VSIZE" = 760000000000)
options(future.globals.maxSize = 760000 * 1024^2)   # ~ for a 640 GB RAM request

# --- Paths (edit to match your environment) -----------------------------------
data_dir    <- "./data/matrices"   # per-sample Cell Ranger / QC subdirectories (named by sample ID)
results_dir <- "./results"         # integration outputs

# Metadata column names
condition <- "condition"
location  <- "location"

# Integration reference: three pairs of proximal/distal HC samples
ref_samples <- c("HC3-P", "HC3-D",
                 "HC6-P", "HC6-D",
                 "HC5-P", "HC5-D")

# Sample excluded for low quality (39/40 retained)
excluded_samples <- c("SSc10-P")


####  FUNCTIONS  ###################################################################################

loadRData <- function(fileName) {
  ## Load a single object from an .RData file and return it (regardless of the
  ## name it was saved under).
  load(fileName)
  get(ls()[ls() != "fileName"])
}


####  MAIN  ########################################################################################

# Sample names = subdirectory names under the data directory, minus exclusions
samples <- list.dirs(data_dir, recursive = FALSE, full.names = FALSE)
samples <- samples[!samples %in% excluded_samples]

# Load each QC-passing object and annotate condition / biopsy location from the sample ID
scrna.list <- list()
for (i in seq_along(samples)) {
  sample <- samples[i]
  print(sample)

  obj <- loadRData(paste0(data_dir, "/", sample,
                          "/analysis/sample_QC/seuratObj_", sample, ".RData"))

  # Condition from the sample-ID prefix (HC / GERD / SSc); location from the -P / -D suffix
  obj[[condition]] <- sub("[0-9].*$", "", sample)
  obj[[location]]  <- ifelse(grepl("-P$", sample), "Proximal", "Distal")

  scrna.list[[i]] <- obj
}
names(scrna.list) <- samples

# Indices of the reference samples within the object list
ref_list <- which(names(scrna.list) %in% ref_samples)

# SCTransform each sample (v2 regularization)
scrna.list <- future_lapply(X = scrna.list, FUN = SCTransform, vst.flavor = "v2")

# Select integration features and prep for SCT integration
int_features <- SelectIntegrationFeatures(scrna.list, nfeatures = 3000)
save(int_features, file = paste0(results_dir, "/eso_intRefFeatures.RData"))

scrna.list <- PrepSCTIntegration(scrna.list, anchor.features = int_features)
scrna.list <- future_lapply(X = scrna.list, FUN = RunPCA, features = int_features)

# Find reference-based integration anchors
anchors <- FindIntegrationAnchors(object.list = scrna.list, normalization.method = "SCT",
                                  anchor.features = int_features, reference = ref_list)
save(anchors, file = paste0(results_dir, "/eso_intRefAnchors.RData"))

# Integrate (single worker for IntegrateData) and run PCA on the integrated assay
plan("multisession", workers = 1)
integrated_set <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

integrated_set <- RunPCA(object = integrated_set, assay = "integrated")
DefaultAssay(integrated_set) <- "integrated"

save(integrated_set, file = paste0(results_dir, "/eso_integratedRefObj.RData"))
