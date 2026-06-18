####################################################################################################
#  04_reintegrate_epithelium.R
#
#  Epithelial-cell subset + reference-based re-integration.
#
#  Paper:
#   Dapas et al. (2025) "Cellular and molecular dysregulation of the esophageal epithelium in
#   systemic sclerosis." JCI Insight.
#
#  Purpose:
#   Takes the fully annotated all-cell object, subsets the esophageal epithelial cells (EECs;
#   cell-type code "Ep"), and re-runs the same reference-based SCTransform/CCA integration used
#   for the full dataset (02_integrate_allcells.R) on that subset. Cell-cycle difference (CC.diff)
#   is regressed out during SCTransform so that proliferating vs non-proliferating states do not
#   dominate the integrated embedding. The result is the EEC object analyzed in
#   05_epithelial_compartments.R (and downstream).
#
#  Pipeline position:
#   02_integrate_allcells.R -> 03_annotate_celltypes.R (saves All_integratedObj.RData)
#     -> 04_reintegrate_epithelium.R (THIS SCRIPT, saves eso_epi_integratedObj.RData)
#     -> 05_epithelial_compartments.R -> ...
#
#  Inputs:
#   - <results_dir>/All_integratedObj.RData
#       Annotated all-cell Seurat object from 03_annotate_celltypes.R, carrying the per-cell
#       cell-type assignment in meta.data$types.
#
#  Outputs:
#   - <results_dir>/eso_epi_integratedObj.RData
#       Re-integrated EEC Seurat object (integrated assay + PCA), the entry point for the
#       epithelial analyses.
#
#  Usage:
#   Submit via run_reintegration.sh, or:  Rscript --vanilla 04_reintegrate_epithelium.R
#   This is memory-intensive (the original used a ~540-640 Gb high-memory node).
#
#  Notes / flags:
#   - Cell-type subset: types == "Ep" (the stratified squamous EECs). Glandular epithelium
#     ("GEp") is intentionally excluded.
#   - nCount filter: cells with <= NCOUNT_FILTER (3500) UMIs are removed. This matches the red
#     threshold line drawn in . The filter appeared as exploratory (commented) code
#     in the original epi-analysis script; please confirm it was applied to the object used for
#     the published figures, since it changes the EEC cell set.
####################################################################################################


####  LOAD LIBRARIES  ##############################################################################

suppressPackageStartupMessages({
  library(methods)
  library(Seurat)
  library(sctransform)
  library(ggplot2)
  library(glmGamPoi)
  library(cowplot)
  library(scales)
  library(future)
  library(future.apply)
})


####  SET PARAMETERS  ##############################################################################

set.seed(0123456789)
plan('multisession', workers = 4)

## Memory options (sized for the original ~640 Gb RAM request).
Sys.setenv('R_MAX_VSIZE' = 520000000000)
options(future.globals.maxSize = 520000 * 1024^2)

## ---- Paths (relative to repository root) ----
results_dir <- './results'

## ---- Integration reference ----
## Reference-based integration anchors on three healthy-control proximal/distal pairs (same
## reference set used in 02_integrate_allcells.R).
ref_samples <- c('HC3-P', 'HC3-D', 'HC6-P', 'HC6-D', 'HC5-P', 'HC5-D')   # de-identified HC reference pairs

## ---- Subset / integration parameters ----
EPI_CELLTYPE  <- 'Ep'    # cell-type code for the stratified squamous epithelium
NCOUNT_FILTER <- 3500    # minimum UMIs/cell retained (see flag in header)
N_FEATURES    <- 3000    # integration features


####  FUNCTIONS  ###################################################################################

loadRData <- function(fileName) {
  load(fileName)
  get(ls()[ls() != 'fileName'])
}


####  MAIN  ########################################################################################

## Load the annotated all-cell object and subset the epithelial cells.
sc_obj <- loadRData(paste0(results_dir, '/All_integratedObj.RData'))
sc_obj <- subset(sc_obj, subset = types == EPI_CELLTYPE)

## Remove low-UMI cells ( threshold; see header flag).
sc_obj <- subset(sc_obj, subset = nCount_RNA > NCOUNT_FILTER)

## Cell-cycle scoring on the RNA assay so CC.diff is available to regress during SCTransform.
DefaultAssay(sc_obj) <- 'RNA'
sc_obj <- CellCycleScoring(sc_obj, s.features = cc.genes$s.genes,
                           g2m.features = cc.genes$g2m.genes, set.ident = FALSE)
sc_obj$CC.diff <- sc_obj$S.Score - sc_obj$G2M.Score

## Split by sample and SCTransform (v2), regressing the cell-cycle difference.
scrna.list <- SplitObject(sc_obj, split.by = 'orig.ident')
ref_list   <- which(names(scrna.list) %in% ref_samples)

scrna.list <- future_lapply(X = scrna.list, FUN = SCTransform, vst.flavor = 'v2',
                            vars.to.regress = c('CC.diff'))
plan('default')

## Reference-based SCT integration.
int_features <- SelectIntegrationFeatures(scrna.list, nfeatures = N_FEATURES)
scrna.list   <- PrepSCTIntegration(scrna.list, anchor.features = int_features)
scrna.list   <- future_lapply(X = scrna.list, FUN = RunPCA, features = int_features)

anchors <- FindIntegrationAnchors(object.list = scrna.list, normalization.method = 'SCT',
                                  anchor.features = int_features, reference = ref_list)

integrated_set <- IntegrateData(anchorset = anchors, normalization.method = 'SCT')
integrated_set <- RunPCA(object = integrated_set, assay = 'integrated')
DefaultAssay(integrated_set) <- 'integrated'

save(integrated_set, file = paste0(results_dir, '/eso_epi_integratedObj.RData'))
