####################################################################################################
#  10_cosmx_spatial.R
#
#  Spatial (CosMx) analysis of the esophageal epithelium: cell-type annotation by label transfer
#  and cell-cell proximity (neighbor-enrichment) analysis.
#
#  Paper:
#   Dapas et al. (2025) "Cellular and molecular dysregulation of the esophageal epithelium in
#   systemic sclerosis." JCI Insight.
#
#  Purpose / figures produced:
#   - Spatial cell-type maps (representative Proximal / Distal FOVs)
#   - Spatial maps colored by superficial subcluster
#   - Spatial metallothionein-signature maps (superficial cells)
#   - Proportion of neighbors of each epithelial compartment, by condition
#   - Change in epithelial-cell neighbor composition, SSc - HC
#   - Fibroblast neighbor composition (-> fibroblast_neighbor_composition.tsv)
#   - Per-sample H&E-matched spatial cell-type maps + cell-type counts; QC scatters
#   - Full neighbor-proportion bar plots
#   ( and  are CellChat panels produced by 09_cellchat_communication.R.)
#
#  Pipeline position:
#   Stands alone on the CosMx data, but uses the scRNA-seq atlas as the label-transfer reference:
#   03_annotate_celltypes.R (cell types), 05_epithelial_compartments.R (compartments),
#   08_superficial_analysis.R (superficial subclusters)  ->  10_cosmx_spatial.R (THIS SCRIPT)
#
#  Inputs:
#   - <nano_root>/<sample>/seuratObject_<fov>.RDS   the 4 CosMx 6k-panel samples
#                                                   (SSc_11, SSc_12, HC_7, HC_8)
#   - <results_dir>/All_integratedObj.RData         scRNA-seq reference for cell-type labels
#   - <results_dir>/eso_epi_integratedObj.RData     reference for epithelial-compartment labels
#   - <results_dir>/eso_epi_superficial.RData       reference for superficial-subcluster labels
#
#  Outputs (to <results_dir> / <suppdata_dir>):
#   - Eso_CosMx_seuratObj.RData, cosMx_labels.RData, cosMx_epiGepiLabels.RData, cosMx_epiSuperLabels.RData
#   - cosMx_spatialProximityAnalysis.RData; adjacency_df_epi.tsv, all_adj_cond.tsv, fibroblast_neighbor_composition.tsv
#   - spatial-map and proximity .pdf / .jpg panels
#
#  Usage:
#   Rscript 10_cosmx_spatial.R   (run from the repository root; memory-intensive)
#
#  Notes:
#   - Removed from the original working file: the exploratory "differentiation score" block, which
#     re-clustered epithelial/immune/stromal/indeterminate subsets to refine annotation (its objects
#     are never used downstream), plus interactive scratch (View(), ElbowPlot/pheatmap dim-selection,
#     and marker feature-plot QC grids).
#   - FLAG: the original loaded pre-computed SingleR label objects (cosMx_labels / cosMx_epiGepiLabels)
#     and only showed the superficial SingleR call. The cell-type and compartment SingleR calls are
#     reconstructed here from the scRNA-seq references so the script runs end-to-end; confirm the
#     reference objects/labels match what was used originally.
#   - The HC_7 Proximal/Distal sections were mislabeled and are swapped (a sample-sheet correction
#     original).
####################################################################################################


####  LOAD LIBRARIES  ##############################################################################

suppressPackageStartupMessages({
  library(methods)
  library(Seurat)
  library(Banksy)
  library(harmony)
  library(SingleR)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(RANN)            # nn2 nearest neighbors
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(cowplot)
  library(gridExtra)
  library(clustree)
  library(RColorBrewer)
  library(scales)
  library(stringr)
})


####  SET PARAMETERS  ##############################################################################

set.seed(0123456789)
options(future.globals.maxSize = 2000 * 1024^2)

## ---- Paths (relative to repository root) ----
nano_root    <- './data/cosmx'          # raw CosMx per-sample Seurat .RDS objects
results_dir  <- './results'
suppdata_dir <- './results/supporting_data'

## ---- CosMx samples ----
## Raw per-sample run directories (read as-is) and their de-identified sample names.
samples <- c('Cos1_6k_5um_101EP_101ED_1_813_27_09_2024_22_07_20_9',
             'Cos1_6k_5um_176EP_176ED_2_813_28_09_2024_22_22_47_208',
             'Cos4_6k_5um_11UP_11LOW_8_816_28_09_2024_22_27_47_905',
             'Cos3_6k_5um_28UP_28LOW_3_813_28_09_2024_22_24_09_283')
sample_names <- c('SSc_12', 'SSc_11', 'HC_7', 'HC_8')
prox_locs    <- list(15:24, 11:23, 1:6, 8:12)   # FOVs assigned to the Proximal biopsy per sample

## FOV image names (derived from the sample IDs).
fovs <- sapply(samples, function(x) {
  s <- gsub('_', '.', word(x, 1, 7, sep = '_'))
  paste(substr(s, 1, nchar(s) - 2), substr(s, nchar(s) - 1, nchar(s)), sep = '.')
}, USE.NAMES = FALSE)

## ---- QC thresholds ----
MIN_COUNT     <- 20      # minimum counts per cell
NEG_PROBE_MAX <- 0.05    # maximum fraction of counts from negative probes

## ---- Clustering / proximity parameters ----
BANKSY_LAMBDA <- 0.2
BANKSY_KGEOM  <- 15
CLUSTER_RES   <- 0.7
K_NEIGHBORS   <- 6       # k for nearest-neighbor search (self excluded -> 5 neighbors)
MAX_DISTANCE  <- 0.05    # 50 um (coordinates in mm)

## ---- Colors ----
cond_cols    <- c('#8156B3', '#e26b53', '#A8A39d')                          # SSc, GERD, HC
layer_cols_5 <- c('#3288BD', '#66C2A5', '#E6F598', '#e9b85d', '#c72f4c')    # 5 epithelial compartments
cell_type_cols <- c('#CD96CD', '#8B5A00', '#7AC5CD', '#4682B4', '#EEA2AD', '#8B2323')  # En, F, L, Mc, My, P
## Combined palette for the compartment-resolved cell types (En, 5 compartments, then others).
all_cols     <- c('#CD96CD', layer_cols_5, '#8B5A00', '#7AC5CD', '#4682B4', '#EEA2AD', '#8B2323', '#FF7F50', '#76EE00')

anno_genes <- c('PDPN', 'IGFBP3', 'COL17A1', 'KRT15', 'DST', 'KRT4', 'KRT13', 'SERPINB3', 'GBP6', 'IVL',
                'KRT17', 'KRT78', 'FLG', 'CNFN', 'CRCT1', 'KRT6A', 'KRT14', 'MUC5B', 'PTPRC', 'CD3D',
                'LYZ', 'HLA-DRA', 'TPSAB1', 'PECAM1', 'PDGFRA', 'PDGFRB', 'RGS5', 'ACTA2')

cell_types  <- c('En', 'Ep', 'F', 'GEp', 'L', 'Mc', 'My', 'P', 'SMCs', 'SMG')
metallo_genes <- c('MT1A', 'MT1E', 'MT1F', 'MT1G', 'MT1H', 'MT1M', 'MT1X', 'MT2A')


####  FUNCTIONS  ###################################################################################

loadRData <- function(fileName) {
  load(fileName)
  get(ls()[ls() != 'fileName'])
}

## kNN neighbor finder within one FOV: for each query cell (epithelial, non-epithelial, or all),
## return the cell-type of its k nearest neighbors within max_dist.
find_neighbors_all_types <- function(fov_name, k = K_NEIGHBORS, max_dist = MAX_DISTANCE, query_types = NULL, all = FALSE) {
  coords <- GetTissueCoordinates(epi_merged, image = fov_name)
  if (is.null(coords) || nrow(coords) == 0) return(NULL)
  cell_ids  <- coords$cell
  cell_meta <- epi_merged@meta.data[cell_ids, c('is_epi', 'pred.epiLayers', 'pred.cell_types', 'condition', 'location', 'condLoc')]
  cell_meta$pred.epiLayers <- as.character(cell_meta$pred.epiLayers)

  if (all)                       query_idx <- seq_len(nrow(cell_meta))
  else if (!is.null(query_types)) query_idx <- which(cell_meta$is_epi == 0)
  else                            query_idx <- which(cell_meta$is_epi == 1)

  coords_mat <- as.matrix(coords[, c('x', 'y')])
  nn_result  <- nn2(coords_mat, coords_mat[query_idx, ], k = min(k + 1, nrow(coords_mat)))
  nbr_idx    <- nn_result$nn.idx[, -1, drop = FALSE]      # drop self
  nbr_dist   <- nn_result$nn.dists[, -1, drop = FALSE]
  n_nbr      <- ncol(nbr_idx)

  df <- data.frame(
    source_type   = rep(cell_meta$pred.epiLayers[query_idx], each = n_nbr),
    neighbor_type = as.vector(t(cell_meta$pred.epiLayers[as.vector(t(nbr_idx))])),
    distance      = as.vector(t(nbr_dist)),
    condition     = rep(cell_meta$condition[query_idx], each = n_nbr),
    location      = rep(cell_meta$location[query_idx], each = n_nbr),
    condLoc       = rep(cell_meta$condLoc[query_idx], each = n_nbr),
    fov           = fov_name, stringsAsFactors = FALSE)
  df[df$distance <= max_dist, ]
}


####  CREATE OBJECT (read + QC the 4 CosMx samples)  ##############################################

obj.list <- list()
for (i in 1:4) {
  nano.obj <- readRDS(paste0(nano_root, '/', samples[i], '/seuratObject_', fovs[i], '.RDS'))
  Idents(nano.obj) <- sample_names[i]

  ## Drop pre-computed RNA/cluster/nn columns from the vendor object.
  for (pat in c('^RNA', '^spatialclust', '^nn_'))
    nano.obj@meta.data <- nano.obj@meta.data[, -grep(pat, colnames(nano.obj@meta.data)), drop = FALSE]

  ## Cell QC.
  nano.obj$QC_minCount               <- nano.obj$nCount_RNA >= MIN_COUNT
  nano.obj$QC_negProbProp            <- nano.obj$nCount_negprobes < nano.obj$nCount_RNA * NEG_PROBE_MAX
  nano.obj$QC_impossibleFeatureCount <- nano.obj$nFeature_RNA <= nano.obj$nCount_RNA
  nano.obj$QC                        <- ifelse(nano.obj$qcCellsFlagged == 1, 'FAIL', 'PASS')

  ## Doublet detection on the passing cells.
  nano.scDbl <- NormalizeData(subset(nano.obj, QC == 'PASS'), scale.factor = median(colSums(nano.obj@assays$RNA$counts)))
  nano.scDbl <- scDblFinder(GetAssayData(nano.scDbl, layer = 'counts'), dbr.sd = 0)
  nano.obj@meta.data[nano.obj$QC == 'PASS', 'scDblFinder.class'] <- nano.scDbl$scDblFinder.class
  nano.obj@meta.data[nano.obj$QC == 'PASS', 'QC'] <- ifelse(nano.scDbl$scDblFinder.class == 'doublet', 'FAIL', 'PASS')
  rm(nano.scDbl)

  ## QC scatter panels ( QC).
  nano.obj$pass_n <- ifelse(nano.obj$QC == 'PASS', 'Pass', 'Fail')
  qc_plots <- lapply(list(c('nCount_RNA', 'nFeature_RNA'), c('nCount_RNA', 'Area.um2'),
                          c('Area.um2', 'complexity'), c('nCount_RNA', 'complexity')), function(fp) {
    FeatureScatter(nano.obj, feature1 = fp[1], feature2 = fp[2], group.by = 'pass_n') +
      ggtitle(paste('r =', round(cor(nano.obj[[fp[1]]][, 1], nano.obj[[fp[2]]][, 1]), 2))) + theme(legend.position = 'none')
  })
  pdf(file = paste0(results_dir, '/QC_scatters.', sample_names[i], '.pdf'), width = 8, height = 6, onefile = FALSE)
  print(plot_grid(plotlist = qc_plots, ncol = 2, labels = c('A', 'B', 'C', 'D')))
  dev.off()

  ## Keep passing cells; assign biopsy location by FOV.
  nano.obj <- subset(nano.obj, QC == 'PASS')
  nano.obj$location <- ifelse(nano.obj$fov %in% prox_locs[[i]], 'Proximal', 'Distal')
  nano.obj$sample   <- sample_names[i]
  obj.list[[sample_names[i]]] <- nano.obj
}


####  NORMALIZATION + SPATIAL EMBEDDING (Banksy / Harmony)  #######################################

obj.list <- lapply(obj.list, function(x) { DefaultAssay(x) <- 'RNA'; NormalizeData(x) })
hvgs <- lapply(obj.list, function(x) VariableFeatures(FindVariableFeatures(x, selection.method = 'vst', nfeatures = 250)))

epi_merged <- Reduce(merge, obj.list)
DefaultAssay(epi_merged) <- 'RNA'
VariableFeatures(epi_merged) <- Reduce(union, hvgs)
epi_merged$sampleFOVs <- paste(epi_merged$sample, epi_merged$fov, sep = '_')
epi_merged$sampleLoc  <- paste(epi_merged$sample, epi_merged$location, sep = '_')

epi_merged <- RunBanksy(epi_merged, lambda = BANKSY_LAMBDA, dimx = 'x_slide_mm', dimy = 'y_slide_mm',
                        assay = 'RNA', slot = 'data', features = 'variable', split.scale = TRUE,
                        use_agf = FALSE, k_geom = BANKSY_KGEOM, group = 'sampleLoc', verbose = TRUE)
epi_merged <- RunPCA(epi_merged, assay = 'BANKSY', reduction.name = 'pca.banksy', features = rownames(epi_merged), npcs = 40)
epi_merged <- RunHarmony(epi_merged, group.by.vars = 'sampleFOVs', reduction = 'pca.banksy', reduction.save = 'harmony.banksy', verbose = FALSE)
epi_merged <- RunUMAP(epi_merged, reduction = 'harmony.banksy', dims = 1:30)
epi_merged <- FindNeighbors(epi_merged, k = 21, reduction = 'harmony.banksy', dims = 1:30)
epi_merged <- FindClusters(epi_merged, cluster.name = 'banksy_cluster', resolution = CLUSTER_RES, algorithm = 1, n.iter = 50, n.start = 50)

for (fov in fovs) DefaultBoundary(epi_merged[[fov]]) <- 'centroids'
save(epi_merged, file = paste0(results_dir, '/Eso_CosMx_seuratObj.RData'))


####  SingleR LABEL TRANSFER  #####################################################################
## Transfer cell-type, epithelial-compartment, and superficial-subcluster labels from the
## scRNA-seq atlas onto the CosMx cells.

cosmx_sce <- as.SingleCellExperiment(epi_merged, assay = 'RNA')

## (1) Main cell types (reference: the annotated all-cell object).
ref_all <- loadRData(paste0(results_dir, '/All_integratedObj.RData'))
ref_all_sce <- as.SingleCellExperiment(ref_all, assay = 'RNA')
pred.labels.cosmx.eso <- SingleR(test = cosmx_sce, ref = ref_all_sce, labels = ref_all$types)
save(pred.labels.cosmx.eso, file = paste0(results_dir, '/cosMx_labels.RData'))
epi_merged$pred.cell_types <- factor(pred.labels.cosmx.eso$pruned.labels, levels = cell_types)
## Fold the (rare) glandular-epithelial calls into Ep for the compartment step.
epi_merged$pred.cell_types <- ifelse(epi_merged$pred.cell_types == 'GEp', 'Ep', as.character(epi_merged$pred.cell_types))

## (2) Epithelial compartments (reference: the EEC object, layers_rep).
ref_epi <- loadRData(paste0(results_dir, '/eso_epi_integratedObj.RData'))
ref_epi_sce <- as.SingleCellExperiment(ref_epi, assay = 'RNA')
cosmx_epi_sce <- as.SingleCellExperiment(subset(epi_merged, pred.cell_types == 'Ep'), assay = 'RNA')
pred.labels.cosmx.epiGepi <- SingleR(test = cosmx_epi_sce, ref = ref_epi_sce, labels = ref_epi$layers_rep)
pred.labels.cosmx.epiGepi$pruned.labels <- mapvalues(pred.labels.cosmx.epiGepi$pruned.labels,
                                                     from = c('superficial', 'suprabasal', 'replicating suprabasal', 'replicating basal', 'basal'),
                                                     to   = c('Ep: superficial', 'Ep: suprabasal', 'Ep: prolif. suprabasal', 'Ep: prolif. basal', 'Ep: basal'))
save(pred.labels.cosmx.epiGepi, file = paste0(results_dir, '/cosMx_epiGepiLabels.RData'))

epi_merged$pred.epiLayers <- as.character(epi_merged$pred.cell_types)
epi_merged@meta.data[rownames(pred.labels.cosmx.epiGepi), 'pred.epiLayers'] <- pred.labels.cosmx.epiGepi$pruned.labels
epi_merged$pred.epiLayers <- factor(epi_merged$pred.epiLayers,
                                    levels = c('Ep: basal', 'Ep: prolif. basal', 'Ep: prolif. suprabasal', 'Ep: suprabasal',
                                               'Ep: superficial', 'SMG', 'En', 'F', 'L', 'Mc', 'My', 'P', 'SMCs'))

## (3) Superficial subclusters (reference: super_set seurat_clusters).
super_set <- loadRData(paste0(results_dir, '/eso_epi_superficial.RData'))
super_sce <- as.SingleCellExperiment(super_set, assay = 'RNA')
cosmx_super_sce <- as.SingleCellExperiment(subset(epi_merged, pred.epiLayers == 'Ep: superficial'), assay = 'RNA')
pred.labels.cosmx.super <- SingleR(test = cosmx_super_sce, ref = super_sce, labels = super_set$seurat_clusters)
save(pred.labels.cosmx.super, file = paste0(results_dir, '/cosMx_epiSuperLabels.RData'))
epi_merged$pred.superCluster <- NA
epi_merged@meta.data[rownames(pred.labels.cosmx.super), 'pred.superCluster'] <- as.character(pred.labels.cosmx.super$pruned.labels)

## Sample-sheet correction: the HC_7 Proximal/Distal sections were mislabeled and are swapped.
hc7 <- epi_merged$sample == 'HC_7'
epi_merged@meta.data[hc7, 'location'] <- mapvalues(epi_merged@meta.data[hc7, 'location'], c('Proximal', 'Distal'), c('Distal', 'Proximal'))
epi_merged$sampleLoc <- paste(epi_merged$sample, epi_merged$location, sep = '_')
epi_merged$condition <- ifelse(grepl('^SSc', epi_merged$sample), 'SSc', 'HC')
epi_merged$condLoc   <- paste(epi_merged$condition, epi_merged$location, sep = '_')

## Metallothionein signature; restricted to superficial cells for display.
epi_merged <- AddModuleScore(epi_merged, features = list(metallo_genes[metallo_genes %in% rownames(epi_merged)]),
                             name = 'metallothioneins', assay = 'RNA')
epi_merged$metallothioneins_super <- epi_merged$metallothioneins1
epi_merged@meta.data[epi_merged$pred.epiLayers != 'Ep: superficial', 'metallothioneins_super'] <- min(epi_merged$metallothioneins1)

save(epi_merged, file = paste0(results_dir, '/Eso_CosMx_seuratObj.predLabels.RData'))

## Ordered slices for the per-sample spatial panels (Proximal, Distal per sample).
slices <- paste(rep(sample_names, each = 2), c('Proximal', 'Distal'), sep = '_')


####     -  spatial cell-type maps + counts  ###################################

## Per-sample spatial maps colored by cell type (representative panels).
plot_list <- list()
for (i in 1:4) {
  for (h in 1:2) {
    sl <- slices[(i - 1) * 2 + h]
    plot_list[[(i - 1) * 2 + h]] <- ggplotGrob(
      ImageDimPlot(subset(epi_merged, sampleLoc == sl), fov = fovs[i], size = 1.2, group.by = 'pred.cell_types',
                   cols = cell_type_cols, dark.background = FALSE, border.color = 'black', border.size = 0.025, crop = TRUE) +
        theme(legend.position = 'none') + ggtitle(sl))
  }
}
jpeg(file = paste0(results_dir, '/epiSpatial_allSamples_cellType.jpg'), width = 10, height = 14, units = 'in', res = 450)
grid.arrange(grobs = plot_list, ncol = 2)
dev.off()

## Cell-type counts per sample and overall proportions.
write.table(as.data.frame(table(epi_merged$sampleLoc, epi_merged$pred.cell_types)),
            paste0(suppdata_dir, '/Suppsuperficial_marker_expression.tsv'), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')

## Representative Proximal / Distal compartment maps.
pdf(file = paste0(results_dir, '/epiSpatial_compartments.proxEx.pdf'), width = 7, height = 7)
print(ImageDimPlot(subset(epi_merged, fov %in% 20), fov = fovs[1], size = 1.2, group.by = 'pred.epiLayers',
                   cols = all_cols, dark.background = FALSE, border.color = 'black', border.size = 0.05, crop = TRUE) +
        labs(fill = 'Cell type') + scale_y_reverse())
dev.off()
pdf(file = paste0(results_dir, '/epiSpatial_compartments.distEx.pdf'), width = 7, height = 7)
print(ImageDimPlot(subset(epi_merged, fov %in% c(6, 7)), fov = fovs[2], size = 1.2, group.by = 'pred.epiLayers',
                   cols = all_cols, dark.background = FALSE, border.color = 'black', border.size = 0.05, crop = TRUE) +
        labs(fill = 'Cell type') + scale_x_reverse())
dev.off()


####  Spatial superficial-subcluster maps  #############################################
super_cols <- c('#1F9E89', '#E76254', '#67AD5B', '#E3B23C', '#9C8CC4')
plot_list <- list()
for (i in 1:4) {
  for (h in 1:2) {
    sl <- slices[(i - 1) * 2 + h]
    plot_list[[(i - 1) * 2 + h]] <- ggplotGrob(
      ImageDimPlot(subset(epi_merged, sampleLoc == sl), fov = fovs[i], size = 1.2, group.by = 'pred.superCluster',
                   cols = super_cols, dark.background = FALSE, border.color = 'black', border.size = 0.025, crop = TRUE) +
        theme(legend.position = 'none') + ggtitle(sl))
  }
}
jpeg(file = paste0(results_dir, '/epiSpatial_allSamples_superCluster.jpg'), width = 10, height = 14, units = 'in', res = 450)
grid.arrange(grobs = plot_list, ncol = 2)
dev.off()


####  Spatial metallothionein maps  ####################################################
metal_cols <- c('#E5EDD7', 'grey85', 'grey75', brewer.pal(n = 9, 'OrRd')[4:9], 'black')
plot_list <- list()
for (i in 1:4) {
  for (h in 1:2) {
    sl <- slices[(i - 1) * 2 + h]
    plot_list[[(i - 1) * 2 + h]] <- ggplotGrob(
      ImageFeaturePlot(subset(epi_merged, sampleLoc == sl), fov = fovs[i], size = 1.2, features = 'metallothioneins_super',
                       cols = metal_cols, dark.background = FALSE, border.color = 'black', border.size = 0.025, crop = TRUE) +
        theme(legend.position = 'none') + ggtitle(sl))
  }
}
jpeg(file = paste0(results_dir, '/epiSpatial_allSamples_metallothionein.jpg'), width = 10, height = 14, units = 'in', res = 450)
grid.arrange(grobs = plot_list, ncol = 2)
dev.off()


####  Cell proximity analysis  ###################################
## For each cell, tabulate the cell-type of its k nearest neighbors within MAX_DISTANCE, then
## summarize neighbor composition by source compartment, condition, and location.

epi_merged$is_epi <- ifelse(epi_merged$pred.cell_types == 'Ep', 1, 0)

adj_epi <- list(); adj_other <- list(); adj_all <- list()
for (fov_name in fovs) {
  adj_all[[fov_name]]   <- find_neighbors_all_types(fov_name, all = TRUE)
  adj_epi[[fov_name]]   <- find_neighbors_all_types(fov_name)
  adj_other[[fov_name]] <- find_neighbors_all_types(fov_name, query_types = 'non-epithelial')
}
adjacency_df_epi   <- do.call(rbind, adj_epi)
adjacency_df_other <- do.call(rbind, adj_other)
adjacency_df_all   <- do.call(rbind, adj_all)

lev <- levels(epi_merged$pred.epiLayers)
adjacency_df_epi$source_type   <- factor(adjacency_df_epi$source_type, levels = lev)
adjacency_df_epi$neighbor_type <- factor(adjacency_df_epi$neighbor_type, levels = lev)
adjacency_df_other$neighbor_type <- factor(adjacency_df_other$neighbor_type, levels = lev)
write.table(adjacency_df_epi, paste0(results_dir, '/adjacency_df_epi.tsv'), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')

## Count adjacency pairs and convert to proportions of neighbors within each grouping.
count_adj <- function(df, ...) {
  extra <- c(...)
  grp   <- c('source_type', 'neighbor_type', extra)
  out <- aggregate(as.formula(paste('distance ~', paste(grp, collapse = '+'))), data = df, FUN = length)
  names(out)[names(out) == 'distance'] <- 'n'
  ## Proportion of each neighbor_type within source_type (x condition/location, if present).
  grp_cols <- c('source_type', extra)
  groups   <- if (length(grp_cols) == 1) out[[grp_cols]] else do.call(interaction, c(out[grp_cols], list(drop = TRUE)))
  out$prop <- ave(out$n, groups, FUN = function(x) x / sum(x))
  out
}
all_adj_epi           <- count_adj(adjacency_df_epi)
all_adj_epi_cond      <- count_adj(adjacency_df_epi, 'condition')
all_adj_epi_condLoc   <- count_adj(adjacency_df_epi, 'condition', 'location')
all_adj_other         <- count_adj(adjacency_df_other)
all_adj_other_cond    <- count_adj(adjacency_df_other, 'condition')
all_adj_other_condLoc <- count_adj(adjacency_df_other, 'condition', 'location')
all_adj_cond          <- count_adj(adjacency_df_all, 'condition')
write.table(all_adj_cond, paste0(results_dir, '/all_adj_cond.tsv'), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')

##  source: fibroblast neighbor composition (non-self).
write.table(all_adj_other_cond[all_adj_other_cond$source_type != all_adj_other_cond$neighbor_type & all_adj_other_cond$source_type == 'F', ],
            paste0(suppdata_dir, '/fibroblast_neighbor_composition.tsv'), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')

##  8D: proportion of neighbors of each epithelial compartment.
neighbor_bar <- function(dat, facet = NULL, xlab = 'Epithelial Subtype') {
  g <- ggplot(dat, aes(x = source_type, y = prop, fill = neighbor_type)) + geom_col() +
    scale_y_continuous(labels = scales::percent) + scale_fill_manual(values = all_cols) +
    labs(x = xlab, y = 'Proportion of Neighbors', fill = 'Neighbor Type') +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  if (!is.null(facet)) g <- g + facet_grid(facet)
  g
}
pdf(file = paste0(results_dir, '/epi_neighbors.cond.pdf'), width = 10, height = 4, onefile = FALSE)
print(neighbor_bar(all_adj_epi_cond, facet = as.formula('~condition')))
dev.off()
pdf(file = paste0(results_dir, '/other_neighbors.pdf'), width = 6, height = 4, onefile = FALSE)
print(neighbor_bar(all_adj_other, xlab = 'Cell Type'))
dev.off()
pdf(file = paste0(results_dir, '/epi_neighbors.condLoc.pdf'), width = 10, height = 8, onefile = FALSE)
print(neighbor_bar(all_adj_epi_condLoc, facet = as.formula('condition ~ location')))
dev.off()

## : change in epithelial-cell neighbor composition, SSc - HC.
epiAll_cond <- aggregate(n ~ neighbor_type + condition, data = all_adj_epi_cond, FUN = sum)
epiAll_cond$prop <- ave(epiAll_cond$n, epiAll_cond$condition, FUN = function(x) x / sum(x))
changes <- data.frame(neighbor_type = epiAll_cond$neighbor_type[epiAll_cond$condition == 'SSc'],
                      change_ssc_hc = epiAll_cond$prop[epiAll_cond$condition == 'SSc'] - epiAll_cond$prop[epiAll_cond$condition == 'HC'])
changes <- changes[order(changes$change_ssc_hc, decreasing = TRUE), ]
changes$neighbor_type <- factor(changes$neighbor_type, levels = changes$neighbor_type)
pdf(file = paste0(results_dir, '/epi_neighbors_change_SSc_HC.pdf'), width = 5, height = 4)
print(ggplot(changes, aes(x = neighbor_type, y = change_ssc_hc, fill = neighbor_type)) + geom_col() +
        scale_fill_manual(values = all_cols) + scale_y_continuous(labels = scales::percent) +
        labs(x = '', y = expression(Delta * ' % EEC neighbors (SSc - HC)')) + theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none'))
dev.off()

save(all_adj_epi, all_adj_epi_condLoc, all_adj_other, all_adj_other_condLoc,
     file = paste0(results_dir, '/cosMx_spatialProximityAnalysis.RData'))
