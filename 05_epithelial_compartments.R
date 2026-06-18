####################################################################################################
#  05_epithelial_compartments.R
#
#  Esophageal epithelial cell (EEC) characterization: differentiation compartments and their
#  proportional shifts by disease condition and esophageal location.
#
#  Paper:
#   Dapas et al. (2025) "Cellular and molecular dysregulation of the esophageal epithelium in
#   systemic sclerosis." JCI Insight.
#
#  Purpose / figures produced:
#   - EEC UMAP colored by differentiation compartment (5 compartments)
#   - Stacked violin of compartment markers across the 5 compartments
#   - Compartment proportions by condition (HC/GERD/SSc) and location
#                   (Proximal/Distal), with MASC + propeller significance testing
#   - EEC QC feature plots, the nCount (UMI) >3500 filter violin (E),
#                   compartment-marker module scores (G), clustering/cluster UMAP (I),
#                   and cell-cycle phase (L)
#   Also computes the per-compartment single-cell DEG lists that are consumed downstream by
#   07_gsea_enrichment.R (GSEA / TF enrichment).
#
#  Pipeline position:
#   02_integrate_allcells.R -> 03_annotate_celltypes.R (cell typing, saves All_integratedObj.RData)
#     -> 04_reintegrate_epithelium.R (epithelial Ep cells subset out and re-integrated)
#     -> 05_epithelial_compartments.R (THIS SCRIPT)
#     -> 06_pseudobulk_dgea.R / 07_gsea_enrichment.R / 08_superficial_analysis.R
#
#  Inputs:
#   - <results_dir>/eso_epi_integratedObj.RData
#       A Seurat object of the epithelial (Ep) cells only, produced by 04_reintegrate_epithelium.R
#       (subset of the annotated All_integratedObj.RData, low-UMI cells removed, then re-integrated
#       with the same reference-based SCTransform/CCA workflow as 02_integrate_allcells.R, regressing
#       the cell-cycle difference). The object already carries the published 35-PC UMAP embedding.
#
#  Outputs (written to <results_dir>):
#   - eso_epi_integratedObj.RData          re-saved object with cell-cycle scores, compartment
#                                          labels (layers / layers_rep), and differentiation score
#   - epi_cluster_props.tsv                per-sample cluster proportions
#   - epi_layer_props.tsv                  per-sample 3-compartment proportions
#   - epi_layerRep_props.tsv               per-sample 5-compartment proportions
#   - diffScoreStats.tsv                   per-sample differentiation-score summary (for 11)
#   - epi_layerRep_props_stats.tsv         MASC + propeller proportion statistics
#   - epi_degs_ssc_hc_layers.RData         per-compartment single-cell DEGs (SSc vs HC)
#   - epi_degs_gerd_hc_layers.RData        per-compartment single-cell DEGs (GERD vs HC)
#   - epi_degs_ssc_gerd_layers.RData       per-compartment single-cell DEGs (SSc vs GERD)
#   - epi_degs_<comp>_layersRep_<Proximal|Distal>.RData
#                                          5-compartment single-cell DEGs split by location
#                                          (inputs to the  GSEA in 07_gsea_enrichment.R)
#   - assorted .pdf figure panels
#
#  Usage:
#   Run from the repository root (so the relative paths below resolve), e.g.:
#     Rscript 05_epithelial_compartments.R
#   on a high-memory node (the object is large; ~160 Gb was used originally).
####################################################################################################


####  LOAD LIBRARIES  ##############################################################################

suppressPackageStartupMessages({
  library(methods)
  library(Seurat)
  library(SeuratWrappers)
  library(ggplot2)
  library(ggpubr)
  library(glmGamPoi)
  library(cowplot)
  library(plyr)
  library(RColorBrewer)
  library(future)
  library(future.apply)
  library(scCustomize)
  library(R.utils)
  library(patchwork)
  library(scales)
  library(beeswarm)
  library(lme4)      # required by the MASC() mixed-model proportion test (was missing originally)
  library(speckle)   # provides propeller() for the proportion test (was missing originally)
  library(stringr)
})


####  SET PARAMETERS  ##############################################################################

set.seed(0123456789)
plan('multisession', workers = 2)

## Memory options (sized for the original ~160 Gb RAM allocation)
Sys.setenv('R_MAX_VSIZE' = 192000000000)
options(future.globals.maxSize = 192000 * 1024^2)

## ---- Paths (relative to repository root) ----------------------------------------------------
## Point these at wherever the integrated objects live in your checkout.
results_dir <- './results'

## ---- Metadata column names -------------------------------------------------------------------
cond <- 'condition'
loc  <- 'location'

## ---- Cohort / grouping -----------------------------------------------------------------------
conditions <- c('HC', 'GERD', 'SSc')   # display order (HC reference -> GERD -> SSc)
locations  <- c('Proximal', 'Distal')

## ---- Color palettes --------------------------------------------------------------------------
## cond_cols is stored in (SSc, GERD, HC) order; most plots use rev(cond_cols) to get HC->SSc.
cond_cols <- c('#8156B3', '#e26b53', '#A8A39d')        # SSc (purple), GERD (orange), HC (grey)

layer_cols_3 <- c('#4ca5b1', '#e9b85d', '#c72f4c')     # basal, suprabasal, superficial
layer_cols_5 <- brewer.pal(n = 11, 'Spectral')[c(10, 9, 7, 4, 2)]
layer_cols_5[4:5] <- layer_cols_3[2:3]                 # match 3- and 5-compartment palettes

## ---- Differentiation-compartment marker modules ---------------------------------------------
## These three modules define the basal / suprabasal / superficial axes used both for the
## differentiation score and for the k-means compartment assignment below. The suprabasal module
## uses GBP6 (as displayed in   and described in the Methods). An earlier
## working version of the code used DSC2 here; that was exploratory and is not the published set.
basal_markers       <- c('PDPN', 'IGFBP3', 'COL17A1', 'KRT15', 'DST')
suprabasal_markers  <- c('KRT4', 'KRT13', 'SERPINB3', 'GBP6', 'IVL')
superficial_markers <- c('KRT17', 'KRT78', 'FLG', 'CNFN', 'CRCT1')

## Markers shown in the   stacked violin (note GBP6 is displayed here).
compartment_violin_markers <- c('DST', 'COL17A1', 'KRT15', 'KRT13',
                                 'PCNA', 'MKI67', 'GBP6', 'IVL', 'CNFN', 'FLG')

## ---- QC / clustering parameters --------------------------------------------------------------
NCOUNT_FILTER   <- 3500     # min UMIs/cell applied when the EEC subset was built ( line)
UMAP_DIMS       <- 1:35     # PCs used for the EEC UMAP (Methods)
CLUSTER_K       <- 100      # FindNeighbors k.param
CLUSTER_RES     <- 0.25     # FindClusters resolution
N_DIFF_BINS     <- 10       # differentiation-score deciles


####  FUNCTIONS  ###################################################################################

## Load a single object from a .RData file under an arbitrary name.
loadRData <- function(fileName) {
  load(fileName)
  get(ls()[ls() != 'fileName'])
}

## FeaturePlot of an AddModuleScore module, with a Min/Max color bar.
module_plot <- function(obj, name, markers) {
  marker_list <- list(markers)
  obj <- AddModuleScore(obj, features = marker_list, name = 'mod_score', assay = 'RNA')
  FeaturePlot(obj, features = 'mod_score1', order = TRUE) & labs(x = 'UMAP 1', y = 'UMAP 2') &
    scale_colour_gradientn(colours = brewer.pal(n = 11, 'Spectral')[11:1], name = 'Avg. Gene\nExpression\n',
                           breaks = c(min(obj@meta.data$mod_score1), max(obj@meta.data$mod_score1)),
                           labels = c('Min', 'Max'), guide = guide_colorbar(frame.colour = 'black', ticks = FALSE)) &
    ggtitle(name, subtitle = paste(markers, collapse = ', ')) &
    theme(plot.title = element_text(face = 'bold', size = 10),
          plot.subtitle = element_text(face = 'italic', size = 8, hjust = 0.5),
          axis.text = element_blank(), axis.title = element_text(size = 8), axis.ticks = element_blank(),
          legend.title = element_text(size = 7, face = 'bold'), legend.text = element_text(size = 6.5))
}

## Default ggplot hue palette (used for cluster colors).
gg_color_hue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

## Grouped proportion barplot (mean +/- sd with jittered per-sample points).
prop_bar_cond <- function(dat, pal, groupBy, splitBy) {
  ggbarplot(dat, x = groupBy, y = 'Proportion', add = c('mean_sd', 'jitter'), fill = splitBy,
            position = position_dodge(0.8), add.params = list(color = splitBy), palette = pal) +
    scale_color_manual(values = rep('black', 3)) +
    theme(axis.title.x = element_blank(), axis.text.y = element_text(size = 8),
          plot.title = element_text(face = 'bold', hjust = 0.5))
}

## MASC: mixed-effects association of single cells (Fonseca et al., Nat Commun 2018).
## Included verbatim as an external/vendored function (logistic mixed model per cluster).
MASC <- function(dataset, cluster, contrast, random_effects = NULL, fixed_effects = NULL, verbose = FALSE) {
  # Check inputs
  cluster <- cluster[!is.na(dataset[[contrast]])]
  dataset <- dataset[!is.na(dataset[[contrast]]), ]

  if (is.factor(dataset[[contrast]]) == FALSE) {
    dataset[[contrast]] <- factor(dataset[[contrast]])
  }

  # Generate design matrix from cluster assignments
  cluster <- as.character(cluster)
  designmat <- model.matrix(~ cluster + 0, data.frame(cluster = cluster))
  dataset <- cbind(designmat, dataset)

  # Convert cluster assignments to string
  cluster <- as.character(cluster)
  # Prepend design matrix generated from cluster assignments
  designmat <- model.matrix(~ cluster + 0, data.frame(cluster = cluster))
  dataset <- cbind(designmat, dataset)
  # Create output list to hold results
  res <- vector(mode = "list", length = length(unique(cluster)))
  names(res) <- attributes(designmat)$dimnames[[2]]

  # Create model formulas
  if (!is.null(fixed_effects) && !is.null(random_effects)) {
    model_rhs <- paste0(c(paste0(fixed_effects, collapse = "+"),
                          paste0("(1|", random_effects, ")", collapse = "+")),
                        collapse = " + ")
    message(paste("Using null model:", "cluster~", model_rhs))
  } else if (is.null(fixed_effects) && !is.null(random_effects)) {
    model_rhs <- paste0("(1|", random_effects, ")", collapse = "+")
    message(paste("Using null model:", "cluster~", model_rhs))
  } else {
    model_rhs <- "1" # only includes intercept
    message(paste("Using null model:", "cluster~", model_rhs))
    stop("No random or fixed effects specified")
  }

  # Initialize list to store model objects for each cluster
  cluster_models <- vector(mode = "list",
                           length = length(attributes(designmat)$dimnames[[2]]))
  names(cluster_models) <- attributes(designmat)$dimnames[[2]]

  # Run nested mixed-effects models for each cluster
  for (i in seq_along(attributes(designmat)$dimnames[[2]])) {
    test_cluster <- attributes(designmat)$dimnames[[2]][i]
    message(paste("Creating logistic mixed models for", test_cluster))
    null_fm <- as.formula(paste0(c(paste0(test_cluster, "~1+"), model_rhs), collapse = ""))
    full_fm <- as.formula(paste0(c(paste0(test_cluster, "~", contrast, "+"), model_rhs), collapse = ""))
    # Run null and full mixed-effects models
    null_model <- lme4::glmer(formula = null_fm, data = dataset,
                              family = binomial, nAGQ = 1, verbose = 0,
                              control = glmerControl(optimizer = "bobyqa"))
    full_model <- lme4::glmer(formula = full_fm, data = dataset,
                              family = binomial, nAGQ = 1, verbose = 0,
                              control = glmerControl(optimizer = "bobyqa"))
    model_lrt <- anova(null_model, full_model)
    # calculate confidence intervals for contrast term beta
    contrast_lvl2 <- paste0(contrast, levels(dataset[[contrast]])[2])
    contrast_ci <- confint.merMod(full_model, method = "Wald", parm = contrast_lvl2)
    # Save model objects to list
    cluster_models[[i]]$null_model <- null_model
    cluster_models[[i]]$full_model <- full_model
    cluster_models[[i]]$model_lrt <- model_lrt
    cluster_models[[i]]$confint <- contrast_ci
  }

  # Organize results into output dataframe
  output <- data.frame(cluster = attributes(designmat)$dimnames[[2]], size = colSums(designmat))
  output$model.pvalue <- sapply(cluster_models, function(x) x$model_lrt[["Pr(>Chisq)"]][2])
  output[[paste(contrast_lvl2, "OR", sep = ".")]] <- sapply(cluster_models, function(x) exp(fixef(x$full)[[contrast_lvl2]]))
  output[[paste(contrast_lvl2, "OR", "95pct.ci.lower", sep = ".")]] <- sapply(cluster_models, function(x) exp(x$confint[contrast_lvl2, "2.5 %"]))
  output[[paste(contrast_lvl2, "OR", "95pct.ci.upper", sep = ".")]] <- sapply(cluster_models, function(x) exp(x$confint[contrast_lvl2, "97.5 %"]))

  # Return MASC results
  return(output)
}


####  LOAD EPITHELIAL OBJECT  ######################################################################
## Epithelial (Ep) cells only, subset from All_integratedObj.RData and re-integrated (see header).
## The object already carries the published 35-PC UMAP embedding.

epi_set <- loadRData(paste0(results_dir, '/eso_epi_integratedObj.RData'))

## PC elbow used to choose the 35 dims (Methods /  H).
pdf(file = paste0(results_dir, '/epi_PC_plot_integratedAssay.pdf'), width = 5, height = 6)
ElbowPlot(epi_set, ndims = 50)
dev.off()


####  CELL CYCLE SCORING  ##########################################################################
##  (phase). CC.diff is also regressed/used later for the "proliferating" split.

s.genes   <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
DefaultAssay(epi_set) <- 'RNA'
epi_set <- CellCycleScoring(epi_set, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
epi_set$CC.diff <- epi_set$S.Score - epi_set$G2M.Score

pdf(file = paste0(results_dir, '/Epi_UMAP_phase.pdf'), width = 7, height = 6)
DimPlot(epi_set, reduction = 'umap', label = FALSE, group.by = 'Phase', shuffle = TRUE, raster = TRUE,
        pt.size = 0.01, cols = c('#85AD51', '#59B9BF', '#A081B5')) + ggtitle('') &
  theme(plot.title = element_text(face = 'bold', size = 10),
        axis.text = element_blank(), axis.title = element_text(size = 8), axis.ticks = element_blank(),
        legend.title = element_text(size = 8, face = 'bold'), legend.text = element_text(size = 8))
dev.off()


####  QC / OVERVIEW PLOTS  #########################################################################

## UMAP overviews by condition and location.
for (grp in c('condition', 'location')) {
  pdf(file = paste0(results_dir, '/Epi_UMAP_', grp, '.pdf'), width = 7, height = 6)
  print(DimPlot(epi_set, reduction = 'umap', label = FALSE, group.by = grp, shuffle = TRUE,
                raster = TRUE, pt.size = 0.01, cols = c('#FF00FF', '#00FFFF', '#FFFF00')) + ggtitle('') &
          theme(plot.title = element_text(face = 'bold', size = 10),
                axis.text = element_blank(), axis.title = element_text(size = 8), axis.ticks = element_blank(),
                legend.title = element_text(size = 8, face = 'bold'), legend.text = element_text(size = 8)))
  dev.off()
}

## Per-sample distribution within the integrated UMAP (QC: checks no single sample dominates).
samples <- names(table(epi_set@meta.data$orig.ident))
for (sample in samples) {
  epi_set@meta.data$sampleBin <- ifelse(epi_set@meta.data$orig.ident == sample, 1, 0)
  pdf(file = paste0(results_dir, '/Epi_UMAP_SampleDistributions/Epi_dimPlot_', sample, '.pdf'), width = 5, height = 5)
  print(DimPlot(epi_set, reduction = 'umap', label = FALSE, group.by = 'sampleBin', order = TRUE,
                cols = c('grey80', 'black')) + ggtitle(sample) + theme(legend.position = 'none'))
  dev.off()
}

## UMI/cell distribution by cluster with the >3500 filter line.
pdf(file = paste0(results_dir, '/Epi_violin_nCount.pdf'), width = 8, height = 5)
print(ggplot(epi_set@meta.data, aes(x = integrated_snn_res.0.5, y = nCount_RNA, fill = integrated_snn_res.0.5)) +
        geom_violin(scale = 'width') + labs(y = 'nCount_RNA', x = 'Cluster') + ggtitle('Epithelial Cell Clusters') +
        theme_classic() +
        theme(axis.title = element_text(face = 'bold', size = 8), axis.ticks = element_blank(),
              legend.position = 'none', panel.grid.major.y = element_line(), panel.grid.minor.y = element_line()) +
        geom_hline(yintercept = NCOUNT_FILTER, color = 'red'))
dev.off()


####  COMPARTMENT-MARKER MODULE SCORES  #############################################

module_plots <- list(
  module_plot(epi_set, 'Basal',       basal_markers),
  module_plot(epi_set, 'Suprabasal',  suprabasal_markers),
  module_plot(epi_set, 'Superficial', superficial_markers))

pdf(file = paste0(results_dir, '/Epi_modulePlots.pdf'), width = 12, height = 4)
print(plot_grid(plotlist = module_plots, ncol = 3))
dev.off()


####  CLUSTERING ( H-I)  ################################################################
## SNN clustering on the integrated assay; resolution 0.25 yields the 8 EEC clusters.

plan('multisession', workers = 4)
DefaultAssay(epi_set) <- 'integrated'
epi_set <- FindNeighbors(epi_set, reduction = 'pca', dims = UMAP_DIMS, k.param = CLUSTER_K, n.trees = 50)
epi_set <- FindClusters(epi_set, resolution = CLUSTER_RES, algorithm = 1, n.iter = 10, n.start = 10)

## Re-label clusters 0-7 into a biologically ordered 1-8 (basal -> superficial). The compartment
## assignment below references these ordered cluster numbers, so this step must precede it.
reorder <- c(6, 1, 8, 7, 3, 4, 5, 2)
epi_set@meta.data$seurat_clusters <- factor(mapvalues(epi_set@meta.data$integrated_snn_res.0.25,
                                                      from = as.character(0:7), to = reorder),
                                            levels = 1:8)

pdf(file = paste0(results_dir, '/Epi_UMAP_clusters_k50_res0.25_ordered.pdf'), width = 7, height = 6)
DimPlot(epi_set, reduction = 'umap', label = TRUE, group.by = 'seurat_clusters', order = as.character(30:0)) + ggtitle('') &
  theme(plot.title = element_text(face = 'bold', size = 10),
        axis.text = element_blank(), axis.title = element_text(size = 8), axis.ticks = element_blank(),
        legend.title = element_text(size = 7, face = 'bold'), legend.text = element_text(size = 6.5))
dev.off()


####  DIFFERENTIATION SCORE  #######################################################################
## Per-cell module scores for the three compartments, combined into a single differentiation axis.

epi_set <- AddModuleScore(epi_set, features = list(basal_markers),       name = 'basal',       assay = 'RNA')
epi_set <- AddModuleScore(epi_set, features = list(suprabasal_markers),  name = 'suprabasal',  assay = 'RNA')
epi_set <- AddModuleScore(epi_set, features = list(superficial_markers), name = 'superficial', assay = 'RNA')

epi_set$diff.score <- epi_set$superficial1 + epi_set$suprabasal1 - epi_set$basal1

pdf(file = paste0(results_dir, '/Epi_UMAP_diff.score.pdf'), width = 7, height = 7)
print(FeaturePlot(epi_set, features = 'diff.score', order = TRUE) & labs(x = 'UMAP 1', y = 'UMAP 2') &
        scale_colour_gradientn(colours = brewer.pal(n = 11, 'Spectral')[11:1], name = 'Differentiation\nScore\n',
                               breaks = c(min(epi_set@meta.data$diff.score), max(epi_set@meta.data$diff.score)),
                               labels = c('Min', 'Max'), guide = guide_colorbar(frame.colour = 'black', ticks = FALSE)) &
        ggtitle('Differentiation Score') &
        theme(plot.title = element_text(face = 'bold', size = 10),
              axis.text = element_blank(), axis.title = element_text(size = 8), axis.ticks = element_blank(),
              legend.title = element_text(size = 7, face = 'bold'), legend.text = element_text(size = 6.5)))
dev.off()


####  COMPARTMENT ASSIGNMENT  #############################################################
## Cells are assigned to basal / suprabasal / superficial using the ordered clusters plus a k-means
## split (on basal+suprabasal markers) of the ambiguous proliferating clusters, then refined for the
## superficial/suprabasal boundary. Finally the basal and suprabasal compartments are split into
## cycling ("proliferating") vs non-cycling using cell-cycle phase, giving the 5 compartments.

DefaultAssay(epi_set) <- 'RNA'

## Split the proliferating clusters (3,4,5) along the basal<->suprabasal axis.
rep_subset <- subset(epi_set, seurat_clusters %in% c(3, 4, 5))
DefaultAssay(rep_subset) <- 'RNA'
layer_km <- kmeans(FetchData(rep_subset, vars = c(basal_markers, suprabasal_markers)), 2,
                   iter.max = 1000, nstart = 10)
rep_subset@meta.data$layer_cluster <- layer_km$cluster

## Initial 3-compartment assignment.
epi_set$layers <- ifelse(epi_set@meta.data$seurat_clusters %in% c(1, 2), 'basal',
                  ifelse(epi_set@meta.data$seurat_clusters %in% c(8), 'superficial',
                  ifelse(epi_set@meta.data$seurat_clusters %in% c(6, 7), 'suprabasal',
                  ifelse(rownames(epi_set@meta.data) %in%
                           rownames(rep_subset@meta.data[rep_subset@meta.data$layer_cluster == 2, ]),
                         'basal', 'suprabasal'))))

## Refine the superficial<->suprabasal boundary (clusters 6,7,8) with a k-means on module scores.
superficial_subset <- subset(epi_set, seurat_clusters %in% c(6, 7, 8))
layer_km <- kmeans(superficial_subset@meta.data[, c('superficial1', 'suprabasal1')], 2,
                   iter.max = 1000, nstart = 10)
superficial_subset@meta.data$layer_cluster <- layer_km$cluster

epi_set$layers <- ifelse(rownames(epi_set@meta.data) %in%
                           rownames(superficial_subset@meta.data[superficial_subset@meta.data$layer_cluster == 1, ]),
                         'superficial', epi_set$layers)
epi_set$layers <- factor(epi_set$layers, levels = c('basal', 'suprabasal', 'superficial'))

## Split basal/suprabasal into cycling vs non-cycling -> 5 compartments (layers_rep).
epi_set$layers_rep <- factor(
  ifelse((epi_set$layers == 'basal')      & (epi_set$Phase %in% c('S', 'G2M')), 'replicating basal',
  ifelse((epi_set$layers == 'suprabasal') & (epi_set$Phase %in% c('S', 'G2M')), 'replicating suprabasal',
  ifelse(epi_set$layers == 'basal', 'basal',
  ifelse(epi_set$layers == 'suprabasal', 'suprabasal', 'superficial')))),
  levels = c('basal', 'replicating basal', 'replicating suprabasal', 'suprabasal', 'superficial'))

## : EEC UMAP colored by the 5 compartments.
pdf(file = paste0(results_dir, '/Epi_UMAP_layersRep_kmeans.pdf'), width = 8, height = 6)
DimPlot(epi_set, reduction = 'umap', label = FALSE, group.by = 'layers_rep', raster = FALSE,
        cols = layer_cols_5) + ggtitle('') &
  theme(plot.title = element_text(face = 'bold', size = 10),
        axis.text = element_blank(), axis.title = element_text(size = 8), axis.ticks = element_blank(),
        legend.title = element_text(size = 7, face = 'bold'), legend.text = element_text(size = 12))
dev.off()

## Also save the 3-compartment UMAP (used as an intermediate / supplement panel).
pdf(file = paste0(results_dir, '/Epi_UMAP_layers_kmeans.pdf'), width = 7, height = 6)
DimPlot(epi_set, reduction = 'umap', label = FALSE, group.by = 'layers', cols = layer_cols_3) + ggtitle('') &
  theme(plot.title = element_text(face = 'bold', size = 10),
        axis.text = element_blank(), axis.title = element_text(size = 8), axis.ticks = element_blank(),
        legend.title = element_text(size = 7, face = 'bold'), legend.text = element_text(size = 8))
dev.off()

## Differentiation-score deciles (used for binned comparisons / DEGs).
epi_set$diff.bin <- cut(epi_set$diff.score,
                        breaks = unique(quantile(epi_set$diff.score, probs = seq.int(0, 1, by = 1 / N_DIFF_BINS))),
                        include.lowest = TRUE)


####  COMPARTMENT MARKER STACKED VIOLIN  ##################################################

DefaultAssay(epi_set) <- 'RNA'
pdf(file = paste0(results_dir, '/Epi_layers_markers_violin.pdf'), width = 3, height = 8)
print(Stacked_VlnPlot(epi_set, compartment_violin_markers, pt.size = 0, group.by = 'layers_rep',
                      raster = FALSE, plot_spacing = 0, colors = layer_cols_5, x_lab_rotate = TRUE) &
        theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()))
dev.off()

## Save the object now that it carries cell-cycle scores, compartment labels, and the diff score.
save(epi_set, file = paste0(results_dir, '/eso_epi_integratedObj.RData'))


####  COMPARTMENT PROPORTIONS BY CONDITION/LOCATION  ######################################

## Per-sample 3-compartment proportions.
epi_layer_df <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(epi_layer_df) <- c('Sample', 'Location', 'Condition', 'Layer', 'Proportion')
for (sample in samples) {
  dat   <- epi_set@meta.data[epi_set@meta.data$orig.ident == sample, ]
  loc_s <- ifelse(substr(sample, nchar(sample), nchar(sample)) == 'P', 'Proximal', 'Distal')
  c_s   <- dat[1, cond]
  total <- nrow(dat)
  for (l in c('basal', 'suprabasal', 'superficial')) {
    epi_layer_df <- rbind(epi_layer_df,
                          data.frame(Sample = sample, Location = loc_s, Condition = as.character(c_s),
                                     Layer = l, Proportion = round(nrow(dat[dat$layers == l, ]) / total, 2)))
  }
}
epi_layer_df$Condition <- factor(epi_layer_df$Condition, levels = conditions)
epi_layer_df$Layer     <- capitalize(epi_layer_df$Layer)
write.table(epi_layer_df, paste0(results_dir, '/epi_layer_props.tsv'),
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')

## Per-sample 5-compartment proportions.
epi_layerRep_df <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(epi_layerRep_df) <- c('Sample', 'Location', 'Condition', 'LayerRep', 'Proportion')
for (sample in samples) {
  dat   <- epi_set@meta.data[epi_set@meta.data$orig.ident == sample, ]
  loc_s <- ifelse(substr(sample, nchar(sample), nchar(sample)) == 'P', 'Proximal', 'Distal')
  c_s   <- dat[1, cond]
  total <- nrow(dat)
  for (l in names(table(dat$layers_rep))) {
    epi_layerRep_df <- rbind(epi_layerRep_df,
                             data.frame(Sample = sample, Location = loc_s, Condition = as.character(c_s),
                                        LayerRep = l, Proportion = round(nrow(dat[dat$layers_rep == l, ]) / total, 2)))
  }
}
epi_layerRep_df$Condition <- factor(epi_layerRep_df$Condition, levels = conditions)
epi_layerRep_df$LayerRep  <- gsub(' ', '\n', capitalize(epi_layerRep_df$LayerRep))
write.table(epi_layerRep_df, paste0(results_dir, '/epi_layerRep_props.tsv'),
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')

## Per-sample differentiation-score summary (consumed by 11_clinical_associations.R).
diffScoreStats <- aggregate(diff.score ~ orig.ident, data = epi_set@meta.data, FUN = mean)
names(diffScoreStats)[2] <- 'mean_diff_score'
diffScoreStats$median_diff_score <- aggregate(diff.score ~ orig.ident, data = epi_set@meta.data, FUN = median)$diff.score
write.table(diffScoreStats, paste0(results_dir, '/diffScoreStats.tsv'),
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')

## Per-sample cluster proportions (supplementary).
epi_cluster_df <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(epi_cluster_df) <- c('Sample', 'Location', 'Condition', 'Cluster', 'nCells', 'totCells', 'Proportion')
for (sample in samples) {
  dat   <- epi_set@meta.data[epi_set@meta.data$orig.ident == sample, ]
  loc_s <- ifelse(substr(sample, nchar(sample), nchar(sample)) == 'P', 'Proximal', 'Distal')
  c_s   <- dat[1, cond]
  total <- nrow(dat)
  for (l in names(table(dat$seurat_clusters))) {
    nCells <- nrow(dat[dat$seurat_clusters == l, ])
    epi_cluster_df <- rbind(epi_cluster_df,
                            data.frame(Sample = sample, Location = loc_s, Condition = as.character(c_s),
                                       Cluster = l, nCells = nCells, totCells = total,
                                       Proportion = round(nCells / total, 2)))
  }
}
write.table(epi_cluster_df, paste0(results_dir, '/epi_cluster_props.tsv'),
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')

## : proportion barplots, split by location.
p1 <- prop_bar_cond(epi_layerRep_df[epi_layerRep_df$Location == 'Proximal', ], cond_cols[3:1], 'LayerRep', 'Condition') + ggtitle('Proximal')
p2 <- prop_bar_cond(epi_layerRep_df[epi_layerRep_df$Location == 'Distal', ],   cond_cols[3:1], 'LayerRep', 'Condition') + ggtitle('Distal')
legend  <- get_legend(p1 + theme(legend.box.margin = margin(0, 0, 0, 12), legend.position = 'right',
                                 legend.title = element_text(face = 'bold')))
p_bars  <- plot_grid(p1 + theme(legend.position = 'none', plot.margin = unit(c(0, 0, 1, 0), 'lines')),
                     p2 + theme(legend.position = 'none', plot.margin = unit(c(1, 0, 0, 0), 'lines')), ncol = 1)
pdf(file = paste0(results_dir, '/Epi_layerRepProp_vs_condition_Region.pdf'), width = 7, height = 7)
plot_grid(p_bars, legend, rel_widths = c(1, 0.2), scale = 0.95)
dev.off()

## Proportion significance: MASC (mixed-effects, sample as random effect) + propeller, per
## pairwise contrast, overall and within each location, FDR-corrected within each location stratum.
## Drives the significance annotations in .
test_df <- epi_set@meta.data[, c('orig.ident', 'location', 'condition', 'layers', 'layers_rep')]
test_df$layers_rep <- gsub(' ', '_', test_df$layers_rep)
test_df$ssc_hc   <- ifelse(test_df$condition == 'SSc',  1, ifelse(test_df$condition == 'HC',   0, NA))
test_df$gerd_hc  <- ifelse(test_df$condition == 'GERD', 1, ifelse(test_df$condition == 'HC',   0, NA))
test_df$ssc_gerd <- ifelse(test_df$condition == 'SSc',  1, ifelse(test_df$condition == 'GERD', 0, NA))

prop_comps <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(prop_comps) <- c('cluster', 'size', 'masc.p', 'OR', 'OR.95.ci.lower', 'OR.95.ci.upper',
                          'comp', 'location', 'propeller.p')
for (comp in c('ssc_hc', 'gerd_hc', 'ssc_gerd')) {
  for (l in c('All', 'Proximal', 'Distal')) {
    temp_df <- if (l == 'All') test_df else test_df[test_df$location == l, ]
    masc <- MASC(temp_df, cluster = temp_df$layers_rep, contrast = comp, random_effects = 'orig.ident')
    masc$cluster  <- gsub('cluster', '', masc$cluster)
    masc$comp     <- comp
    masc$location <- l
    propel <- propeller(clusters = temp_df$layers_rep, sample = temp_df$orig.ident, group = temp_df[, comp])
    masc <- merge(masc, propel[, c('BaselineProp.clusters', 'P.Value')],
                  by.x = 'cluster', by.y = 'BaselineProp.clusters')
    prop_comps <- rbind(prop_comps, setNames(masc, names(prop_comps)))
  }
}
prop_comps <- prop_comps[, c(7, 8, 1, 2, 4:6, 3, 9)]
prop_comps[prop_comps$location == 'All', 'masc.p.fdr']  <- p.adjust(prop_comps[prop_comps$location == 'All', 'masc.p'], method = 'fdr')
prop_comps[prop_comps$location != 'All', 'masc.p.fdr']  <- p.adjust(prop_comps[prop_comps$location != 'All', 'masc.p'], method = 'fdr')
prop_comps[prop_comps$location == 'All', 'propeller.p.fdr'] <- p.adjust(prop_comps[prop_comps$location == 'All', 'propeller.p'], method = 'fdr')
prop_comps[prop_comps$location != 'All', 'propeller.p.fdr'] <- p.adjust(prop_comps[prop_comps$location != 'All', 'propeller.p'], method = 'fdr')
write.table(prop_comps, paste0(results_dir, '/epi_layerRep_props_stats.tsv'),
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')


####  PER-COMPARTMENT SINGLE-CELL DEGs  ############################################################
## Wilcoxon DEGs per 3-level compartment for each pairwise contrast. Saved for downstream GSEA / TF
## enrichment (07_gsea_enrichment.R). (The original loaded these list objects from disk;
## the computation that produces them is kept here and the results are saved with standard names.)

plan('multisession', workers = 2)
DefaultAssay(epi_set) <- 'RNA'
Idents(epi_set) <- 'condition'

layers <- c('basal', 'suprabasal', 'superficial')
deg_list_ssc_hc <- setNames(lapply(layers, function(l) {
  FindMarkers(subset(epi_set, layers == l), ident.1 = 'SSc', ident.2 = 'HC', logfc.threshold = 0.05, min.pct = 0.01) }), layers)
deg_list_ssc_gerd <- setNames(lapply(layers, function(l) {
  FindMarkers(subset(epi_set, layers == l), ident.1 = 'SSc', ident.2 = 'GERD', logfc.threshold = 0.05, min.pct = 0.01) }), layers)
deg_list_gerd_hc <- setNames(lapply(layers, function(l) {
  FindMarkers(subset(epi_set, layers == l), ident.1 = 'GERD', ident.2 = 'HC', logfc.threshold = 0.05, min.pct = 0.01) }), layers)

save(deg_list_ssc_hc,   file = paste0(results_dir, '/epi_degs_ssc_hc_layers.RData'))
save(deg_list_gerd_hc,  file = paste0(results_dir, '/epi_degs_gerd_hc_layers.RData'))
save(deg_list_ssc_gerd, file = paste0(results_dir, '/epi_degs_ssc_gerd_layers.RData'))

## Per-compartment (5-level layers_rep) DEGs, computed separately within each biopsy location.
## Tested broadly (no logFC threshold, min.pct = 0.01) so the downstream GSEA ranking statistic
## (-log10 p * |log2FC|) is well populated. These feed 07_gsea_enrichment.R.
layersRep <- levels(epi_set$layers_rep)
for (region in c('Proximal', 'Distal')) {
  loc_set <- subset(epi_set, location == region)
  deg_list_ssc_hc <- setNames(lapply(layersRep, function(l) {
    FindMarkers(subset(loc_set, layers_rep == l), ident.1 = 'SSc', ident.2 = 'HC', logfc.threshold = 0, min.pct = 0.01) }), layersRep)
  deg_list_gerd_hc <- setNames(lapply(layersRep, function(l) {
    FindMarkers(subset(loc_set, layers_rep == l), ident.1 = 'GERD', ident.2 = 'HC', logfc.threshold = 0, min.pct = 0.01) }), layersRep)
  deg_list_ssc_gerd <- setNames(lapply(layersRep, function(l) {
    FindMarkers(subset(loc_set, layers_rep == l), ident.1 = 'SSc', ident.2 = 'GERD', logfc.threshold = 0, min.pct = 0.01) }), layersRep)
  save(deg_list_ssc_hc,   file = paste0(results_dir, '/epi_degs_ssc_hc_layersRep_', region, '.RData'))
  save(deg_list_gerd_hc,  file = paste0(results_dir, '/epi_degs_gerd_hc_layersRep_', region, '.RData'))
  save(deg_list_ssc_gerd, file = paste0(results_dir, '/epi_degs_ssc_gerd_layersRep_', region, '.RData'))
}
