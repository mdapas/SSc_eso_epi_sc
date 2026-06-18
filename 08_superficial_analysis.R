####################################################################################################
#  08_superficial_analysis.R
#
#  Subclustering and characterization of the superficial esophageal epithelial cells (EECs).
#
#  Paper:
#   Dapas et al. (2025) "Cellular and molecular dysregulation of the esophageal epithelium in
#   systemic sclerosis." JCI Insight.
#
#  Purpose / figures produced:
#   - UMAP of the superficial compartment, colored by subcluster (5 clusters)
#   - Stacked violin of subcluster markers (FLG, CTSV, MT1H, GJB6, FGFBP1)
#   - Superficial module score by subcluster
#   - Subcluster proportions by condition, with MASC significance
#   - Metallothionein module score on the UMAP (Proximal: GERD&HC vs SSc)
#   - Metallothionein module score by subcluster and condition
#   - Metallothionein module score by biopsy location and condition (HC vs SSc)
#   - resolution-scan clustree
#   - subcluster-marker + superficial-module feature plots
#   - transcription-factor (IRF1, MYC, E2F4, NFE2L2) and their target-gene module
#                    scores (feature plots split by condition, and violins by subcluster), plus the
#                    TF-vs-target expression correlations reported in the text
#   - SSc-associated immune genes (CD44, CD74, C10orf99, HBEGF, HLA-B, ...) in
#                    subcluster 5 vs subclusters 1-4
#   ( and 4H are spatial panels produced by 10_cosmx_spatial.R, not here.)
#
#  Pipeline position:
#   05_epithelial_compartments.R (assigns compartments; saves eso_epi_integratedObj.RData with layers)
#     -> 08_superficial_analysis.R (THIS SCRIPT)
#
#  Inputs:
#   - <results_dir>/eso_epi_integratedObj.RData   EEC object carrying the compartment labels
#                                                 (meta.data$layers / layers_rep) and diff.score.
#   - <results_dir>/pb_deg_res_layerRep_byRegion.RData  (06_pseudobulk_dgea.R) ct.genes + by-location
#                                                 pseudobulk DEGs, used by the FLI1 module below.
#   - <pathway_maps_dir>/ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.map   TF target map (FLI1).
#
#  Outputs (to <results_dir> / <suppdata_dir>):
#   - eso_epi_superficial.RData                   superficial subset, subclustered (5 clusters)
#   - epiSuper_clusters_props_stats.tsv           subcluster proportion statistics (MASC)
#   - epiSuper_cluster_props.tsv                  per-sample subcluster proportions (for 11)
#   - epi_super_TF_target_correlations.tsv        TF-vs-target expression correlations
#   - superficial_umap_coords.tsv, superficial_marker_expression.tsv, superficial_module_score.tsv, metallothionein_umap.tsv, metallothionein_by_subcluster.tsv, superficial_marker_module_umap.tsv, tf_<tf>_targets.tsv, cluster5_immune_genes.tsv, fli1_target_module.tsv  source data
#   - assorted .pdf figure panels
#
#  Usage:
#   Rscript 08_superficial_analysis.R   (run from the repository root)
#
#  Notes:
#   - The superficial subset is created here by subsetting the 'superficial' compartment from the
#     epithelial object and re-running PCA/UMAP/clustering on the existing integrated assay
#     (no new integration). In the original working code this subset was loaded pre-made; the
#     creation step is folded in here so the script is self-contained.
#   - Removed from the original working file: a wholesale copy of the epi_analysis.R compartment
#     analysis and its Venn / diff-score-bin / volcano / pseudotime / Monocle tail (none of which
#     belong in the superficial script); FOXM1 exploration; and numerous broken/duplicate scratch
#     plots (e.g. an undefined super_set.prox used before assignment, an undefined pb_super_obj,
#     and a diff.score plot referencing an undefined loc).
#   - FLAG:  (FLI1 target-gene module) is computed here because it operates on the
#     superficial cells, but it depends on 06_pseudobulk_dgea.R outputs (the ENCODE/ChEA TF map and the
#     by-location pseudobulk DEGs) and FLI1 is otherwise part of the endothelial . The
#     reorganization pass may relocate it. It references the superficial pseudobulk comparison by
#     index (qlf.res[[30]]); confirm that is the intended contrast.
#   - The subcluster-proportion MASC is computed within each biopsy location
#     (Proximal and Distal).
####################################################################################################


####  LOAD LIBRARIES  ##############################################################################

suppressPackageStartupMessages({
  library(methods)
  library(Seurat)
  library(ggplot2)
  library(ggpubr)
  library(scCustomize)   # FeaturePlot_scCustom (split feature plots)
  library(RColorBrewer)
  library(cowplot)
  library(patchwork)
  library(vioplot)
  library(scales)
  library(plyr)
  library(dplyr)
  library(stringr)
  library(future)
  library(lme4)          # MASC mixed-model proportion test
  library(speckle)       # propeller()
})


####  SET PARAMETERS  ##############################################################################

set.seed(0123456789)
plan('multisession', workers = 4)

## ---- Paths (relative to repository root) ----
results_dir  <- './results'
suppdata_dir <- './results/supporting_data'
pathway_maps_dir <- './results/pathway_gene_maps'   # Enrichr gene-set maps (shared with 06_pseudobulk_dgea.R)

## ---- Colors ----
cond_cols  <- c('#8156B3', '#e26b53', '#A8A39d')                          # SSc, GERD, HC
super_cols <- c('#1F9E89', '#E76254', '#67AD5B', '#E3B23C', '#9C8CC4')    # subclusters 1-5 (subcluster UMAP palette)

## ---- Marker / module gene sets ----
superficial_markers <- c('KRT17', 'KRT78', 'FLG', 'CNFN', 'CRCT1')
superficial_marker_genes       <- c('FLG', 'CTSV', 'MT1H', 'GJB6', 'FGFBP1')
metallo_genes       <- c('MT1A', 'MT1E', 'MT1F', 'MT1G', 'MT1H', 'MT1M', 'MT1X', 'MT2A')
immune_genes        <- c('CD44', 'CD74', 'C10orf99', 'HBEGF', 'HLA-B',
                         'PSMB8', 'PSMB9', 'TAP1', 'SGK1', 'LTB4R', 'CLEC2B', 'BST2')

## Transcription-factor target-gene modules. The IRF1 set is the small curated
## list used originally; the MYC / E2F4 / NFE2L2 sets are the MSigDB/Enrichr target-gene sets.
irf1_target_genes  <- c('RAB10', 'PPP2R2A', 'LCMT1', 'SKA2', 'ZDHHC14')
myc_target_genes   <- c('SRM','MRTO4','HNRNPR','SRSF10','ZNF593','RPA2','SFPQ','YBX1','EBNA1BP2','SERBP1','RPL5','PRPF38B','CKS1B','GPATCH4','IVNS1ABP','PARP1','ARV1','RPS27A','MPHOSPH10','CCT7','EIF5B','ZC3H8','NIFK','NOP58','EEF1B2','FASTKD2','NCL','PTMA','THUMPD3','SEC13','RPSA','UMPS','TFRC','PAICS','OSTC','GAR1','ABCE1','DCTD','CCT5','BRIX1','RPS23','LYRM7','TCERG1','NPM1','HNRNPAB','EEF1E1','MRPS18B','TUBB','RPF2','GJA1','RPS12','AIMP2','BZW2','TBRG4','CCT6A','TBL2','PRKDC','MRPS28','MTDH','DCAF13','IARS','MRPL50','SET','CDC123','DDX21','NOLC1','RPLP2','CHID1','MRPL23','IPO7','SLC3A2','BANF1','CCND1','RPS3','CLNS1A','CEP57','STT3A','COPS7A','EMG1','MRPS35','PA2G4','SHMT2','PWP1','NAA25','RAN','TFDP1','ATXN3','INO80','PKM','RPS2','RNPS1','NUBP1','RSL1D1','RPS15A','DCTPP1','GOT2','NOB1','USP10','TSR1','SLC25A11','EIF5A','TP53','RAB34','TLCD1','AATF','MLLT6','RPL19','RPL27','NME1','NME2','USP14','NARS','PTBP1','DUS3L','ILF3','TNPO2','DDX39A','ATP13A1','URI1','TIMM50','SNRPA','RPS19','ZNF428','GRWD1','BAX','PRMT1','PRPF31','HSPBP1','PCNA','RPS21','SOD1','PSMG1','RANBP1','EWSR1','TXN2','RANGAP1','RPL36A','RBMX','DKC1','FKBP11','YDJC')
e2f4_target_genes  <- c('SRSF10','RPA2','RBBP4','YARS','MAGOH','USP1','ITGB3BP','SERBP1','SRSF11','ANP32E','GABPB2','SHC1','CKS1B','DHX9','TPR','UBE2T','NUCKS1','NSL1','TATDN3','CENPF','ADSS','RRM2','SMC6','SRSF7','POLE4','RIF1','PRPF40A','ARL6IP6','HAT1','NOP58','PTMA','DTYMK','XPC','MBD4','TOPBP1','GMPS','PAICS','LARP7','MAD2L1','MSH3','TCERG1','NPM1','HNRNPAB','DEK','GMNN','DDX39B','OARD1','CENPW','MEPCE','CALU','PRKDC','TERF1','ATAD2','MYC','CHRAC1','NOL8','SMC2','GTF3C5','ZMYND19','MASTL','CDK1','DNAJC9','KIF20B','NOLC1','NDUFS3','SSRP1','ZFP91','SNHG1','SLC3A2','SAC3D1','ANAPC15','PRCP','CEP57','H2AFX','MRPL51','USP5','TUBA1B','CBX5','YEATS4','NUP37','ANAPC5','HMGB1','RNASEH2B','CKAP2','MZT1','TFDP1','SUPT16H','LRR1','TMX1','CDKN3','MTHFD1','EFCAB11','UBR7','SIVA1','NUSAP1','DUT','MRPL46','MRPS11','PRC1','RNPS1','NFATC2IP','CFAP20','CKLF','EXOSC6','CMC2','RPA1','TMEM107','COPS3','TOP2A','MRPL27','RAD51C','TK1','SNRPD1','PTBP1','ALDH16A1','FAM110A','PCNA','SNX5','RANBP1','MCM5','LMF2','SMC1A')
nfe2l2_target_genes <- c('KIF1B','PRDX1','F3','GOLPH3L','PRDX6','SLC30A1','SRP9','RBKS','YPEL5','STRN','MEIS1','RAB3GAP1','MMADHC','RND3','MARCH7','ZFAND2B','MKRN2','ABHD5','PRKCD','SLMAP','TBC1D23','SSR3','TBL1XR1','TBC1D14','ARAP2','DAPP1','GOLPH3','RAI14','RICTOR','ARHGEF28','PAM','CTNNA1','PXDC1','MARCKS','HINT3','PTPRK','NHSL1','HECA','CITED2','VTA1','SASH1','SYNJ2','MAP3K4','C1GALT1','SNX13','TAX1BP1','YWHAG','BUD31','HBP1','MDFIC','DUSP4','SARAF','SLC20A2','ZNF706','OXR1','EMC2','EIF3H','NSMCE2','ANXA1','MEGF9','FNBP1','RXRA','ABI1','MAPK8','IPMK','ZNF365','MICU1','ADK','SMNDC1','RNF141','EHF','CCDC90B','MAML2','LAYN','SIK2','UBE4A','ARHGAP32','PLEKHG6','OSBPL8','PAWR','APPL2','ABCB9','ZDHHC20','UBL3','KLF5','LMO7','FBXL3','ARHGAP5','SPTSSA','FRMD6','TTC9','MAP3K9','ZNF410','CPPED1','TMEM159','AKTIP','SUMO2','MAFG','TGIF1','ZNF519','BLVRB','FTL','OSER1','PPDPF','PACSIN2','APOOL','GCNT1','TBCEL')

## ---- Clustering parameters ----
RES_SCAN  <- seq(0.05, 0.5, by = 0.05)   # resolution scan for the  clustree
CLUSTER_RES <- 0.25                       # resolution giving the 5 superficial subclusters


####  FUNCTIONS  ###################################################################################

loadRData <- function(fileName) {
  load(fileName)
  get(ls()[ls() != 'fileName'])
}

## AddModuleScore FeaturePlot with a Min/Max color bar.
module_plot <- function(obj, name, markers) {
  obj <- AddModuleScore(obj, features = list(markers), name = 'mod_score', assay = 'RNA')
  FeaturePlot(obj, features = 'mod_score1', order = TRUE) & labs(x = 'UMAP 1', y = 'UMAP 2') &
    scale_colour_gradientn(colours = brewer.pal(n = 11, 'Spectral')[11:1], name = 'Avg. Gene\nExpression\n',
                           breaks = c(min(obj$mod_score1), max(obj$mod_score1)), labels = c('Min', 'Max'),
                           guide = guide_colorbar(frame.colour = 'black', ticks = FALSE)) &
    ggtitle(name, subtitle = paste(markers, collapse = ', ')) &
    theme(plot.title = element_text(face = 'bold', size = 10),
          plot.subtitle = element_text(face = 'italic', size = 8, hjust = 0.5),
          axis.text = element_blank(), axis.title = element_text(size = 8), axis.ticks = element_blank())
}

## MASC: mixed-effects association of single cells (Fonseca et al., Nat Commun 2018), included
## verbatim as a vendored external function.
MASC <- function(dataset, cluster, contrast, random_effects = NULL, fixed_effects = NULL, verbose = FALSE) {
  cluster <- cluster[!is.na(dataset[[contrast]])]
  dataset <- dataset[!is.na(dataset[[contrast]]), ]
  if (is.factor(dataset[[contrast]]) == FALSE) dataset[[contrast]] <- factor(dataset[[contrast]])

  cluster <- as.character(cluster)
  designmat <- model.matrix(~ cluster + 0, data.frame(cluster = cluster))
  dataset <- cbind(designmat, dataset)

  cluster <- as.character(cluster)
  designmat <- model.matrix(~ cluster + 0, data.frame(cluster = cluster))
  dataset <- cbind(designmat, dataset)
  res <- vector(mode = "list", length = length(unique(cluster)))
  names(res) <- attributes(designmat)$dimnames[[2]]

  if (!is.null(fixed_effects) && !is.null(random_effects)) {
    model_rhs <- paste0(c(paste0(fixed_effects, collapse = "+"),
                          paste0("(1|", random_effects, ")", collapse = "+")), collapse = " + ")
    message(paste("Using null model:", "cluster~", model_rhs))
  } else if (is.null(fixed_effects) && !is.null(random_effects)) {
    model_rhs <- paste0("(1|", random_effects, ")", collapse = "+")
    message(paste("Using null model:", "cluster~", model_rhs))
  } else {
    model_rhs <- "1"
    message(paste("Using null model:", "cluster~", model_rhs))
    stop("No random or fixed effects specified")
  }

  cluster_models <- vector(mode = "list", length = length(attributes(designmat)$dimnames[[2]]))
  names(cluster_models) <- attributes(designmat)$dimnames[[2]]

  for (i in seq_along(attributes(designmat)$dimnames[[2]])) {
    test_cluster <- attributes(designmat)$dimnames[[2]][i]
    message(paste("Creating logistic mixed models for", test_cluster))
    null_fm <- as.formula(paste0(c(paste0(test_cluster, "~1+"), model_rhs), collapse = ""))
    full_fm <- as.formula(paste0(c(paste0(test_cluster, "~", contrast, "+"), model_rhs), collapse = ""))
    null_model <- lme4::glmer(formula = null_fm, data = dataset, family = binomial, nAGQ = 1, verbose = 0,
                              control = glmerControl(optimizer = "bobyqa"))
    full_model <- lme4::glmer(formula = full_fm, data = dataset, family = binomial, nAGQ = 1, verbose = 0,
                              control = glmerControl(optimizer = "bobyqa"))
    model_lrt <- anova(null_model, full_model)
    contrast_lvl2 <- paste0(contrast, levels(dataset[[contrast]])[2])
    contrast_ci <- confint.merMod(full_model, method = "Wald", parm = contrast_lvl2)
    cluster_models[[i]]$null_model <- null_model
    cluster_models[[i]]$full_model <- full_model
    cluster_models[[i]]$model_lrt <- model_lrt
    cluster_models[[i]]$confint <- contrast_ci
  }

  output <- data.frame(cluster = attributes(designmat)$dimnames[[2]], size = colSums(designmat))
  output$model.pvalue <- sapply(cluster_models, function(x) x$model_lrt[["Pr(>Chisq)"]][2])
  output[[paste(contrast_lvl2, "OR", sep = ".")]] <- sapply(cluster_models, function(x) exp(fixef(x$full)[[contrast_lvl2]]))
  output[[paste(contrast_lvl2, "OR", "95pct.ci.lower", sep = ".")]] <- sapply(cluster_models, function(x) exp(x$confint[contrast_lvl2, "2.5 %"]))
  output[[paste(contrast_lvl2, "OR", "95pct.ci.upper", sep = ".")]] <- sapply(cluster_models, function(x) exp(x$confint[contrast_lvl2, "97.5 %"]))
  return(output)
}

## Feature plots of a gene/module split by condition (Proximal cells), used for .
tf_featurePlot <- function(obj, feats, cols) {
  suppressWarnings(lapply(feats, function(feat) {
    FeaturePlot_scCustom(obj, features = feat, order = TRUE, split.by = 'condition', na_cutoff = NA,
                         split_collect = TRUE, pt.size = 1, combine = FALSE, colors_use = cols) &
      labs(x = 'UMAP 1', y = 'UMAP 2') &
      theme(plot.title = element_text(face = 'bold', size = 10), axis.text = element_blank(),
            axis.title = element_text(size = 8), axis.ticks = element_blank(),
            legend.title = element_blank(), legend.text = element_text(size = 6.5))
  }))
}

## Violin of a gene's expression in subcluster 5 vs subclusters 1-4, by condition.
gene_cluster5vioPlot <- function(obj, feat) {
  pd <- FetchData(obj, vars = c(feat, 'seurat_clusters', 'condition'), layer = 'data')
  pd$cluster <- ifelse(pd$seurat_clusters == 5, '5', '1-4')
  ggplot(pd, aes(x = condition, y = .data[[feat]], fill = condition)) +
    geom_violin(trim = TRUE, scale = 'width', draw_quantiles = c(0.25, 0.5, 0.75)) +
    facet_wrap(~cluster, nrow = 1) + scale_fill_manual(values = rev(cond_cols)) + theme_bw() +
    ggtitle(paste(feat, '- Superficial, Proximal')) +
    theme(strip.text = element_text(face = 'bold'), axis.title.x = element_blank(),
          axis.text.x = element_blank(), panel.grid.major.x = element_blank(), axis.ticks = element_blank())
}


####  CREATE SUPERFICIAL SUBSET  ###################################################################
## Subset the superficial compartment from the epithelial object and re-embed/recluster on the
## existing integrated assay (no new integration).

epi_set   <- loadRData(paste0(results_dir, '/eso_epi_integratedObj.RData'))
super_set <- subset(epi_set, layers == 'superficial')
rm(epi_set)

DefaultAssay(super_set) <- 'integrated'
super_set <- RunPCA(super_set, assay = 'integrated')
super_set <- RunUMAP(super_set, assay = 'integrated', dims = 1:10, n.neighbors = 15L, min.dist = 0.3,
                     spread = 1, repulsion.strength = 1, negative.sample.rate = 10L, n.epochs = 1000)
super_set <- FindNeighbors(super_set, reduction = 'pca', dims = 1:35, k.param = 100, n.trees = 50)

## Resolution scan ( clustree source) + the chosen resolution.
for (r in RES_SCAN) super_set <- FindClusters(super_set, resolution = r, algorithm = 1, n.iter = 10, n.start = 10)
super_set$seurat_clusters <- factor(as.numeric(super_set[[paste0('integrated_snn_res.', CLUSTER_RES)]][, 1]), levels = 1:5)

save(super_set, file = paste0(results_dir, '/eso_epi_superficial.RData'))


####  MODULE SCORES  ###############################################################################

DefaultAssay(super_set) <- 'RNA'
super_set <- AddModuleScore(super_set, features = list(superficial_markers), name = 'superficial',   assay = 'RNA')
super_set <- AddModuleScore(super_set, features = list(metallo_genes),       name = 'metal_score',    assay = 'RNA')
super_set <- AddModuleScore(super_set, features = list(irf1_target_genes),   name = 'irf1_target_genes',   assay = 'RNA')
super_set <- AddModuleScore(super_set, features = list(myc_target_genes),    name = 'myc_target_genes',    assay = 'RNA')
super_set <- AddModuleScore(super_set, features = list(e2f4_target_genes),   name = 'e2f4_target_genes',   assay = 'RNA')
super_set <- AddModuleScore(super_set, features = list(nfe2l2_target_genes), name = 'nfe2l2_target_genes', assay = 'RNA')

super_set$condition <- factor(super_set$condition, levels = c('HC', 'GERD', 'SSc'))
super_set.prox <- subset(super_set, location == 'Proximal')


####    ############################################################################

pdf(file = paste0(results_dir, '/epi_super_UMAP_clusters.pdf'), width = 7, height = 6)
DimPlot(super_set, reduction = 'umap', label = TRUE, group.by = 'seurat_clusters', cols = super_cols) + ggtitle('') &
  theme(axis.text = element_blank(), axis.title = element_text(size = 8), axis.ticks = element_blank())
dev.off()

pdf(file = paste0(results_dir, '/epi_super_markers_violin.pdf'), width = 4, height = 6)
print(scCustomize::Stacked_VlnPlot(super_set, superficial_marker_genes, pt.size = 0, group.by = 'seurat_clusters',
                                   colors = super_cols, x_lab_rotate = TRUE))
dev.off()

pdf(file = paste0(results_dir, '/epi_super_superficialScore_violin.pdf'), width = 5, height = 3)
print(VlnPlot(super_set, features = 'superficial1', group.by = 'seurat_clusters', cols = super_cols, pt.size = 0) +
        ggtitle('Superficial Module Score') + theme(legend.position = 'none'))
dev.off()


####  Marker + module feature plots  #############################################
DefaultAssay(super_set) <- 'RNA'
pdf(file = paste0(results_dir, '/epi_super_UMAP_markers.pdf'), width = 12, height = 7)
print(FeaturePlot(super_set, features = superficial_marker_genes, order = TRUE, ncol = 3) &
        scale_colour_gradientn(colours = brewer.pal(n = 9, 'YlGnBu')[3:9]) &
        theme(axis.text = element_blank(), axis.title = element_text(size = 9), axis.ticks = element_blank()))
dev.off()

pdf(file = paste0(results_dir, '/epi_super_modulePlot_superficial.pdf'), width = 5, height = 4)
print(module_plot(super_set, 'Superficial', superficial_markers))
dev.off()


####  Metallothioneins  ######################################################
## : metallothionein module score on the UMAP, Proximal, GERD&HC vs SSc.
metal_plots <- lapply(c('Proximal', 'Distal'), function(region) {
  lapply(list(c('GERD', 'HC'), 'SSc'), function(grp) {
    FeaturePlot(subset(super_set, location == region & condition %in% grp), features = 'metal_score1',
                order = TRUE, pt.size = 1) & labs(x = 'UMAP 1', y = 'UMAP 2') &
      scale_colour_gradientn(colours = c('grey85', 'grey75', brewer.pal(n = 9, 'OrRd')[4:9], 'black'),
                             name = 'Relative\nExpression\n',
                             breaks = c(min(super_set$metal_score1), max(super_set$metal_score1)),
                             labels = c('Min', 'Max'), guide = guide_colorbar(frame.colour = 'black', ticks = FALSE)) &
      ggtitle(paste(region, ',', paste(grp, collapse = ' & ')), subtitle = 'Metallothioneins') &
      theme(plot.title = element_text(face = 'bold', size = 10),
            plot.subtitle = element_text(face = 'italic', size = 8, hjust = 0.5),
            axis.text = element_blank(), axis.title = element_text(size = 8), axis.ticks = element_blank())
  })
})
pdf(file = paste0(results_dir, '/epi_super_metallothionein_UMAP.pdf'), width = 7, height = 7)
print(ggarrange(plotlist = c(metal_plots[[1]], metal_plots[[2]]), nrow = 2, ncol = 2, common.legend = TRUE, legend = 'right'))
dev.off()

## : metallothionein module score by subcluster and condition (Proximal).
pdf(file = paste0(results_dir, '/epi_super_metallothionein_violin_byCluster.pdf'), width = 7, height = 3)
print(VlnPlot(super_set.prox, features = 'metal_score1', group.by = 'seurat_clusters', split.by = 'condition',
              cols = cond_cols[c(3, 2, 1)], pt.size = 0))
dev.off()

## : metallothionein module score by location and condition (HC vs SSc).
metallothionein_loc_df <- super_set@meta.data[super_set$condition %in% c('HC', 'SSc'), c('location', 'condition', 'metal_score1')]
metallothionein_loc_df$condition <- factor(metallothionein_loc_df$condition, levels = c('HC', 'SSc'))
metallothionein_loc_df$location  <- factor(metallothionein_loc_df$location,  levels = c('Proximal', 'Distal'))
pdf(file = paste0(results_dir, '/epi_super_metallothionein_module_byLocation.pdf'), width = 5, height = 4)
print(ggplot(metallothionein_loc_df, aes(x = condition, y = metal_score1, fill = condition)) +
        geom_violin(scale = 'width', draw_quantiles = 0.5) + facet_wrap(~location) +
        scale_fill_manual(values = cond_cols[c(3, 1)]) + theme_bw() +
        labs(y = 'Metallothionein Module Score', subtitle = paste(metallo_genes, collapse = ', ')) +
        theme(legend.position = 'none', axis.title.x = element_blank(),
              plot.subtitle = element_text(face = 'italic', size = 7)))
dev.off()


####  Transcription factors and target genes  ##################################
tf_panels <- list(
  IRF1   = list(feats = c('IRF1',   'irf1_target_genes1'),   cols = c('grey85', 'grey75', brewer.pal(9, 'YlOrRd')[c(2:7, 9)])),
  MYC    = list(feats = c('MYC',    'myc_target_genes1'),    cols = c('grey85', 'grey75', brewer.pal(9, 'YlGnBu')[4:9], 'black')),
  E2F4   = list(feats = c('E2F4',   'e2f4_target_genes1'),   cols = c('grey85', 'grey75', brewer.pal(9, 'PuRd')[4:9],  'black')),
  NFE2L2 = list(feats = c('NFE2L2', 'nfe2l2_target_genes1'), cols = c('grey85', 'grey75', brewer.pal(9, 'YlGn')[2:7],  'black')))

for (tf in names(tf_panels)) {
  plots <- tf_featurePlot(super_set.prox, tf_panels[[tf]]$feats, tf_panels[[tf]]$cols)
  pdf(file = paste0(results_dir, '/epi_super_modulePlots_', tolower(tf), 'Genes.pdf'), width = 10, height = 6)
  print(patchwork::wrap_plots(plots, nrow = 2))
  dev.off()
}

## TF target-gene module scores by subcluster and condition (Proximal), via vioplot.
for (set in c('irf1_target_genes1', 'myc_target_genes1', 'e2f4_target_genes1', 'nfe2l2_target_genes1')) {
  pdf(file = paste0(results_dir, '/epiSuperClusters_', word(set, 1, sep = '\\_'), 'TargetGenes.pdf'), width = 4, height = 2)
  par(mfrow = c(1, 5), mar = c(2, 0.4, 2, 0.4))
  pd <- super_set.prox@meta.data
  pd$condition <- factor(pd$condition, c('HC', 'GERD', 'SSc'))
  for (l in 1:5) {
    rng  <- range(pd[, set]); ymin <- min(rng) - 0.1 * sum(abs(rng)); ymax <- max(rng) + 0.1 * sum(abs(rng))
    vioplot(as.formula(paste(set, '~condition')), data = pd[pd$seurat_clusters == l, ], ylab = '', xlab = '',
            col = rev(cond_cols), yaxt = 'n', xaxt = 'n', colMed = 'black', colMed2 = rev(cond_cols),
            rectCol = 'white', ylim = c(ymin, ymax),
            panel.first = { axis(2, tck = 1, col.ticks = 'gray80', lty = 3) }, pchMed = 23, lwd = 1.05)
    mtext(l, side = 3, padj = -1, las = 1, cex = 0.8, col = 'grey30')
  }
  dev.off()
}

## TF expression vs target-module correlations (reported in the text).
tf_cors <- sapply(c('IRF1', 'MYC', 'E2F4', 'NFE2L2'), function(tf) {
  cor.test(super_set@assays$RNA$data[tf, ], super_set[[paste0(tolower(tf), '_target_genes1')]][, 1])$estimate
})
write.table(data.frame(TF = names(tf_cors), r = round(tf_cors, 3)),
            paste0(results_dir, '/epi_super_TF_target_correlations.tsv'),
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')


####  FLI1 target-gene module  ###################################################
## FLI1 (an endothelial transcription factor) target genes that are differentially expressed in the
## superficial pseudobulk comparison. Operates on the superficial cells, so it lives here for now;
## it depends on 06_pseudobulk_dgea.R outputs (the ENCODE/ChEA TF map and the by-location pseudobulk DEGs),
## so this script must run after 06_pseudobulk_dgea.R. (The reorganization pass may relocate this, e.g.
## alongside the rest of the endothelial  content.)

annot_tf <- read.table(paste0(pathway_maps_dir, '/ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.map'),
                       header = FALSE, sep = '\t', quote = '')
load(paste0(results_dir, '/pb_deg_res_layerRep_byRegion.RData'))   # brings in ct.genes + qlf.res

## FLI1 ENCODE targets that are in the superficial background and significant in the superficial
## by-location pseudobulk DEG result. qlf.res[[30]] is the last of the 30 by-location entries =
## the superficial Proximal SSc-vs-GERD comparison (FLAG: confirm this is the intended contrast).
fli1_target_genes <- annot_tf[annot_tf$V1 == 'FLI1 ENCODE', 'V2']
fli1_target_genes <- fli1_target_genes[fli1_target_genes %in% ct.genes[['Superficial']]]
fli1_target_genes <- fli1_target_genes[fli1_target_genes %in% rownames(qlf.res[[30]][qlf.res[[30]]$p_adj < 0.05, ])]
super_set <- AddModuleScore(super_set, features = list(fli1_target_genes), name = 'fli1_target_genes', assay = 'RNA')

write.table(super_set@meta.data[, c('fli1_target_genes1', 'condition'), drop = FALSE],
            paste0(suppdata_dir, '/fli1_target_module.tsv'), quote = FALSE, row.names = TRUE, col.names = TRUE, sep = '\t')

pdf(file = paste0(results_dir, '/epi_super_fli1_violin.pdf'), width = 5, height = 6)
vioplot(fli1_target_genes1 ~ condition, data = super_set@meta.data[super_set$location == 'Distal', ],
        col = rev(cond_cols), ylab = 'FLI1 Target-Gene Module Score', xlab = '')
dev.off()


####  SSc-associated immune genes (subcluster 5 vs 1-4)  #########################
immune_plots <- lapply(immune_genes, function(g) gene_cluster5vioPlot(super_set.prox, g))
pdf(file = paste0(results_dir, '/epi_super_immuneGenes_cluster5.pdf'), width = 12, height = 8)
print(cowplot::plot_grid(plotlist = immune_plots, ncol = 4))
dev.off()


####  Subcluster proportions by condition (MASC)  ######################################
## Per-sample subcluster proportions feed the  bars; MASC gives the significance,
## computed within each biopsy location.
test_df <- super_set@meta.data[, c('orig.ident', 'location', 'condition', 'seurat_clusters')]
test_df$ssc_hc   <- ifelse(test_df$condition == 'SSc',  1, ifelse(test_df$condition == 'HC',   0, NA))
test_df$gerd_hc  <- ifelse(test_df$condition == 'GERD', 1, ifelse(test_df$condition == 'HC',   0, NA))
test_df$ssc_gerd <- ifelse(test_df$condition == 'SSc',  1, ifelse(test_df$condition == 'GERD', 0, NA))

prop_comps <- data.frame()
for (comp in c('ssc_hc', 'gerd_hc', 'ssc_gerd')) {
  for (l in c('Proximal', 'Distal')) {
    temp_df <- test_df[test_df$location == l, ]
    masc <- MASC(temp_df, cluster = temp_df$seurat_clusters, contrast = comp, random_effects = 'orig.ident')
    masc$cluster  <- gsub('cluster', '', masc$cluster)
    masc$comp     <- comp
    masc$location <- l
    prop_comps <- rbind(prop_comps, setNames(masc, names(masc)))
  }
}
## FDR-correct within each location stratum.
for (l in c('Proximal', 'Distal')) {
  prop_comps[prop_comps$location == l, 'masc.p.fdr'] <- p.adjust(prop_comps[prop_comps$location == l, 'model.pvalue'], method = 'fdr')
}
write.table(prop_comps, paste0(results_dir, '/epiSuper_clusters_props_stats.tsv'),
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')

## Per-sample superficial-subcluster proportions (consumed by 11_clinical_associations.R).
super_props <- as.data.frame(prop.table(table(super_set$orig.ident, super_set$seurat_clusters), margin = 1))
names(super_props) <- c('Sample', 'Cluster', 'Proportion')
super_props$Cluster  <- as.integer(as.character(super_props$Cluster))
super_props$Location <- ifelse(grepl('-P$', super_props$Sample), 'Proximal', 'Distal')
write.table(super_props, paste0(results_dir, '/epiSuper_cluster_props.tsv'),
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')


####  FIGURE SOURCE DATA  #########################################################################

write.table(FetchData(super_set, c('UMAP_1', 'UMAP_2', 'seurat_clusters')),
            paste0(suppdata_dir, '/superficial_umap_coords.tsv'), quote = FALSE, row.names = TRUE, col.names = TRUE, sep = '\t')
write.table(FetchData(super_set, c('seurat_clusters', superficial_marker_genes)),
            paste0(suppdata_dir, '/superficial_marker_expression.tsv'), quote = FALSE, row.names = TRUE, col.names = TRUE, sep = '\t')
write.table(super_set@meta.data[, c('seurat_clusters', 'superficial1'), drop = FALSE],
            paste0(suppdata_dir, '/superficial_module_score.tsv'), quote = FALSE, row.names = TRUE, col.names = TRUE, sep = '\t')
write.table(FetchData(super_set.prox, c('UMAP_1', 'UMAP_2', 'condition', 'metal_score1')),
            paste0(suppdata_dir, '/metallothionein_umap.tsv'), quote = FALSE, row.names = TRUE, col.names = TRUE, sep = '\t')
write.table(super_set@meta.data[, c('seurat_clusters', 'condition', 'metal_score1'), drop = FALSE],
            paste0(suppdata_dir, '/metallothionein_by_subcluster.tsv'), quote = FALSE, row.names = TRUE, col.names = TRUE, sep = '\t')
write.table(super_set@meta.data[, paste0('integrated_snn_res.', RES_SCAN), drop = FALSE],
            paste0(suppdata_dir, '/superficial_resolution_scan.tsv'), quote = FALSE, row.names = TRUE, col.names = TRUE, sep = '\t')
write.table(FetchData(super_set, c('UMAP_1', 'UMAP_2', superficial_marker_genes, 'superficial1')),
            paste0(suppdata_dir, '/superficial_marker_module_umap.tsv'), quote = FALSE, row.names = TRUE, col.names = TRUE, sep = '\t')
for (tf in c('irf1', 'myc', 'e2f4', 'nfe2l2')) {
  write.table(FetchData(super_set, c('UMAP_1', 'UMAP_2', 'condition', 'seurat_clusters',
                                     toupper(tf), paste0(tf, '_target_genes1'))),
              paste0(suppdata_dir, '/tf_', tf, '_targets.tsv'), quote = FALSE, row.names = TRUE, col.names = TRUE, sep = '\t')
}
write.table(FetchData(super_set, c('condition', 'seurat_clusters', 'CD44', 'CD74', 'C10orf99', 'HBEGF', 'HLA-B')),
            paste0(suppdata_dir, '/cluster5_immune_genes.tsv'), quote = FALSE, row.names = TRUE, col.names = TRUE, sep = '\t')
