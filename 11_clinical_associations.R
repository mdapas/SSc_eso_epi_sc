####################################################################################################
#  11_clinical_associations.R
#
#  Association of esophageal epithelial phenotypes with clinical motility/sensory measures.
#
#  Paper:
#   Dapas et al. (2025) "Cellular and molecular dysregulation of the esophageal epithelium in
#   systemic sclerosis." JCI Insight.
#
#  Purpose / figures produced:
#   - Correlation heatmap of the quantitative clinical traits  (source -> clinical_trait_correlations.tsv)
#   - PCA biplot of the quantitative clinical traits, colored by HRM classification
#   - Mean superficial metallothionein module score vs clinical-trait PC1 / PC2 (Proximal)
#   - Superficial-EEC proportion vs PC1 (and the superficial subcluster proportions vs PC1)
#   - Superficial-EEC proportion by motility phenotype (FLIP / HRM / NM)
#   - compartment proportions by HRM / FLIP phenotype
#   - superficial proportion vs SSc disease duration
#   - EEC proportions vs PC1 / PC2 (all compartments)  (source -> compartment_props_clinical_pcs.tsv)
#   - TRIM11 pseudobulk expression by motility phenotype
#   Statistics: ordinal logistic regression (rms::orm) of phenotype on each measure, combined across
#   locations with ACAT; linear regression on the clinical-trait PCs.
#
#  Pipeline position:
#   06_pseudobulk_dgea.R (pseudobulk obj + DEGs) + 05_epithelial_compartments.R (compartment proportions) +
#   08_superficial_analysis.R (superficial subclusters + module scores)  ->  11_clinical_associations.R
#
#  Inputs:
#   - <pb_dir>/pb_deg_res_layerRep_byRegion.RData   (06_pseudobulk_dgea.R) ct.genes + qlf.res
#   - <pb_dir>/pbObj_layerRep_byRegion.RData         (06_pseudobulk_dgea.R) the pseudobulk DGEList (pb_obj)
#   - <results_dir>/eso_epi_superficial.RData        (epiSuper) super_set with module scores
#   - <results_dir>/epi_layerRep_props.tsv           (epi_analysis) per-sample compartment proportions
#   - <results_dir>/epiSuper_cluster_props.tsv       (epiSuper) per-sample superficial-cluster props
#   - <clin_dir>/clinical_data.txt                   per-participant clinical traits + motility class
#   - <clin_dir>/ssc_clin_dat.txt                    SSc disease duration + superficial proportions
#   - <results_dir>/diffScoreStats.tsv               per-sample differentiation-score stats
#
#  Outputs (to <results_dir>): the .pdf panels above, clinical_trait_correlations.tsv, compartment_props_clinical_pcs.tsv, pca_biplot_data.csv
#
#  Usage:
#   Rscript 11_clinical_associations.R   (run from the repository root)
#
#  Notes:
#   - Achalasia samples (SSc4, SSc2) are excluded from the ordinal
#     motility-phenotype regressions.
#   - Sample/patient IDs are the published de-identified IDs (see config/sample_map.txt); the
#     clinical_data.txt / ssc_clin_dat.txt tables are expected to key on those de-identified IDs.
#   - FLAG:  PCs use FactoMineR::PCA on the raw trait matrix (object `pca`). The Methods
#     describe mean-imputation of the few missing values; an imputed PCA (`pca_imp`) is also computed
#     for reference. Confirm which was used for the published figure (they are very similar given
#     only ~6/90 values are missing).
#   - Removed from the original working file: an SSc cutaneous-subtype PCA (rebuttal-only), and
#     exploratory gene-vs-PC / apoptosis / BAK:BCL2-ratio analyses not shown in the paper.
####################################################################################################


####  LOAD LIBRARIES  ##############################################################################

suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(edgeR)
  library(rms)            # orm() ordinal regression
  library(ACAT)           # ACAT() combined p-values
  library(FactoMineR)     # PCA()
  library(factoextra)     # fviz_pca_biplot()
  library(UpSetR)         # upset()
  library(beeswarm)
  library(scales)
  library(stringr)
  library(pals)           # alphabet()
})


####  SET PARAMETERS  ##############################################################################

set.seed(0123456789)

## ---- Paths (relative to repository root) ----
results_dir <- './results'
pb_dir      <- './results'                 # pseudobulk objects (06_pseudobulk_dgea.R outputs)
clin_dir    <- './results/clinical_data'   # clinical data tables

## ---- Colors ----
cond_cols    <- c('#8156B3', '#e26b53', '#A8A39d')                          # SSc, GERD, HC
layer_cols_5 <- c('#4ca5b1', '#7fbf7b', '#e9f5a3', '#e9b85d', '#c72f4c')    # 5 compartments

## ---- Contrasts / compartments ----
comps      <- list(c('SSc', 'HC'), c('GERD', 'HC'), c('SSc', 'GERD'))
layers_rep <- factor(c('Basal', 'ProliferatingBasal', 'ProliferatingSuprabasal', 'Suprabasal', 'Superficial'),
                     levels = c('Basal', 'ProliferatingBasal', 'ProliferatingSuprabasal', 'Suprabasal', 'Superficial'))

## ---- Clinical traits ----
## clin_vars: all quantitative traits; pca_vars: the 9 traits entered into the PCA.
clin_vars <- c('median_irp', 'mean_DCI', 'pressure_60ml', 'egd_di_60ml', 's_basal_ee_egjp', 'egj_cont_index',
               'max_acid_exposure', 'gerdq_impact_score', 'bedq_score', 'neqol_score')
pca_vars  <- c('median_irp', 'mean_DCI', 'pressure_60ml', 'egd_di_60ml', 's_basal_ee_egjp', 'egj_cont_index',
               'gerdq_impact_score', 'bedq_score', 'neqol_score')

## Ordered motility-phenotype levels.
hrm_vals  <- c('Normal', 'Ineffective', 'Absent')
flip_vals <- c('Normal', 'Diminished/Impaired', 'Absent')
nm_vals   <- c('Normal', 'Ineffective', 'Absent')


####  FUNCTIONS  ###################################################################################

loadRData <- function(fileName) {
  load(fileName)
  get(ls()[ls() != 'fileName'])
}


####  LOAD DATA  ##################################################################################

load(paste0(pb_dir, '/pb_deg_res_layerRep_byRegion.RData'))   # ct.genes + qlf.res
load(paste0(pb_dir, '/pbObj_layerRep_byRegion.RData'))        # pb_obj (pseudobulk DGEList)
super_set <- loadRData(paste0(results_dir, '/eso_epi_superficial.RData'))

clin_dat <- read.delim(paste0(clin_dir, '/clinical_data.txt'), na.strings = c('NA', 'NULL', 'null', ''))

## Attach clinical traits + motility classifications to the pseudobulk sample metadata.
pb_obj$samples$indv <- sapply(pb_obj$samples$sample, function(x) strsplit(x, '-')[[1]][1])
pb_obj$samples <- merge(pb_obj$samples,
                        clin_dat[, c('Participant_ID', clin_vars, 'cc_v4_class_hrm', 'cc_v4_class_peristalsis',
                                     'flip_contract_pattern_v2', 'NM_model_motility')],
                        by.x = 'indv', by.y = 'Participant_ID')
pb_obj$samples$cc_v4_class_hrm <- toupper(pb_obj$samples$cc_v4_class_hrm)

## Recode motility phenotypes onto the ordered 3-level scales.
pb_obj$samples$cc_v4_class_peristalsis[pb_obj$samples$cc_v4_class_peristalsis == 'premature'] <- NA
pb_obj$samples$cc_v4_class_peristalsis <- mapvalues(pb_obj$samples$cc_v4_class_peristalsis,
                                                    c('normal', 'ineffective', 'absent'), hrm_vals)
pb_obj$samples$flip_contract_pattern_v2 <- sub(' CONTRACTILE RESPONSE', '', pb_obj$samples$flip_contract_pattern_v2)
pb_obj$samples$flip_contract_pattern_v2 <- mapvalues(pb_obj$samples$flip_contract_pattern_v2,
                                                     c('NORMAL', 'BORDERLINE/DIMINISHED', 'IMPAIRED/DISORDERED', 'ABSENT'),
                                                     flip_vals[c(1, 2, 2, 3)])
pb_obj$samples$NM_model_motility <- mapvalues(pb_obj$samples$NM_model_motility,
                                              c('Normal', 'Stage I  - Ineffective', 'Stage II - Ineffective', 'Stage III - Absent'),
                                              nm_vals[c(1, 2, 2, 3)])

## Differentiation-score stats and per-sample compartment proportions.
diffScore_dat <- read.delim(paste0(results_dir, '/diffScoreStats.tsv'))
pb_obj$samples$sample <- sapply(pb_obj$samples$sample, function(x) strsplit(x, '_')[[1]][1])
pb_obj$samples <- merge(pb_obj$samples, diffScore_dat, by.x = 'sample', by.y = 'orig.ident')

prop_dat <- read.delim(paste0(results_dir, '/epi_layerRep_props.tsv'))
prop_dat$LayerRep <- mapvalues(prop_dat$LayerRep,
                               c('Basal', 'Proliferating basal', 'Proliferating suprabasal', 'Suprabasal', 'Superficial'),
                               levels(layers_rep))
pb_obj$samples <- merge(pb_obj$samples, prop_dat[, c('Sample', 'LayerRep', 'Proportion')],
                        by.x = c('sample', 'layer'), by.y = c('Sample', 'LayerRep'))

## Mean superficial module scores per sample (require the epiSuper module scores).
sup_means <- super_set@meta.data
add_mean <- function(col) aggregate(as.formula(paste0(col, '~orig.ident')), sup_means, FUN = mean)
for (m in c('metal_score1', 'irf1_target_genes1', 'myc_target_genes1', 'e2f4_target_genes1', 'nfe2l2_target_genes1')) {
  agg <- add_mean(m)
  pb_obj$samples[[sub('1$', '', sub('_score', '', sub('_genes', '', m)))]] <-
    agg[match(pb_obj$samples$sample, agg$orig.ident), 2]
}
names(pb_obj$samples)[names(pb_obj$samples) == 'metal'] <- 'metal_mean'

## Exclude achalasia for the ordinal motility regressions.
pb_obj_noAc <- pb_obj[, !pb_obj$samples$indv %in% c('SSc4', 'SSc2')]   # exclude achalasia samples

layer    <- factor(pb_obj$samples$layer, levels = levels(layers_rep))
condition <- factor(pb_obj$samples$condition, levels = c('SSc', 'GERD', 'HC'))
location  <- factor(pb_obj$samples$location, levels = c('Distal', 'Proximal'))
layer_noAc    <- factor(pb_obj_noAc$samples$layer, levels = levels(layers_rep))
condition_noAc <- factor(pb_obj_noAc$samples$condition, levels = c('SSc', 'GERD', 'HC'))
location_noAc  <- factor(pb_obj_noAc$samples$location, levels = c('Distal', 'Proximal'))

pheno_flip <- factor(pb_obj_noAc$samples$flip_contract_pattern_v2, ordered = TRUE, levels = flip_vals)
pheno_hrm  <- factor(pb_obj_noAc$samples$cc_v4_class_peristalsis,  ordered = TRUE, levels = hrm_vals)
pheno_nm   <- factor(pb_obj_noAc$samples$NM_model_motility,        ordered = TRUE, levels = nm_vals)
phenos     <- list(flip = pheno_flip, hrm = pheno_hrm, nm = pheno_nm)


####   10A-B  -  compartment proportions x motility phenotype  ######################

pdf(file = paste0(results_dir, '/superProp_motilPhenos.pdf'), width = 9, height = 6)
par(mfrow = c(2, 3))
temp.dat <- cbind.data.frame(prop = pb_obj_noAc$samples$Proportion, layer_noAc, location_noAc, condition_noAc,
                             pheno_flip, pheno_hrm, pheno_nm)
for (loc in c('Proximal', 'Distal')) {
  for (p in c('pheno_flip', 'pheno_hrm', 'pheno_nm')) {
    dat <- temp.dat[temp.dat$layer == 'Superficial' & temp.dat$location == loc, ]
    boxplot(as.formula(paste0('prop~', p)), dat, main = paste(loc, sub('pheno_', '', p)))
    beeswarm(as.formula(paste0('prop~', p)), dat, pch = 21, add = TRUE, cex = 1.4, pwbg = cond_cols[dat$condition])
  }
}
dev.off()

## Ordinal regression of phenotype on compartment proportion (SSc samples), per compartment x
## location, combined across locations by ACAT. (Proliferating compartments expressed relative to
## their parent compartment.)
ord.dat <- temp.dat[temp.dat$condition == 'SSc', ]
ord.dat[ord.dat$layer == 'ProliferatingBasal', 'prop']      <- ord.dat[ord.dat$layer == 'ProliferatingBasal', 'prop'] / ord.dat[ord.dat$layer == 'Basal', 'prop']
ord.dat[ord.dat$layer == 'ProliferatingSuprabasal', 'prop'] <- ord.dat[ord.dat$layer == 'ProliferatingSuprabasal', 'prop'] / ord.dat[ord.dat$layer == 'Suprabasal', 'prop']

ord_reg.prop.res <- lapply(setNames(levels(layers_rep), levels(layers_rep)), function(l) {
  as.data.frame(do.call(rbind, lapply(seq_along(phenos), function(p) {
    pheno_name <- paste0('pheno_', names(phenos)[p])
    res <- c()
    for (loc in c('Proximal', 'Distal')) {
      suffix <- paste(names(phenos)[p], substr(loc, 1, 1), sep = '.')
      stats <- tryCatch({
        prop.test <- orm(as.formula(paste0(pheno_name, '~prop')), ord.dat[ord.dat$layer == l & ord.dat$location == loc, ])
        b_test <- coef(prop.test)['prop']
        p_test <- pchisq(b_test^2 / vcov(prop.test)[2, 2], 1, lower.tail = FALSE)
        c(b_test, p_test, sign(b_test) * qnorm(p_test / 2))
      }, error = function(e) rep(NA, 3))
      names(stats) <- paste0(c('b_', 'p_', 'z_'), suffix)
      res <- c(res, stats)
    }
    res <- t(as.data.frame(res)); rownames(res) <- pheno_name; res
  })))
})


####  TRIM11 pseudobulk expression x motility phenotype  ############################
## Transcription factors of interest (from the TF enrichment) plus TRIM11.

tf_genes   <- c('IRF1', 'MYC', 'E2F4', 'NFE2L2', 'TRIM11')
count_norm <- as.data.frame(log(cpm(pb_obj_noAc[tf_genes, ]) + 1))

pdf(file = paste0(results_dir, '/PseudoDEG_super_motilPhenos_TRIM11.pdf'), width = 7, height = 4)
par(mfrow = c(1, 3), mar = c(3, 2, 4.5, 1), oma = c(0.5, 2.5, 0, 0))
for (test in c('FLIP', 'HRM', 'NM')) {
  pheno <- switch(test, FLIP = pheno_flip, HRM = pheno_hrm, NM = pheno_nm)
  data <- cbind.data.frame(gene = as.numeric(t(count_norm['TRIM11', ])), layer_noAc, pheno, condition_noAc, location_noAc)
  data <- data[data$layer_noAc == 'Superficial' & data$location == 'Proximal', ]
  boxplot(gene ~ pheno, data = data, ylab = '', yaxt = 'n', xlab = '', pch = 16)
  beeswarm(gene ~ pheno, data = data, add = TRUE, pwbg = cond_cols[data$condition_noAc], pch = 21, cex = 1.5)
  axis(2, at = 0:3, labels = 0:3)
  title(main = test, line = 3, cex.main = 1.5, font.main = 2)
  if (test == 'FLIP') mtext('log(CPM+1)', 2, line = 2.8)
  mtext(bquote(bolditalic('TRIM11')), side = 3, padj = -1.5, las = 1, cex = 0.8, col = 'grey30')
}
dev.off()

## Ordinal regression of phenotype on each TF gene (Proximal superficial), combined by ACAT.
gene_ord <- do.call(rbind, lapply(tf_genes, function(g) {
  res <- c(gene = g)
  for (p in seq_along(phenos)) {
    suffix <- paste(names(phenos)[p], 'P', sep = '.')
    gene.dat <- cbind.data.frame(gene = as.numeric(t(count_norm[g, ])), layer_noAc, location_noAc, pheno = phenos[[p]])
    stats <- tryCatch({
      gt <- orm(pheno ~ gene, gene.dat[gene.dat$layer == 'Superficial' & gene.dat$location == 'Proximal', ])
      b <- coef(gt)['gene']; pv <- pchisq(b^2 / vcov(gt)[2, 2], 1, lower.tail = FALSE)
      c(b, pv)
    }, error = function(e) rep(NA, 2))
    names(stats) <- paste0(c('b_', 'p_'), suffix)
    res <- c(res, stats)
  }
  data.frame(t(res), stringsAsFactors = FALSE)
}))
write.table(gene_ord, paste0(results_dir, '/clinical_TFgene_ordinalReg.tsv'), quote = FALSE, row.names = FALSE, sep = '\t')


####  Clinical-trait PCA  ##############################################################
## Missing-data overview ( imputation panels).
upset_df <- as.data.frame(!is.na(clin_dat[clin_dat$cond == 'SSc', pca_vars])) * 1
names(upset_df) <- c('Median IRP', 'Mean DCI', 'Pressure 60ml', 'EGJ distensibility\nindex', 'End-expiratory\nEGJ Pressure',
                     'EGJ contractile\nindex', 'GerdQ', 'BEDQ', 'NEQOL')
pdf(file = paste0(results_dir, '/missingDat_imputation_upset.pdf'), width = 4, height = 4.5)
upset(upset_df, nsets = 10, mb.ratio = c(0.3, 0.7), point.size = 2.5, text.scale = 1.2, set_size.show = TRUE,
      sets.x.label = 'N Values', mainbar.y.label = 'N Patients', set_size.scale_max = 13)
dev.off()

## SSc, Distal, Superficial rows define one row per individual for the trait PCA.
ssc_rows <- pb_obj$samples$layer == 'Superficial' & pb_obj$samples$condition == 'SSc' & pb_obj$samples$location == 'Distal'
pca     <- PCA(pb_obj$samples[ssc_rows, pca_vars], graph = FALSE)
pca_imp <- PCA(pb_obj$samples[ssc_rows, pca_vars] %>%                      # imputed reference (see header FLAG)
                 mutate(across(everything(), ~ifelse(is.na(.), mean(., na.rm = TRUE), .))), graph = FALSE)

pdf(file = paste0(results_dir, '/SSc_quantClinVar_PCA_hrm.pdf'), width = 8, height = 5.5)
print(fviz_pca_biplot(pca, col.ind = pb_obj$samples[ssc_rows, 'cc_v4_class_hrm'],
                      addEllipses = FALSE, label = 'var', col.var = 'black', repel = TRUE,
                      legend.title = 'HRM', pointsize = 4, invisible = c('quali'), xlim = c(-3.5, 3.5), ylim = c(-3.2, 3.2),
                      xlab = paste0('PC1 (', round(pca$eig[1, 2], 2), '%)'),
                      ylab = paste0('PC2 (', round(pca$eig[2, 2], 2), '%)'), title = '') +
        theme_minimal() + theme(panel.background = element_rect(fill = rgb(254, 252, 250, maxColorValue = 255))))
dev.off()

##  source data (biplot coordinates).
ind_coords <- as.data.frame(pca$ind$coord[, 1:2]); ind_coords$name <- pb_obj$samples[ssc_rows, 'indv']
ind_coords$hrm_class <- pb_obj$samples[ssc_rows, 'cc_v4_class_hrm']
var_coords <- as.data.frame(pca$var$coord[, 1:2]); var_coords$name <- rownames(var_coords); var_coords$hrm_class <- NA
write.csv(rbind(ind_coords, var_coords), paste0(results_dir, '/pca_biplot_data.csv'), row.names = FALSE)

## Attach the PCs back to the sample metadata for the regressions below.
pcs <- as.data.frame(pca$ind$coord[, 1:2]); colnames(pcs) <- c('PC1', 'PC2')
pcs$indv <- pb_obj$samples[ssc_rows, 'indv']
pb_obj$samples <- merge(pb_obj$samples, pcs, by = 'indv', all.x = TRUE)
write.table(pb_obj$samples, paste0(results_dir, '/compartment_props_clinical_pcs.tsv'), quote = FALSE, row.names = FALSE, sep = '\t')


####  Clinical-trait correlation heatmap  ##############################################
trait_cor <- cor(pb_obj$samples[, pca_vars], use = 'pairwise.complete.obs', method = 'spearman')
write.table(trait_cor, paste0(results_dir, '/clinical_trait_correlations.tsv'), quote = FALSE, row.names = TRUE, col.names = TRUE, sep = '\t')
pdf(file = paste0(results_dir, '/clinTrait_correlation_heatmap.pdf'), width = 7, height = 6)
pheatmap::pheatmap(trait_cor, display_numbers = TRUE, breaks = seq(-1, 1, length.out = 101),
                   color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'RdBu')))(100))
dev.off()


####   10D  -  EEC proportion vs clinical PCs  ######################################

names(layer_cols_5) <- levels(layers_rep)
pdf(file = paste0(results_dir, '/clinQuantTraitPCs_superProp.pdf'), width = 8, height = 10)
par(mfrow = c(5, 4), mar = c(2, 3, 2, 1), mgp = c(0.5, 1, 0))
for (l in levels(layers_rep)) {
  for (quant_var in c('PC1', 'PC2')) {
    for (loc in c('Proximal', 'Distal')) {
      temp.dat <- pb_obj$samples[pb_obj$samples$layer == l & pb_obj$samples$location == loc, ]
      model <- lm(as.formula(paste('Proportion ~', quant_var)), data = temp.dat)
      p <- summary(model)$coefficients[2, 4]
      mod_seq <- seq(min(temp.dat[, quant_var], na.rm = TRUE), max(temp.dat[, quant_var], na.rm = TRUE), length.out = 200)
      new.dat <- setNames(data.frame(mod_seq), quant_var)
      preds <- predict(model, newdata = new.dat, interval = 'confidence')
      plot(as.formula(paste('Proportion ~', quant_var)), data = temp.dat, ylab = '', xlab = quant_var,
           pch = 16, main = loc, cex = 0.1, xaxt = 'n')
      polygon(c(rev(mod_seq), mod_seq), c(rev(preds[, 3]), preds[, 2]), col = alpha(layer_cols_5[l], 0.2), border = NA)
      abline(model$coefficients[1], model$coefficients[2], lty = 3, col = 'pink3', lwd = 2)
      points(as.formula(paste('Proportion ~', quant_var)), data = temp.dat, pch = 21, lwd = 1.2, bg = layer_cols_5[l])
      text(max(temp.dat[, quant_var], na.rm = TRUE), max(temp.dat$Proportion, na.rm = TRUE) * 0.99,
           paste('P =', round(p, 3)), adj = 1, col = 'darkred', cex = 0.75)
    }
  }
}
dev.off()


####  Superficial metallothionein vs clinical PCs  #####################################
plot_dat <- pb_obj$samples[pb_obj$samples$layer == 'Superficial' & pb_obj$samples$location == 'Proximal' &
                             pb_obj$samples$condition == 'SSc', ]
pdf(file = paste0(results_dir, '/metal_clinVarPCs.pdf'), width = 8, height = 4)
xrange <- c(-3.5, 3.5)
par(mfrow = c(1, 2), mar = c(4, 4, 1, 1), mgp = c(2.2, 0.8, 0))
for (pc in c('PC1', 'PC2')) {
  model <- lm(as.formula(paste0('metal_mean ~ ', pc)), data = plot_dat)
  p <- summary(model)$coefficients[2, 4]
  mod_seq <- seq(xrange[1] - 0.5, xrange[2] + 0.5, length.out = 200)
  preds <- predict(model, newdata = setNames(data.frame(mod_seq), pc), interval = 'confidence')
  plot(as.formula(paste0('metal_mean ~ ', pc)), data = plot_dat, pch = 16, cex = 0.1,
       xlab = pc, ylab = 'Mean metallothionein module score', xlim = xrange)
  polygon(c(rev(mod_seq), mod_seq), c(rev(preds[, 3]), preds[, 2]), col = alpha('mediumpurple', 0.2), border = NA)
  abline(model$coefficients[1], model$coefficients[2], lty = 3, col = 'purple', lwd = 2)
  points(as.formula(paste0('metal_mean ~ ', pc)), data = plot_dat, pch = 21, bg = cond_cols[1])
  text(xrange[2], max(plot_dat$metal_mean, na.rm = TRUE) * 0.99, paste('P =', round(p, 3)), adj = 1, col = 'darkred', cex = 0.8)
}
dev.off()


####  Superficial proportion vs SSc disease duration  ###############################

ssc_clin_dat <- read.delim(paste0(clin_dir, '/ssc_clin_dat.txt'))
pdf(file = paste0(results_dir, '/diseaseDur_superProp.pdf'), width = 10, height = 5)
par(mfrow = c(1, 2))
for (var in c('superProp_P', 'superProp_D')) {
  loc <- ifelse(var == 'superProp_P', 'Proximal', 'Distal')
  model <- lm(as.formula(paste(var, '~ Disease.Duration..mo.')), data = ssc_clin_dat)
  p <- summary(model)$coefficients[2, 4]
  mod_seq <- seq(0, 600, length.out = 200)
  preds <- predict(model, newdata = data.frame(Disease.Duration..mo. = mod_seq), interval = 'confidence')
  plot(as.formula(paste(var, '~ Disease.Duration..mo.')), data = ssc_clin_dat, ylab = 'Proportion Superficial',
       xlab = 'SSc Disease Duration (mo)', pch = 16, col = 'grey20', main = loc, cex = 1.4)
  polygon(c(rev(mod_seq), mod_seq), c(rev(preds[, 3]), preds[, 2]), col = alpha('lightpink', 0.2), border = NA)
  abline(model$coefficients[1], model$coefficients[2], lty = 3, col = 'pink3', lwd = 2)
  points(as.formula(paste(var, '~ Disease.Duration..mo.')), data = ssc_clin_dat, pch = 18, col = cond_cols[3])
  text(500, max(ssc_clin_dat[, var], na.rm = TRUE) * 0.99, paste('P =', round(p, 3)), adj = 1, col = 'darkred', cex = 0.7)
}
dev.off()


####  Superficial subcluster proportions vs PC1  #######################################
## Per-sample superficial-subcluster proportions associated with the clinical-trait PCs.

superClust_dat <- read.delim(paste0(results_dir, '/epiSuper_cluster_props.tsv'))
superClust_dat <- merge(superClust_dat, pb_obj$samples[, c('sample', 'location', 'PC1', 'PC2')],
                        by.x = c('Sample', 'Location'), by.y = c('sample', 'location'))

super_cols <- c('#1F9E89', '#E76254', '#67AD5B', '#E3B23C', '#9C8CC4')
pdf(file = paste0(results_dir, '/clinQuantTraitPCs_superClusterProp.pdf'), width = 7, height = 3.5)
par(mfrow = c(1, 2), mar = c(2, 3, 2, 1), mgp = c(0.5, 1, 0))
for (loc in c('Proximal', 'Distal')) {
  n <- 1
  for (cl in 1:5) {
    temp.dat <- superClust_dat[superClust_dat$Cluster == cl & superClust_dat$Location == loc, ]
    model <- lm(Proportion ~ PC1, data = temp.dat)
    p <- summary(model)$coefficients[2, 4]
    mod_seq <- seq(-3.5, 3.6, length.out = 200)
    preds <- predict(model, newdata = data.frame(PC1 = mod_seq), interval = 'confidence')
    if (n == 1) {
      plot(Proportion ~ PC1, data = temp.dat, main = loc, pch = 16, cex = 0.1, xlim = c(-3.5, 3.6), ylim = c(0, 0.55), xaxt = 'n')
    } else {
      plot(Proportion ~ PC1, data = temp.dat, main = '', pch = 16, cex = 0.1, xlim = c(-3.5, 3.6), ylim = c(0, 0.55), xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
    }
    lw <- ifelse(cl == 3, 0.6, 0.2)
    polygon(c(rev(mod_seq), mod_seq), c(rev(preds[, 3]), preds[, 2]), col = alpha(super_cols[cl], lw), border = NA)
    abline(model$coefficients[1], model$coefficients[2], lty = 3, col = super_cols[cl], lwd = 2)
    points(Proportion ~ PC1, data = temp.dat, pch = 21, bg = super_cols[cl], lwd = ifelse(cl == 3, 1.4, 0.8))
    if (cl == 3) text(3.5, 0.54, paste('SC3 P =', signif(p, 2)), adj = 1, col = 'darkgreen', cex = 0.75)
    par(new = TRUE); n <- n + 1
  }
}
dev.off()
