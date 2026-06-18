####################################################################################################
#  06_pseudobulk_dgea.R
#
#  Pseudobulk differential gene-expression analysis of the esophageal epithelial compartments.
#
#  Paper:
#   Dapas et al. (2025) "Cellular and molecular dysregulation of the esophageal epithelium in
#   systemic sclerosis." JCI Insight.
#
#  Purpose / figures produced:
#   - N differentially-expressed genes per compartment (Proximal/Distal),
#                   categorized SSc-only / GERD-only / shared. The category counts are computed
#                   here (comp_cols classification + the DEG-count summary); the stacked bars
#                   themselves were assembled in Illustrator from these counts.
#   - Proximal-vs-Distal expression correlation (variable genes) per compartment
#                   and condition (boxplots); source data -> proximal_distal_correlation.tsv.
#   - GERD-vs-HC vs SSc-vs-HC log2FC scatter for the Superficial compartment
#                   (Proximal & Distal), with disease-gene labels.
#   - Per-compartment log2FC jitter / percent-DE strip plots; source -> compartment_logfc_jitter.tsv.
#   - Per-comparison pseudobulk DEG tables (Excel).
#   - Transcription-factor (ENCODE/ChEA) over-representation among the
#                        Superficial-compartment DEGs (IRF1, MYC, E2F4, NFE2L2, ...), plus GO/KEGG
#                        over-representation tables.
#   ( GSEA dotplot, is produced by 07_gsea_enrichment.R, not here.)
#
#  Pipeline position:
#   04_reintegrate_epithelium.R -> 05_epithelial_compartments.R (defines compartments: meta.data$layers_rep)
#     -> 06_pseudobulk_dgea.R (THIS SCRIPT)
#
#  Inputs:
#   - <results_dir>/eso_epi_integratedObj.RData    EEC object with compartment labels (layers_rep)
#   - <pathway_maps_dir>/{GO_Biological_Process_2023,KEGG_2021_Human,
#       ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X}.map   2-column (term,gene) Enrichr gene-set maps
#
#  Outputs (to <results_dir>):
#   - pb_deg_res_layerRep_byRegion.RData                    ct.genes + qlf.res (by-location DEGs)
#   - pbObj_layerRep_byRegion.RData                         pseudobulk DGEList (for 11_clinical_associations.R)
#   - ssc_eso.epi.pb_deg.single_df.xlsx                     by-location DEG tables
#   - ssc_eso.epi.pb_deg.single_df.combinedL.xlsx           combined-location DEG tables
#   - PseudoDEG_enrich_Superficial_*.tsv                    GO/KEGG/TF over-representation
#   - <suppdata_dir>/proximal_distal_correlation.tsv, <suppdata_dir>/compartment_logfc_jitter.tsv   figure source data
#   - assorted .pdf figure panels
#
#  Usage:
#   Rscript 06_pseudobulk_dgea.R   (run from the repository root so the relative paths resolve)
#
#  Notes:
#   - This script was originally run locally on the aggregated/pseudobulk data (its paths pointed
#     at a local machine); they are parameterized to relative paths below.
#   - Removed from the original working file (not part of the manuscript): SSc cutaneous-subtype
#     comparisons and a parallel limma/voom DEG analysis (both rebuttal-only), SSc-vs-GERD volcano
#     plots (no volcano panels in the paper), a downsampling/permutation DEG power analysis, the
#     differentiation-bin pseudobulk exploration, and a clinical/motility tail that duplicates
#     11_clinical_associations.R. The GO/KEGG/TF over-representation enrichment IS retained (it produces
#     the over-representation tables and is distinct from the GSEA in 07_gsea_enrichment.R).
####################################################################################################


####  LOAD LIBRARIES  ##############################################################################

suppressPackageStartupMessages({
  library(methods)
  library(Seurat)
  library(edgeR)
  library(ggplot2)
  library(ggpubr)
  library(ggpattern)
  library(RColorBrewer)
  library(pheatmap)
  library(gridExtra)
  library(ggplotify)
  library(cowplot)
  library(openxlsx)
  library(scales)
  library(stringr)
  library(clusterProfiler)   # enricher() for the GO/KEGG/TF over-representation
  library(pals)       # alphabet.colors() for the per-sample MDS palette
})


####  SET PARAMETERS  ##############################################################################

set.seed(0123456789)

## ---- Paths (relative to repository root) ----
results_dir  <- './results'                 # pseudobulk objects + DEG outputs
suppdata_dir <- './results/supporting_data'  # figure source-data tables (proximal_distal_correlation.tsv, etc.)
pathway_maps_dir <- './results/pathway_gene_maps'  # Enrichr-style gene-set maps (GO/KEGG/TF)

## ---- Metadata column names ----
cond <- 'condition'
loc  <- 'location'

## ---- Colors ----
cond_cols    <- c('#8156B3', '#e26b53', '#A8A39d')                          # SSc, GERD, HC
layer_cols_3 <- c('#4ca5b1', '#e9b85d', '#c72f4c')                          # basal, suprabasal, superficial
layer_cols_5 <- c('#4B86B8', '#7DC0A6', '#E9F5A3', '#E1BA6C', '#B73D4F')    # 5 compartments

## ---- Contrasts / compartments ----
comps      <- list(c('SSc', 'HC'), c('GERD', 'HC'), c('SSc', 'GERD'))
layers     <- c('basal', 'suprabasal', 'superficial')
layers_rep <- c('Basal', 'ProliferatingBasal', 'ProliferatingSuprabasal', 'Suprabasal', 'Superficial')

## ---- Pseudobulk gene-filter thresholds ----
MIN_CPM         <- 1      # expressed = CPM > MIN_CPM ...
MIN_SAMPLE_PROP <- 0.75   # ... in > MIN_SAMPLE_PROP of one condition's samples
FDR_CUT         <- 0.05
N_VAR_GENES     <- 2000   # : prox-distal correlation uses the 2,000 most variable genes (Methods)

## Disease-associated genes labeled in .
metallothioneins <- c('MT1A', 'MT1E', 'MT1F', 'MT1G', 'MT1H', 'MT1M', 'MT1X', 'MT2A')
disease_genes_labeled <- c('PTGES', 'MFGE8', 'FCGBP', 'BST2', 'CD44', 'APOBEC3A', 'C10orf99', 'LTB4R', 'ACKR3',
                 'HLA-B', 'CD74', 'TAP1', 'PSMB8', 'PSMB9', 'CLEC2B', 'SGK1', 'HBEGF',
                 metallothioneins, 'H19', 'MUC22', 'SLC8A1-AS1')


####  FUNCTIONS  ###################################################################################

loadRData <- function(fileName) {
  load(fileName)
  get(ls()[ls() != 'fileName'])
}


####  BUILD PSEUDOBULK OBJECT  #####################################################################

epi_set <- loadRData(paste0(results_dir, '/eso_epi_integratedObj.RData'))

## 2,000 most variable genes (Methods: the  proximal-vs-distal correlation uses the 2,000
## most variable genes). Computed inline from the RNA assay rather than loaded from a saved file.
epi_set <- FindVariableFeatures(epi_set, assay = 'RNA', selection.method = 'vst', nfeatures = N_VAR_GENES)
var_genes_all <- VariableFeatures(epi_set, assay = 'RNA')


## Normalize the compartment labels to single tokens (Basal, ProliferatingBasal, ...) and order them.
epi_set@meta.data$layers_rep <- str_to_title(gsub('replicating', 'proliferating', epi_set@meta.data$layers_rep))
epi_set@meta.data$layers_rep <- gsub(' ', '', epi_set@meta.data$layers_rep)
epi_set@meta.data$layers_rep <- factor(epi_set@meta.data$layers_rep,
                                       levels = names(table(epi_set@meta.data$layers_rep))[c(1, 2, 3, 5, 4)])

## Aggregate raw counts to pseudobulk (sample x compartment) and build the DGEList.
epi_aggro <- AggregateExpression(epi_set, group.by = c('orig.ident', 'layers_rep'),
                                 return.seurat = TRUE, assays = 'RNA', slot = 'counts')
pb_obj <- DGEList(epi_aggro@assays$RNA@counts)

## Sample-level metadata parsed from the pseudobulk column names ("<sample>_<compartment>").
pb_obj$samples$sample    <- rownames(pb_obj$samples)
pb_obj$samples$condition <- ifelse(grepl('^I', pb_obj$samples$sample), 'GERD',
                            ifelse(grepl('^P', pb_obj$samples$sample), 'SSc', 'HC'))
pb_obj$samples$layer     <- word(rownames(pb_obj$samples), 2, sep = '\\_')
pb_obj$samples$location  <- ifelse(grepl('D_', pb_obj$samples$sample), 'Distal', 'Proximal')
pb_obj$samples$group     <- paste(pb_obj$samples$condition, pb_obj$samples$layer, sep = '_')

## Grouping factors. NOTE: `layer` is defined here (it is used by the gene filter immediately
## below). In the original working script this factor was not declared until ~700 lines later,
## so the filter only ran because `layer` happened to already exist in the session.
samples     <- factor(sub('\\-.*$', '', pb_obj$samples$sample))
sample_cols <- sample(alphabet.colors(26), 20)
condition   <- factor(pb_obj$samples$condition, levels = c('SSc', 'GERD', 'HC'))
location    <- factor(pb_obj$samples$location, levels = c('Proximal', 'Distal'))
layer       <- factor(pb_obj$samples$layer, levels = names(table(pb_obj$samples$layer))[c(1, 2, 3, 5, 4)])


####  GENE FILTER  #################################################################################
## Per compartment, keep genes expressed (CPM > MIN_CPM) in > MIN_SAMPLE_PROP of at least one
## condition's samples, and drop single-/two-sample outliers (driven by one donor).

count_cpm  <- cpm(pb_obj)
count_lcpm <- cpm(pb_obj, log = TRUE)

ct.genes <- lapply(setNames(levels(layer), levels(layer)),
                   function(ct) {
                     ct.df.l   <- as.data.frame(count_lcpm[, layer == ct])
                     ct.ssc.df  <- as.data.frame(count_cpm[, layer == ct & condition == 'SSc'])
                     ct.gerd.df <- as.data.frame(count_cpm[, layer == ct & condition == 'GERD'])
                     ct.hc.df   <- as.data.frame(count_cpm[, layer == ct & condition == 'HC'])
                     ct.ssc.genes  <- rownames(ct.ssc.df[rowSums(ct.ssc.df > MIN_CPM) > ncol(ct.ssc.df) * MIN_SAMPLE_PROP, ])
                     ct.gerd.genes <- rownames(ct.gerd.df[rowSums(ct.gerd.df > MIN_CPM) > ncol(ct.gerd.df) * MIN_SAMPLE_PROP, ])
                     ct.hc.genes   <- rownames(ct.hc.df[rowSums(ct.hc.df > MIN_CPM) >= ncol(ct.hc.df) * MIN_SAMPLE_PROP, ])

                     ct.outlier.genes  <- rownames(ct.df.l[rowSums(ct.df.l > (rowMeans(ct.df.l) + 3 * rowSds(as.matrix(ct.df.l)))) == 1, ])
                     ct.outlier.genes2 <- rownames(ct.df.l[rowSums(ct.df.l > (rowMeans(ct.df.l) + 3 * rowSds(as.matrix(ct.df.l)))) == 2, ])
                     for (gene in ct.outlier.genes2) {
                       top_two <- sub('\\-.*$', '', names(ct.df.l[gene, order(as.numeric(ct.df.l[gene, ]), decreasing = TRUE)][1:2]))
                       if (top_two[1] == top_two[2]) {
                         ct.outlier.genes <- c(ct.outlier.genes, gene)
                       }
                     }
                     ct.genes <- unique(c(ct.ssc.genes, ct.gerd.genes, ct.hc.genes))
                     ct.genes <- ct.genes[!ct.genes %in% ct.outlier.genes]
                     return(ct.genes)
                   })

pb_obj <- pb_obj[unique(unlist(ct.genes)), , keep = FALSE]


####  MDS OVERVIEW  ################################################################################

x_range <- c(-3, 6.2)
pdf(file = paste0(results_dir, '/PseudoDEG_epi_layerRep_MDS.pdf'), width = 11, height = 11)
par(mfrow = c(2, 2))
plotMDS(pb_obj, pch = 21, bg = sample_cols[samples],   main = 'MDS by Sample', cex = 1.5, xlim = x_range)
legend('bottomright', legend = levels(samples), pch = 21, pt.bg = sample_cols, cex = 0.65)
plotMDS(pb_obj, pch = 21, bg = layer_cols_5[layer],    main = 'MDS by Epithelial Compartment', cex = 1.5, xlim = x_range)
legend('bottomright', legend = levels(layer), pch = 21, pt.bg = layer_cols_5, cex = 0.8)
plotMDS(pb_obj, pch = 21, bg = cond_cols[condition],   main = 'MDS by Condition', cex = 1.5, xlim = x_range)
legend('bottomright', legend = levels(condition), pch = 21, pt.bg = cond_cols, cex = 0.8)
plotMDS(pb_obj, pch = 21, bg = c(6, 13)[location],     main = 'MDS by Biopsy Location', cex = 1.5, xlim = x_range)
legend('bottomright', legend = levels(location), pch = 21, pt.bg = c(6, 13), cex = 0.8)
dev.off()


####  DGEA BY COMPARTMENT x LOCATION (edgeR QL-F)  #################################################
## Canonical pseudobulk DEGs: condition contrasts within each compartment x location.

pb_obj$samples$group <- paste(pb_obj$samples$condition, pb_obj$samples$layer, pb_obj$samples$location, sep = '_')
group  <- factor(pb_obj$samples$group, levels = unique(pb_obj$samples$group))
design <- model.matrix(~ 0 + group)
colnames(design) <- gsub('group', '', colnames(design))

contrast_comps <- c(); contrast_names <- c()
for (l in layers_rep) {
  for (lo in rev(levels(location))) {
    for (c in comps) {
      contrast_names <- c(contrast_names, paste(c[1], c[2], l, lo, sep = '_'))
      contrast_comps <- c(contrast_comps, paste0(paste(c[1], l, lo, sep = '_'), '-', paste(c[2], l, lo, sep = '_')))
    }
  }
}
contrasts <- makeContrasts(contrasts = contrast_comps, levels = design)
colnames(contrasts) <- contrast_names

qlf.res <- list()
for (i in 1:ncol(contrasts)) {
  ct <- strsplit(colnames(contrasts)[i], '_')[[1]][3]
  if (i == 1 | ct != strsplit(colnames(contrasts)[max(i - 1, 1)], '_')[[1]][3]) {
    pb_obj.ct <- pb_obj[ct.genes[[ct]], , keep.lib.sizes = FALSE]
    pb_obj.ct <- calcNormFactors(pb_obj.ct)
    pb_obj.ct <- estimateDisp(pb_obj.ct, design, robust = TRUE)
  }
  fit <- glmQLFit(pb_obj.ct, design, robust = TRUE)
  qlf <- glmQLFTest(fit, contrast = contrasts[, i])
  qlf$comparison      <- colnames(contrasts)[i]
  qlf$table$p_adj     <- p.adjust(qlf$table$PValue, method = 'fdr')
  qlf.res[[names(contrasts[1, ])[i]]] <- qlf$table
}

save(ct.genes, qlf.res, file = paste0(results_dir, '/pb_deg_res_layerRep_byRegion.RData'))
## Also save the pseudobulk DGEList itself (sample metadata + filtered counts) for 11_clinical_associations.R.
save(pb_obj, file = paste0(results_dir, '/pbObj_layerRep_byRegion.RData'))

## Single combined sheet of by-location DEG tables.
x <- qlf.res
for (i in names(x)) {
  x[[i]]$Gene        <- rownames(x[[i]])
  x[[i]]$Comparison  <- paste(word(i, 1, sep = '_'), word(i, 2, sep = '_'), sep = '_')
  x[[i]]$Compartment <- gsub('Proliferating', 'Prolif', word(i, 3, sep = '_'))
  x[[i]]$Location    <- gsub('Distal', 'Dist', gsub('Proximal', 'Prox', word(i, 4, sep = '_')))
  x[[i]]$logFC  <- signif(x[[i]]$logFC, 4)
  x[[i]]$logCPM <- signif(x[[i]]$logCPM, 4)
  x[[i]]$F      <- signif(x[[i]]$F, 4)
  x[[i]]$PValue <- signif(x[[i]]$PValue, 3)
  x[[i]]$p_adj  <- signif(x[[i]]$p_adj, 3)
  x[[i]] <- x[[i]][order(x[[i]]$PValue), c(6:9, 1:5)]
}
x <- as.data.frame(do.call(rbind, x))
write.xlsx(x, paste0(results_dir, '/ssc_eso.epi.pb_deg.single_df.xlsx'), withFilter = TRUE, firstRow = TRUE,
           colWidths = 'auto', headerStyle = createStyle(textDecoration = 'BOLD'))
rm(x)


####  SSc-only / GERD-only / shared DEG classification  #################################
## Classify each gene per compartment x location by significance in the SSc-vs-HC and GERD-vs-HC
## contrasts. The per-category counts (printed below) are the source of the  stacked bars
## (the bar chart itself was drawn in Illustrator from these counts).

comp_cols <- lapply(setNames(as.vector(outer(layers_rep, levels(location), paste, sep = '_')),
                             as.vector(outer(layers_rep, levels(location), paste, sep = '_'))),
                    function(l) {
                      factor(ifelse(qlf.res[[paste('SSc', 'HC', l, sep = '_')]]$p_adj < FDR_CUT,
                              ifelse(qlf.res[[paste('GERD', 'HC', l, sep = '_')]]$p_adj < FDR_CUT, 'cornflowerblue', cond_cols[1]),
                              ifelse(qlf.res[[paste('GERD', 'HC', l, sep = '_')]]$p_adj < FDR_CUT, cond_cols[2], 'grey90')),
                             levels = c('grey90', cond_cols[1], cond_cols[2], 'cornflowerblue'))
                    })

## DEG-count summary per compartment x location ( source counts).
for (l in names(comp_cols)) {
  tab <- table(comp_cols[[l]])
  cat(l, '| shared(SSc&GERD):', tab['cornflowerblue'],
      '| SSc-only:', tab[cond_cols[1]], '| GERD-only:', tab[cond_cols[2]], '\n')
}


####  GERD-vs-HC vs SSc-vs-HC log2FC scatter  ###########################################
## Shown for the Superficial compartment (Proximal & Distal) in the paper; all compartments are
## plotted here. Points colored by the comp_cols classification; pink line = regression through 0.

pdf(file = paste0(results_dir, '/PseudoDEG_epi_scatterFC_SScGERD_HCloc_5layer_genes.pdf'), width = 21, height = 8)
xrange <- c(-6.7, 6.7); yrange <- c(-6.7, 6.7)
par(mfrow = c(2, 5))
for (lo in levels(location)) {
  for (l in layers_rep) {
    temp_data <- merge(qlf.res[[paste('GERD', 'HC', l, lo, sep = '_')]],
                       qlf.res[[paste('SSc', 'HC', l, lo, sep = '_')]], by = 'row.names')
    model   <- lm(I(logFC.y) ~ 0 + logFC.x, data = temp_data)
    mod_seq <- seq(-8, 8, length.out = 100)
    preds   <- predict(model, newdata = data.frame(logFC.x = mod_seq), interval = 'confidence')
    plot(qlf.res[[paste('GERD', 'HC', l, lo, sep = '_')]][, 'logFC'],
         qlf.res[[paste('SSc', 'HC', l, lo, sep = '_')]][, 'logFC'],
         type = 'n', ylab = '', xlab = '', xlim = xrange, ylim = yrange, xaxt = 'n', yaxt = 'n')
    polygon(c(rev(mod_seq), mod_seq), c(rev(preds[, 3]), preds[, 2]), col = alpha('lightpink', 0.5), border = NA)
    abline(0, model$coefficients, lty = 3, col = 'pink')
    par(new = TRUE)
    plot(qlf.res[[paste('GERD', 'HC', l, lo, sep = '_')]][order(comp_cols[[paste(l, lo, sep = '_')]]), 'logFC'],
         qlf.res[[paste('SSc', 'HC', l, lo, sep = '_')]][order(comp_cols[[paste(l, lo, sep = '_')]]), 'logFC'],
         pch = 16, col = as.character(sort(comp_cols[[paste(l, lo, sep = '_')]])), cex = 0.85,
         ylab = 'SSc vs. HC (logFC)', xlab = 'GERD vs. HC (logFC)', main = paste(lo, l), xlim = xrange, ylim = yrange)
    abline(0, 1, lty = 2, col = 'grey80')
    text(-5.3, 6, paste('r =', round(cor(qlf.res[[paste('GERD', 'HC', l, lo, sep = '_')]]$logFC,
                                         qlf.res[[paste('SSc', 'HC', l, lo, sep = '_')]]$logFC, method = 'spearman'), 2)))
    text(6.7, -5.5, paste0('m = ', round(model$coefficients, 2), '\n',
                           'r\u00b2 = ', round(summary(model)$r.squared, 2)), col = 'palevioletred3', adj = c(1, 1))
  }
}
dev.off()


####  OVER-REPRESENTATION ENRICHMENT (Superficial; GO / KEGG / TF)  ################################
## Over-representation (clusterProfiler::enricher) of the Superficial-compartment pseudobulk DEGs
## against Enrichr-style gene-set maps, separately for Proximal and Distal. The ENCODE/ChEA TF map
## yields the transcription-factor enrichment reported in the paper (e.g. IRF1, MYC, E2F4, NFE2L2)
## Uses the by-location qlf.res computed above (run here, before the
## combined-location DGEA below overwrites qlf.res).
##
## FLAG: the SSc-vs-HC TF call uses tf = FALSE while the GERD-vs-HC and SSc-vs-GERD TF calls use
## tf = TRUE (preserved exactly as in the original). The tf flag controls whether the TF gene-set
## background is restricted to expressed genes, so it affects the SSc-vs-HC TF enrichment p-values;
## confirm this asymmetry was intended.

## Gene-set maps: 2-column (term, gene) tables exported from Enrichr. Place under pathway_maps_dir.
annot_go   <- read.table(paste0(pathway_maps_dir, '/GO_Biological_Process_2023.map'), header = FALSE, sep = '\t', quote = '')
annot_kegg <- read.table(paste0(pathway_maps_dir, '/KEGG_2021_Human.map'), header = FALSE, sep = '\t', quote = '')
annot_tf   <- read.table(paste0(pathway_maps_dir, '/ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.map'), header = FALSE, sep = '\t', quote = '')

## Over-representation per location for one compartment (ct) and one contrast (comp). d selects
## all / up- (+) / down- (-) regulated DEGs; tf restricts the TF map to expressed-gene targets.
SEA <- function(deg.list, annot, ct, d = c('all', '+', '-'), tf = FALSE, comp = 'SSc_GERD') {
  d <- match.arg(d)
  enrich.res <- lapply(setNames(paste(ct, levels(location), sep = '_'), paste(ct, levels(location), sep = '_')),
                       function(grp) {
    background_genes <- ct.genes[[ct]]
    test_dat   <- deg.list[[paste0(comp, '_', grp)]]
    test_genes <- rownames(test_dat[test_dat$p_adj <= 0.05, ])
    if (d == '+') {
      test_genes <- test_genes[test_genes %in% rownames(test_dat[test_dat$logFC > 0, ])]
    } else if (d == '-') {
      test_genes <- test_genes[test_genes %in% rownames(test_dat[test_dat$logFC < 0, ])]
    }
    if (length(test_genes) == 0) return(NULL)
    if (tf) annot <- annot_tf[word(annot_tf[, 1], 1) %in% background_genes, ]
    ego <- enricher(gene = test_genes, universe = background_genes, pAdjustMethod = 'fdr',
                    pvalueCutoff = 1, qvalueCutoff = 1, TERM2GENE = annot)
    ego_df <- as.data.frame(ego@result)
    if (nrow(ego_df) == 0) return(NULL)
    ego_df$p.adjust <- p.adjust(ego_df$pvalue, method = 'bonferroni')
    ego_df <- ego_df[, c(1, 3:ncol(ego_df))]
    rownames(ego_df) <- NULL
    return(ego_df)
  })
  enrich.res[!sapply(enrich.res, is.null)]
}

## Write a SEA result list (keyed by "<compartment>_<location>") to a single TSV.
write_sea <- function(res, fname) {
  if (length(res) == 0) return(invisible())
  out <- do.call(rbind, lapply(names(res), function(g) { d <- res[[g]]; d$group <- g; d }))
  write.table(out, paste0(results_dir, '/', fname), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')
}

sea_results <- list(
  GO_SSc_HC          = SEA(qlf.res, annot_go,   'Superficial',      comp = 'SSc_HC'),
  GO_SSc_HC_up       = SEA(qlf.res, annot_go,   'Superficial', '+', comp = 'SSc_HC'),
  GO_SSc_HC_down     = SEA(qlf.res, annot_go,   'Superficial', '-', comp = 'SSc_HC'),
  GO_SSc_GERD        = SEA(qlf.res, annot_go,   'Superficial',      comp = 'SSc_GERD'),
  GO_SSc_GERD_up     = SEA(qlf.res, annot_go,   'Superficial', '+', comp = 'SSc_GERD'),
  GO_SSc_GERD_down   = SEA(qlf.res, annot_go,   'Superficial', '-', comp = 'SSc_GERD'),
  KEGG_SSc_GERD      = SEA(qlf.res, annot_kegg, 'Superficial'),
  KEGG_SSc_GERD_up   = SEA(qlf.res, annot_kegg, 'Superficial', '+'),
  KEGG_SSc_GERD_down = SEA(qlf.res, annot_kegg, 'Superficial', '-'),
  TF_SSc_HC          = SEA(qlf.res, annot_tf,   'Superficial', tf = FALSE, comp = 'SSc_HC'),   # see FLAG (tf = FALSE)
  TF_GERD_HC         = SEA(qlf.res, annot_tf,   'Superficial', tf = TRUE,  comp = 'GERD_HC'),
  TF_SSc_GERD        = SEA(qlf.res, annot_tf,   'Superficial', tf = TRUE,  comp = 'SSc_GERD'))

for (nm in names(sea_results)) write_sea(sea_results[[nm]], paste0('PseudoDEG_enrich_Superficial_', nm, '.tsv'))


####  DGEA BY COMPARTMENT (combined location)  ####################################################
## Location-pooled DEGs, used for the jitter plots and the combined table.

pb_obj$samples$group <- paste(pb_obj$samples$condition, pb_obj$samples$layer, sep = '_')
group  <- factor(pb_obj$samples$group, levels = unique(pb_obj$samples$group))
design <- model.matrix(~ 0 + group)
colnames(design) <- gsub('group', '', colnames(design))
pb_obj <- estimateDisp(pb_obj, design, robust = TRUE)

contrast_comps <- c(); contrast_names <- c()
for (c in comps) {
  for (l in layers_rep) {
    contrast_names <- c(contrast_names, paste(c[1], c[2], l, sep = '_'))
    contrast_comps <- c(contrast_comps, paste0(paste(c[1], l, sep = '_'), '-', paste(c[2], l, sep = '_')))
  }
}
contrasts <- makeContrasts(contrasts = contrast_comps, levels = design)
colnames(contrasts) <- contrast_names

qlf.res <- lapply(setNames(1:ncol(contrasts), names(contrasts[1, ])),
                  function(i) {
                    pb_obj.ct <- pb_obj[ct.genes[[strsplit(colnames(contrasts)[i], '_')[[1]][3]]], , keep.lib.sizes = FALSE]
                    pb_obj.ct <- calcNormFactors(pb_obj.ct)
                    fit <- glmQLFit(pb_obj.ct, design, robust = TRUE)
                    qlf <- glmQLFTest(fit, contrast = contrasts[, i])
                    qlf$comparison  <- colnames(contrasts)[i]
                    qlf$table$p_adj <- p.adjust(qlf$table$PValue, method = 'fdr')
                    return(qlf$table)
                  })

x <- qlf.res
for (i in names(x)) {
  x[[i]]$Gene        <- rownames(x[[i]])
  x[[i]]$Comparison  <- paste(word(i, 1, sep = '_'), word(i, 2, sep = '_'), sep = '_')
  x[[i]]$Compartment <- gsub('Proliferating', 'Prolif', word(i, 3, sep = '_'))
  x[[i]]$logFC  <- signif(x[[i]]$logFC, 4)
  x[[i]]$logCPM <- signif(x[[i]]$logCPM, 4)
  x[[i]]$F      <- signif(x[[i]]$F, 4)
  x[[i]]$PValue <- signif(x[[i]]$PValue, 3)
  x[[i]]$p_adj  <- signif(x[[i]]$p_adj, 3)
  x[[i]] <- x[[i]][order(x[[i]]$PValue), c(6:8, 1:5)]
}
x <- as.data.frame(do.call(rbind, x))
write.xlsx(x, paste0(results_dir, '/ssc_eso.epi.pb_deg.single_df.combinedL.xlsx'), withFilter = TRUE, firstRow = TRUE,
           colWidths = 'auto', headerStyle = createStyle(textDecoration = 'BOLD'))
rm(x)


####  Per-compartment log2FC jitter + percent-DE  ################################
comp_cols_jit <- c(cond_cols[1:2], '#B26183')
for (c in comps) {
  strip_df <- lapply(layers_rep, function(i) {
    tmp_df <- qlf.res[[paste(c[1], c[2], i, sep = '_')]]
    tmp_df$bin <- word(i, 1, sep = '\\_')
    rownames(tmp_df) <- paste(rownames(tmp_df), i, sep = '_')
    return(tmp_df)
  })
  strip_df <- as.data.frame(do.call(rbind, strip_df))
  strip_df$bin <- factor(strip_df$bin, levels = layers_rep)

  prop_degs <- sapply(layers_rep, function(i) {
    xn <- paste(c[1], c[2], i, sep = '_')
    round(nrow(qlf.res[[xn]][qlf.res[[xn]]$p_adj < FDR_CUT, ]) / nrow(qlf.res[[xn]]), 5)
  })

  pdf(file = paste0(results_dir, '/PseudoDEG_epiLayerRep_jitterFC.', c[1], '_', c[2], '.pdf'), width = 8, height = 5)
  par(mar = c(4.5, 4.5, 2, 4.5))
  stripchart(logFC ~ bin, strip_df[strip_df$p_adj >= FDR_CUT, ], method = 'jitter', ylim = c(-8, 6), cex = 0.5,
             cex.axis = 0.8, cex.lab = 1.1, pch = 16, jitter = 0.3, vertical = TRUE, col = 'grey90',
             main = paste(c[1], 'vs.', c[2], 'DEGs'), cex.main = 1.2, xlab = '', xaxt = 'n', yaxt = 'n',
             offset = 1, xlim = c(0.6, 4.6), at = seq(1, length.out = 5, by = 0.8))
  axis(1, at = seq(1, length.out = 5, by = 0.8), cex.axis = 0.8,
       labels = c('Basal', 'Proliferating\nBasal', 'Proliferating\nSuprabasal', 'Suprabasal', 'Superficial'))
  axis(2, at = seq(-6, 6, 2), labels = seq(-6, 6, 2), las = 2, cex.axis = 0.8)
  lines(x = seq(1, length.out = 5, by = 0.8), y = prop_degs * (14) - 8, col = 'cornflowerblue', lwd = 1.75)
  points(x = seq(1, length.out = 5, by = 0.8), y = prop_degs * (14) - 8, bg = 'cornflowerblue', pch = 23, cex = 1.05)
  par(new = TRUE)
  stripchart(logFC ~ bin, strip_df[strip_df$p_adj < FDR_CUT, ], method = 'jitter', ylim = c(-8, 6), cex = 0.5, offset = 1,
             xaxt = 'n', yaxt = 'n', ylab = '', pch = 16, jitter = 0.3, vertical = TRUE,
             col = comp_cols_jit[which(comps %in% list(c))], xlim = c(0.6, 4.6), at = seq(1, length.out = 5, by = 0.8))
  axis(4, at = seq(-8, 6, by = 2.8), labels = seq(0, 100, length = 6), col = 'cornflowerblue', lwd = 0, lwd.ticks = 1,
       las = 2, cex.axis = 0.8, col.axis = 'cornflowerblue')
  mtext('Percent of Genes DE', side = 4, padj = 5, las = 3, cex = 1.1, col = 'cornflowerblue')
  dev.off()
}

write.table(strip_df[, c(1, 5, 6)], paste0(suppdata_dir, '/compartment_logfc_jitter.tsv'),
            quote = FALSE, row.names = TRUE, col.names = TRUE, sep = '\t')


####  Proximal-vs-Distal expression correlation  #######################################
## Per sample, correlate proximal vs distal pseudobulk profiles (over variable genes) within each
## compartment; summarize by condition. Source data -> proximal_distal_correlation.tsv.

count_norm   <- as.data.frame(cpm(pb_obj, log = TRUE))
count_norm    <- count_norm[var_genes_all, ]   # var_genes_all = 2,000 most variable genes (computed above)

## Detailed proximal-distal correlation heatmaps per compartment (supplementary).
heatmaps <- lapply(setNames(layers_rep, layers_rep), function(layer) {
  set1 <- grepl(paste0('P_', layer), names(count_norm))
  set2 <- grepl(paste0('D_', layer), names(count_norm))
  correlation_matrix <- cor(count_norm[, set1], count_norm[, set2], use = 'pairwise.complete.obs')
  colnames(correlation_matrix) <- c('HC_6', 'HC_5', 'HC_4', 'HC_3', 'HC_2', 'HC_1', 'GERD_4', 'GERD_3', 'GERD_2', 'GERD_1',
                                    'SSc_10', 'SSc_9', 'SSc_8', 'SSc_7', 'SSc_6', 'SSc_5', 'SSc_4', 'SSc_3', 'SSc_2', 'SSc_1')
  rownames(correlation_matrix) <- colnames(correlation_matrix)[c(1:10, 12:20)]
  breaks <- seq(0.8, 1, by = 0.001)
  as.ggplot(pheatmap(correlation_matrix, cluster_rows = FALSE, cluster_cols = FALSE, legend = TRUE, main = layer,
                     breaks = breaks, color = colorRampPalette((brewer.pal(9, 'RdPu')))(length(breaks) - 1), border_color = NA))
})
pdf(file = paste0(results_dir, '/Dist_mat_corrHeatmaps_varGenes_all2k.pdf'), width = 22, height = 4)
do.call(grid.arrange, c(heatmaps, ncol = 5))
dev.off()

## : per-sample proximal-distal correlation, summarized by condition (boxplots).
pdf(file = paste0(results_dir, '/Prox_Dist_varGenes_all2k.pdf'), width = 16, height = 4)
par(mfrow = c(1, 5))
invisible(lapply(setNames(layers_rep, layers_rep), function(l) {
  set1 <- grepl(paste0('P_', l), names(count_norm))
  set2 <- grepl(paste0('D_', l), names(count_norm))
  correlation_matrix <- cor(count_norm[, set1], count_norm[, set2], use = 'pairwise.complete.obs')
  corrs <- c()
  for (i in rownames(correlation_matrix)) corrs <- c(corrs, correlation_matrix[i, gsub('P_', 'D_', i)])
  names(corrs) <- gsub('P_', '', rownames(correlation_matrix))
  plot_list <- list(corrs[1:6], corrs[7:10], corrs[11:19])   # HC (6), GERD (4), SSc (9)
  plot.new()
  rect(par('usr')[1], par('usr')[3], par('usr')[2], par('usr')[4], col = '#fcfaf7')
  par(new = TRUE)
  boxplot(plot_list, col = rev(cond_cols), main = l, names = c('HCs', 'GERD', 'SSc'), pch = 21, cex = 1.5,
          ylim = c(0.865, 0.995), staplewex = 0, bg = rev(cond_cols), lty = 1, outline = FALSE, range = 99)
}))
dev.off()

##  source data.
corr_df <- do.call(rbind, lapply(setNames(layers_rep, layers_rep), function(l) {
  set1 <- grepl(paste0('P_', l), names(count_norm))
  set2 <- grepl(paste0('D_', l), names(count_norm))
  correlation_matrix <- cor(count_norm[, set1], count_norm[, set2], use = 'pairwise.complete.obs')
  corrs <- c()
  for (i in rownames(correlation_matrix)) corrs <- c(corrs, correlation_matrix[i, gsub('P_', 'D_', i)])
  names(corrs) <- gsub('P_', '', rownames(correlation_matrix))
  data.frame(layer = l, sample = names(corrs), correlation = corrs,
             condition = c(rep('HC', 6), rep('GERD', 4), rep('SSc', 9)))
}))
rownames(corr_df) <- NULL
write.table(corr_df, paste0(suppdata_dir, '/proximal_distal_correlation.tsv'), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')


####  DISEASE-GENE EXPRESSION BOXPLOTS (superficial)  #############################################
## Per-gene pseudobulk expression (Superficial compartment) by condition and location, for the
## genes labeled in .

count_norm <- as.data.frame(cpm(pb_obj, log = TRUE))
plot_list  <- list(); i <- 1
for (gene in disease_genes_labeled) {
  data <- cbind.data.frame(gene = as.numeric(t(count_norm[gene, ])), condition, layer, location)
  data <- data[data$layer == 'Superficial', c(1, 2, 4)]
  data$condition <- factor(data$condition, c('HC', 'GERD', 'SSc'))
  rng <- c(min(data$gene), max(data$gene))
  plot_list[[i]] <- ggplot(data, aes(x = factor(condition), gene, fill = condition, pattern = location)) +
    geom_boxplot_pattern(outlier.size = 0.4, pattern_spacing = 0.03, pattern_density = 0.05,
                         position = position_dodge(0.85), width = 0.7) +
    scale_pattern_manual(values = c('stripe', 'none')) +
    theme(legend.position = 'none', axis.title = element_blank(), axis.text.x = element_blank(),
          panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
          panel.background = element_rect(fill = '#FCFAF7'), axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 6), plot.title = element_text(size = 8, face = 4, hjust = 0.5),
          panel.border = element_rect(color = 'black', fill = NA, size = 1)) +
    scale_fill_manual(values = rep(rev(cond_cols), 2)) + labs(title = gene) +
    scale_y_continuous(limits = c(rng[1] - (rng[2] - rng[1]) * 0.02, rng[2] + (rng[2] - rng[1]) * 0.12))
  i <- i + 1
}
pdf(file = paste0(results_dir, '/PseudoDEG_epiLayers_genes.pdf'), width = 12, height = 7)
cowplot::plot_grid(plotlist = plot_list, nrow = 4)
dev.off()
