####################################################################################################
#  07_gsea_enrichment.R
#
#  Gene-set enrichment analysis (GSEA) of the epithelial-compartment single-cell DEGs.
#
#  Paper:
#   Dapas et al. (2025) "Cellular and molecular dysregulation of the esophageal epithelium in
#   systemic sclerosis." JCI Insight.
#
#  Purpose / figures produced:
#   - GSEA dotplot (MSigDB C2 canonical pathways) for SSc-vs-HC and GERD-vs-HC, by
#              compartment (basal / suprabasal / superficial) and biopsy location.
#   - GSEAtable.tsv   full GSEA result table (all compartments and contrasts).
#
#  Pipeline position:
#   05_epithelial_compartments.R (saves the per-compartment by-location single-cell DEG lists)
#     -> 07_gsea_enrichment.R (THIS SCRIPT)
#   (Note: the transcription-factor / GO over-representation enrichment that uses the *pseudobulk*
#    DEGs lives in 06_pseudobulk_dgea.R, not here, because it is coupled to that script's qlf.res.)
#
#  Inputs (from <results_dir>, produced by 05_epithelial_compartments.R):
#   - epi_degs_<comp>_layersRep_Proximal.RData
#   - epi_degs_<comp>_layersRep_Distal.RData
#       Named lists (one element per compartment) of FindMarkers tables, for
#       comp in {ssc_hc, gerd_hc, ssc_gerd}.
#
#  Outputs (to <results_dir>):
#   - epi_DEG.SSc_HC.GERD_HC.GSEA_layersRep.pdf    dotplot
#   - GSEAtable.tsv                                full GSEA results
#
#  Usage:
#   Rscript 07_gsea_enrichment.R   (run from the repository root)
####################################################################################################


####  LOAD LIBRARIES  ##############################################################################

suppressPackageStartupMessages({
  library(fgsea)
  library(dplyr)
  library(ggplot2)
  library(msigdbr)
  library(stringr)
  library(RColorBrewer)
  library(ggh4x)       # facet_grid2 / strip_themed
})


####  SET PARAMETERS  ##############################################################################

set.seed(123)

results_dir <- './results'

cond_cols    <- c('#8156B3', '#e26b53', '#A8A39d')   # SSc, GERD, HC
layer_cols_3 <- c('#4ca5b1', '#e9b85d', '#c72f4c')   # basal, suprabasal, superficial

## GSEA parameters (Methods).
GSEA_MIN_SIZE <- 15
GSEA_MAX_SIZE <- 1000
GSEA_N_PERM   <- 100000
PLOT_P_CUT    <- 0.001   # pathways shown in the  dotplot


####  FUNCTIONS  ###################################################################################

loadRData <- function(fileName) {
  load(fileName)
  get(ls()[ls() != 'fileName'])
}

## Build a ranked gene vector from a FindMarkers table for fgsea:
##   weight = -log10(p_val) * |avg_log2FC|, restricted to genes detected in >1% of cells in both
##   groups, with mitochondrial and ribosomal genes removed.
get_genes <- function(df) {
  df <- df[df$pct.1 > 0.01 & df$pct.2 > 0.01, ]
  df <- df[!is.na(df$p_val), ]
  df[df$p_val == 0, 'p_val'] <- 5e-324
  df$weight <- -log10(df$p_val) * abs(df$avg_log2FC)
  df <- df[order(df$weight, decreasing = TRUE), ]
  v <- df[, 'weight']
  names(v) <- rownames(df)
  v <- v[!is.na(v)]
  v <- v[!grepl('^MT-', names(v))]
  v <- v[!grepl('^RP',  names(v))]
  v <- v[!grepl('^MRP', names(v))]
  return(v)
}


####  LOAD DEG LISTS  ##############################################################################
## Per-compartment single-cell DEG lists, split by biopsy location.

deg_list_ssc_hc_p   <- loadRData(paste0(results_dir, '/epi_degs_ssc_hc_layersRep_Proximal.RData'))
deg_list_gerd_hc_p  <- loadRData(paste0(results_dir, '/epi_degs_gerd_hc_layersRep_Proximal.RData'))
deg_list_ssc_gerd_p <- loadRData(paste0(results_dir, '/epi_degs_ssc_gerd_layersRep_Proximal.RData'))
deg_list_ssc_hc_d   <- loadRData(paste0(results_dir, '/epi_degs_ssc_hc_layersRep_Distal.RData'))
deg_list_gerd_hc_d  <- loadRData(paste0(results_dir, '/epi_degs_gerd_hc_layersRep_Distal.RData'))
deg_list_ssc_gerd_d <- loadRData(paste0(results_dir, '/epi_degs_ssc_gerd_layersRep_Distal.RData'))


####  GENE SETS (MSigDB C2 canonical pathways)  ###################################################

c2_df <- msigdbr(species = 'Homo sapiens', category = 'C2')
c2_df <- c2_df[grepl('^CP', c2_df$gs_subcat), ]          # canonical pathways (Reactome, KEGG, NABA, ...)

fgsea_c2_sets <- c2_df %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_c2_sets <- fgsea_c2_sets[!duplicated(fgsea_c2_sets)]
fgsea_c2_sets <- lapply(fgsea_c2_sets, unique)


####  RUN GSEA  ####################################################################################
## fgsea per compartment, separately for each contrast and biopsy location.

run_gsea <- function(deg_list, comp_label, loc_label) {
  setNames(lapply(names(deg_list), function(layer) {
    res <- fgseaMultilevel(fgsea_c2_sets, stats = get_genes(deg_list[[layer]]),
                           minSize = GSEA_MIN_SIZE, maxSize = GSEA_MAX_SIZE,
                           scoreType = 'pos', nPermSimple = GSEA_N_PERM)
    res$comp  <- comp_label
    res$loc   <- loc_label
    res$layer <- layer
    res
  }), names(deg_list))
}

enrich_ssc_hc_p   <- run_gsea(deg_list_ssc_hc_p,   'SSc_HC',   'Proximal')
enrich_gerd_hc_p  <- run_gsea(deg_list_gerd_hc_p,  'GERD_HC',  'Proximal')
enrich_ssc_gerd_p <- run_gsea(deg_list_ssc_gerd_p, 'SSc_GERD', 'Proximal')
enrich_ssc_hc_d   <- run_gsea(deg_list_ssc_hc_d,   'SSc_HC',   'Distal')
enrich_gerd_hc_d  <- run_gsea(deg_list_gerd_hc_d,  'GERD_HC',  'Distal')
enrich_ssc_gerd_d <- run_gsea(deg_list_ssc_gerd_d, 'SSc_GERD', 'Distal')


####  GSEA dotplot  #####################################################################
## Non-proliferating compartments (basal, suprabasal, superficial = layersRep indices 1, 4, 5),
## SSc-vs-HC and GERD-vs-HC, pathways with raw p < PLOT_P_CUT.

enrich_plot_df <- data.frame()
for (i in c(1, 4, 5)) {
  for (df in list(enrich_ssc_hc_p[[i]], enrich_ssc_hc_d[[i]], enrich_gerd_hc_p[[i]], enrich_gerd_hc_d[[i]])) {
    enrich_plot_df <- rbind(enrich_plot_df, df)
  }
}

enrich_plot_df$padj    <- p.adjust(enrich_plot_df$pval, method = 'fdr')
enrich_plot_df         <- enrich_plot_df[enrich_plot_df$pval < PLOT_P_CUT, ]
enrich_plot_df$db      <- sub('\\_.*', '', enrich_plot_df$pathway)
enrich_plot_df$pathway <- str_extract(enrich_plot_df$pathway, '_.*')
enrich_plot_df$pathway <- gsub('_', ' ', substr(enrich_plot_df$pathway, 2, nchar(enrich_plot_df$pathway)))
enrich_plot_df$layer   <- factor(enrich_plot_df$layer, levels = c('basal', 'suprabasal', 'superficial'))
enrich_plot_df$comp    <- factor(enrich_plot_df$comp, levels = c('SSc_HC', 'GERD_HC'))
enrich_plot_df$loc     <- factor(enrich_plot_df$loc, levels = c('Proximal', 'Distal'))
enrich_plot_df$pathway <- factor(enrich_plot_df$pathway, levels = names(sort(table(enrich_plot_df$pathway))))

pdf(file = paste0(results_dir, '/epi_DEG.SSc_HC.GERD_HC.GSEA_layersRep.pdf'), width = 5.5, height = 6)
print(ggplot(enrich_plot_df, aes(x = comp, y = pathway, color = NES, size = -log10(pval))) +
        geom_point() +
        facet_grid2(loc ~ layer, space = 'free', scale = 'free_y',
                    strip = strip_themed(background_x = elem_list_rect(fill = layer_cols_3),
                                         background_y = elem_list_rect(fill = c('grey35', 'grey50')),
                                         text_x = element_text(face = 'bold', size = 7),
                                         text_y = element_text(face = 'bold', size = 7, color = 'white'))) +
        theme(text = element_text(size = 8),
              axis.text.x = element_text(face = 'bold', size = 7, colour = cond_cols[c(1, 2)]),
              axis.title = element_blank()) +
        scale_color_gradient(low = 'skyblue', high = 'navy'))
dev.off()


####  FULL GSEA RESULT TABLE  ######################################################################
## All compartments and all three contrasts (P & D), FDR-corrected. leadingEdge is a list column,
## so flatten it for writing.

gsea_all <- do.call(rbind, c(enrich_ssc_hc_p, enrich_ssc_hc_d, enrich_gerd_hc_p, enrich_gerd_hc_d,
                             enrich_ssc_gerd_p, enrich_ssc_gerd_d))
gsea_all$padj_global <- p.adjust(gsea_all$pval, method = 'fdr')
gsea_all$leadingEdge <- sapply(gsea_all$leadingEdge, paste, collapse = ',')
write.table(apply(gsea_all, 2, as.character), paste0(results_dir, '/GSEAtable.tsv'),
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')
