####################################################################################################
#  09_cellchat_communication.R
#
#  Cell-cell communication analysis (CellChat) across the esophageal cell types and epithelial
#  compartments.
#
#  Paper:
#   Dapas et al. (2025) "Cellular and molecular dysregulation of the esophageal epithelium in
#   systemic sclerosis." JCI Insight.
#
#  Purpose / figures produced:
#   - Differential interaction-strength heatmap, SSc vs. GERD & HC (symmetric), plus the
#              pairwise differential heatmaps (SSc vs HC, SSc vs GERD, GERD vs HC).
#   - Differential incoming/outgoing signaling-role scatter, per pairwise comparison.
#   - Incoming/outgoing signaling-role heatmaps for SSc (SSc-enriched pathways).
#   - Info flow: rankNet pathway comparison (the signaling pathways enriched in SSc).
#   - Per-condition ligand-receptor interaction tables (.csv).
#
#  Pipeline position:
#   03_annotate_celltypes.R (cell typing) + 05_epithelial_compartments.R (epithelial compartments)
#     -> 09_cellchat_communication.R (THIS SCRIPT)
#   (Originally this lived at the end of a newer version of 03_annotate_celltypes.R; it is
#    pulled into its own script here because it depends on downstream products - the epithelial
#    compartments. The reorganization pass can relocate it.)
#
#  Inputs:
#   - <results_dir>/All_integratedObj.RData        annotated all-cell object (meta.data$types, $t_cell_subtype from 03)
#   - <results_dir>/eso_epi_integratedObj.RData    EEC object with compartment labels (layers_rep)
#
#  Outputs (to <results_dir>):
#   - cellchat_perCondition_list.RData, cellchat_merged_3cond.RData, cellchat_merged_SSc_gerdHC.RData
#   - .pdf panels; cellchat_rankNet*.pdf
#   - cellchat_interactions.<cond>.csv
#
#  Usage:
#   Rscript 09_cellchat_communication.R   (run from the repository root; memory-intensive)
#
#  Notes:
#   - CellChat groups are `types_t_cells`: the base cell types with the epithelium split into its
#     five differentiation compartments and the lymphocytes (L) split into CD4+ T / CD8+ T.
#   - The CD4+/CD8+ split is taken from 03_annotate_celltypes.R (sc_obj$t_cell_subtype), produced by
#     the CD45+ immune re-integration there.
#   - The epithelium is downsampled to 100,000 cells (DOWNSAMPLE_N) before running CellChat, as in
#     the original, to balance the very large epithelial fraction against the other cell types.
####################################################################################################


####  LOAD LIBRARIES  ##############################################################################

suppressPackageStartupMessages({
  library(methods)
  library(Seurat)
  library(CellChat)
  library(patchwork)
  library(plyr)
  library(dplyr)
  library(RColorBrewer)
  library(circlize)
  library(ComplexHeatmap)
  library(R.utils)        # capitalize()
})


####  SET PARAMETERS  ##############################################################################

set.seed(123)

results_dir <- './results'

conditions   <- c('HC', 'GERD', 'SSc')
cond_cols    <- c('#A8A39d', '#e26b53', '#8156B3')   # HC, GERD, SSc

DOWNSAMPLE_N <- 100000   # epithelial cells retained
NBOOT        <- 40       # CellChat permutations
MIN_CELLS    <- 10       # filterCommunication threshold

## Epithelial compartment labels (with line breaks as used for the figure axes).
epi_clusters <- c('Basal', 'Proliferating\nbasal', 'Proliferating\nsuprabasal', 'Suprabasal', 'Superficial')

## Display order of the 14 CellChat groups and their colors.
group_levels <- c('Ep: basal', 'Ep: prolif.\nbasal', 'Ep: prolif.\nsuprabasal', 'Ep: suprabasal',
                  'Ep: superficial', 'CD4+ T', 'CD8+ T', 'Mc', 'My', 'SMG', 'En', 'F', 'P', 'SMCs')
epi_cols_5 <- c('#3288BD', '#66C2A5', '#E6F598', '#e9b85d', '#c72f4c')
all_cols   <- c(epi_cols_5, '#AD00D4', '#4682B4', '#EEA2AD', '#CD96CD', '#76EE00',
                'slategray3', 'brown4', '#FF7F50', 'orange4')
names(all_cols) <- group_levels

## SSc-enriched signaling pathways shown in the  role heatmaps.
ssc_deg_paths <- c('COLLAGEN', 'LAMININ', 'NOTCH', 'FN1', 'Prostaglandin', 'GRN', 'JAM', 'CD46', 'EGF',
                   'GALECTIN', 'PECAM1', 'ADGRA', '12oxoLTB4', 'VISFATIN', 'Cholesterol', 'CD96',
                   'NECTIN', 'CDH5', 'MPZ', 'EPHA', 'CysLTs', 'CD40', 'ANGPT', 'EPHB', 'CD39')


####  FUNCTIONS  ###################################################################################

loadRData <- function(fileName) {
  load(fileName)
  get(ls()[ls() != 'fileName'])
}

## Standard CellChat workflow for one (already-subset) Seurat object.
run_cellchat <- function(obj, group.by = 'types_t_cells') {
  cc <- createCellChat(obj, group.by = group.by, assay = 'RNA')
  cc@DB <- CellChatDB.human
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = 'triMean', nboot = NBOOT)
  cc <- filterCommunication(cc, min.cells = MIN_CELLS)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc, slot.name = 'netP')
  cc
}

## Custom differential-interaction heatmap (vendored from the original; extends CellChat's
## netVisual_heatmap with ylim.top/ylim.right/row.show/col.show controls for the symmetric differential-interaction heatmap).
my_netVisual_heatmap <- function(object, comparison = c(1, 2), measure = c("count", "weight"), signaling = NULL, slot.name = c("netP", "net"), color.use = NULL, color.heatmap = NULL,
                                 title.name = NULL, width = NULL, height = NULL, ylim.top = NULL, ylim.right = NULL,
                                 font.size = 8, font.size.title = 10, cluster.rows = FALSE, cluster.cols = FALSE,
                                 sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, row.show = NULL, col.show = NULL) {
  if (!is.null(measure)) measure <- match.arg(measure)
  slot.name <- match.arg(slot.name)
  if (class(object@net[[1]]) == "list") {
    message("Do heatmap based on a merged object \n")
    if (is.null(color.heatmap)) color.heatmap <- c('#2166ac', '#b2182b')
    obj1 <- object@net[[comparison[1]]][[measure]]
    obj2 <- object@net[[comparison[2]]][[measure]]
    net.diff <- obj2 - obj1
    if (measure == "count") {
      if (is.null(title.name)) title.name = "Differential number of interactions"
    } else if (measure == "weight") {
      if (is.null(title.name)) title.name = "Differential interaction strength"
    }
    legend.name = "Relative values"
  } else {
    message("Do heatmap based on a single object \n")
    if (is.null(color.heatmap)) color.heatmap <- "Reds"
    if (!is.null(signaling)) {
      prob <- slot(object, slot.name)$prob
      if (slot.name == "net") prob[object@net$pval > thresh] <- 0
      net.diff <- prob[, , signaling]
      if (is.null(title.name)) title.name = paste0(signaling, " signaling network")
      legend.name <- "Communication Prob."
    } else if (!is.null(measure)) {
      net.diff <- object@net[[measure]]
      if (measure == "count") {
        if (is.null(title.name)) title.name = "Number of interactions"
      } else if (measure == "weight") {
        if (is.null(title.name)) title.name = "Interaction strength"
      }
      legend.name <- title.name
    }
  }
  net <- net.diff
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source", "target")
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) sources.use <- rownames(net.diff)[sources.use]
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) targets.use <- rownames(net.diff)[targets.use]
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- rownames(net.diff)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0
  if (is.null(color.use)) color.use <- scPalette(ncol(net))
  names(color.use) <- colnames(net)
  color.use.row <- color.use
  color.use.col <- color.use
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    if (length(idx1) > 0) { net <- net[-idx1, ]; color.use.row <- color.use.row[-idx1] }
    if (length(idx2) > 0) { net <- net[, -idx2]; color.use.col <- color.use.col[-idx2] }
  }
  mat <- net
  if (!is.null(row.show)) { mat <- mat[row.show, , drop = FALSE]; color.use.row <- color.use.row[row.show] }
  if (!is.null(col.show)) { mat <- mat[, col.show, drop = FALSE]; color.use.col <- color.use.col[col.show] }
  if (min(mat) < 0) {
    lim <- max(abs(min(mat)), abs(max(mat)))
    color.heatmap.use = colorRamp3(c(-lim, 0, lim), c(color.heatmap[1], "#f7f7f7", color.heatmap[2]))
    colorbar.break <- c(-lim, 0, lim)
  } else {
    if (length(color.heatmap) == 3) {
      color.heatmap.use = colorRamp3(c(0, min(mat), max(mat)), color.heatmap)
    } else if (length(color.heatmap) == 2) {
      color.heatmap.use = colorRamp3(c(min(mat), max(mat)), color.heatmap)
    } else if (length(color.heatmap) == 1) {
      color.heatmap.use = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(100)
    }
    colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*", "\\1", min(mat, na.rm = T))) + 1), round(max(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*", "\\1", max(mat, na.rm = T))) + 1))
  }
  df.col <- data.frame(group = colnames(mat)); rownames(df.col) <- colnames(mat)
  df.row <- data.frame(group = rownames(mat)); rownames(df.row) <- rownames(mat)
  col_annotation <- HeatmapAnnotation(df = df.col, col = list(group = color.use.col), which = "column",
                                      show_legend = FALSE, show_annotation_name = FALSE, simple_anno_size = grid::unit(0.2, "cm"))
  row_annotation <- HeatmapAnnotation(df = df.row, col = list(group = color.use.row), which = "row",
                                      show_legend = FALSE, show_annotation_name = FALSE, simple_anno_size = grid::unit(0.2, "cm"))
  ha1 = rowAnnotation(Strength = anno_barplot(rowSums(abs(mat)), border = FALSE, gp = gpar(fill = color.use.row, col = color.use.row), add_numbers = FALSE, ylim = ylim.right), show_annotation_name = FALSE)
  ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(abs(mat)), border = FALSE, gp = gpar(fill = color.use.col, col = color.use.col), add_numbers = FALSE, ylim = ylim.top), show_annotation_name = FALSE)
  if (sum(abs(mat) > 0) == 1) { color.heatmap.use = c("white", color.heatmap.use) } else { mat[mat == 0] <- NA }
  ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", name = legend.name,
                bottom_annotation = col_annotation, left_annotation = row_annotation, top_annotation = ha2, right_annotation = ha1,
                cluster_rows = cluster.rows, cluster_columns = cluster.rows,
                row_names_side = "left", row_names_rot = 0, row_names_gp = gpar(fontsize = font.size), column_names_gp = gpar(fontsize = font.size),
                column_title = title.name, column_title_gp = gpar(fontsize = font.size.title), column_names_rot = 90,
                row_title = "Sources (Sender)", row_title_gp = gpar(fontsize = font.size.title), row_title_rot = 90,
                heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"), title_position = "leftcenter-rot",
                                            border = NA, legend_height = unit(20, "mm"), labels_gp = gpar(fontsize = 8), grid_width = unit(2, "mm")))
  return(ht1)
}


####  BUILD COMPOSITE CELL-TYPE LABEL  ############################################################
## types_t_cells = base cell types, with the epithelium split into compartments and the
## lymphocytes split into CD4+ T / CD8+ T.

sc_obj  <- loadRData(paste0(results_dir, '/All_integratedObj.RData'))
sc_obj$samples <- sc_obj$orig.ident

## Epithelial compartments from the EEC object.
epi_set <- loadRData(paste0(results_dir, '/eso_epi_integratedObj.RData'))
epi_set$layers_rep <- gsub(' ', '\n', capitalize(gsub('replicating', 'proliferating', epi_set$layers_rep)))
sc_obj$types_epi <- as.character(sc_obj$types)
sc_obj$types_epi[colnames(epi_set)] <- as.character(epi_set$layers_rep)
rm(epi_set); gc()

## Split lymphocytes into CD4+ / CD8+ T cells using the labels resolved in
## 03_annotate_celltypes.R (sc_obj$t_cell_subtype).
sc_obj$types_t_cells <- sc_obj$types_epi
has_tcell <- !is.na(sc_obj$t_cell_subtype)
sc_obj$types_t_cells[has_tcell] <- sc_obj$t_cell_subtype[has_tcell]

## Map internal codes to the display labels used for the communication groups, then order them.
display_map <- c('Basal' = 'Ep: basal', 'Proliferating\nbasal' = 'Ep: prolif.\nbasal',
                 'Proliferating\nsuprabasal' = 'Ep: prolif.\nsuprabasal', 'Suprabasal' = 'Ep: suprabasal',
                 'Superficial' = 'Ep: superficial')
sc_obj$types_t_cells <- ifelse(sc_obj$types_t_cells %in% names(display_map),
                               display_map[sc_obj$types_t_cells], sc_obj$types_t_cells)
sc_obj$types_t_cells <- factor(sc_obj$types_t_cells, levels = group_levels)


####  DOWNSAMPLE EPITHELIUM  ######################################################################
## Retain DOWNSAMPLE_N epithelial cells plus all non-epithelial cells.

epi_display    <- group_levels[1:5]
epi_cells      <- colnames(sc_obj)[sc_obj$types_t_cells %in% epi_display]
other_cells    <- colnames(sc_obj)[!sc_obj$types_t_cells %in% c(epi_display, NA) & !is.na(sc_obj$types_t_cells)]
epi_downsampled <- sample(epi_cells, size = min(DOWNSAMPLE_N, length(epi_cells)), replace = FALSE)
downsampled_obj <- subset(sc_obj, cells = c(other_cells, epi_downsampled))
downsampled_obj$types_t_cells <- droplevels(factor(downsampled_obj$types_t_cells, levels = group_levels))
rm(sc_obj); gc()


####  CELLCHAT PER CONDITION  #####################################################################

## Three-condition run (HC, GERD, SSc) ->  and pairwise differential heatmaps.
cellchat_list <- setNames(lapply(conditions, function(cond) {
  run_cellchat(downsampled_obj[, colnames(downsampled_obj)[downsampled_obj$condition == cond]])
}), conditions)
save(cellchat_list, file = paste0(results_dir, '/cellchat_perCondition_list.RData'))

cellchat_3 <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))
cellchat_3 <- computeNetSimilarityPairwise(cellchat_3, type = 'functional')
cellchat_3 <- rankNetPairwise(cellchat_3)
save(cellchat_3, file = paste0(results_dir, '/cellchat_merged_3cond.RData'))

## Two-group run (SSc vs. GERD & HC) -> symmetric .
downsampled_obj$condition_ssc <- ifelse(downsampled_obj$condition == 'SSc', 'SSc', 'GERD_HC')
cellchat_2list <- setNames(lapply(c('GERD_HC', 'SSc'), function(cond) {
  run_cellchat(downsampled_obj[, colnames(downsampled_obj)[downsampled_obj$condition_ssc == cond]])
}), c('GERD_HC', 'SSc'))
cellchat_2 <- mergeCellChat(cellchat_2list, add.names = names(cellchat_2list))
save(cellchat_2, file = paste0(results_dir, '/cellchat_merged_SSc_gerdHC.RData'))


####  Differential interaction heatmaps  ###############################################
## Symmetric SSc vs. GERD & HC (net.diff = SSc - GERD_HC).
pdf(file = paste0(results_dir, '/cellChat_diffInts.SSc_GERD-HC.symm.pdf'), width = 6, height = 5.5)
print(my_netVisual_heatmap(cellchat_2, measure = 'weight', comparison = c(1, 2),
                           ylim.top = c(0, 2), ylim.right = c(0, 3),
                           title.name = 'Differential interactions: SSc vs. GERD & HC', color.use = all_cols,
                           row.show = 1:14, col.show = 1:14))
dev.off()

## Pairwise differential heatmaps (HC=1, GERD=2, SSc=3).
for (cmp in list(c('SSc_HC', c(1, 3)), c('SSc_GERD', c(2, 3)), c('GERD_HC', c(1, 2)))) {
  nm <- cmp[[1]]; idx <- as.integer(c(cmp[[2]], cmp[[3]]))
  pdf(file = paste0(results_dir, '/cellChat_diffInts.', nm, '.pdf'), width = 6, height = 5.5)
  print(netVisual_heatmap(cellchat_3, measure = 'weight', comparison = idx,
                          title.name = paste('Differential interactions:', sub('_', ' vs. ', nm)),
                          color.use = all_cols))
  dev.off()
}


####  Differential signaling-role scatter  #############################################
for (cmp in list(c('SSc_HC', c(1, 3)), c('SSc_GERD', c(2, 3)), c('GERD_HC', c(1, 2)))) {
  nm <- cmp[[1]]; idx <- as.integer(c(cmp[[2]], cmp[[3]]))
  pdf(file = paste0(results_dir, '/cellChat_diffScatter.', nm, '.pdf'), width = 6, height = 5.5)
  print(netAnalysis_diff_signalingRole_scatter(cellchat_3, comparison = idx, color.use = all_cols))
  dev.off()
}


####  SSc signaling-role heatmaps (incoming / outgoing)  ###############################
pdf(file = paste0(results_dir, '/cellChat_SSc_signalingOutvsIn.pdf'), width = 10, height = 8)
p_in  <- netAnalysis_signalingRole_heatmap(cellchat_list[['SSc']], signaling = ssc_deg_paths, pattern = 'incoming',
                                           title = 'SSc', width = 8, height = 16, color.use = all_cols,
                                           ylim.top = c(0, 8), ylim.right = c(0, 4))
p_out <- netAnalysis_signalingRole_heatmap(cellchat_list[['SSc']], signaling = ssc_deg_paths, pattern = 'outgoing',
                                           title = 'SSc', width = 8, height = 16, color.use = all_cols,
                                           ylim.top = c(0, 8), ylim.right = c(0, 4))
print(p_in + p_out)
dev.off()


####  INFO FLOW (rankNet)  ########################################################################
## Pathways differentially active between SSc and HC (information-flow comparison).

pdf(file = paste0(results_dir, '/cellChat_rankNet.SSc_HC.pdf'), width = 6, height = 8)
print(rankNet(cellchat_3, mode = 'comparison', comparison = c(1, 3), measure = 'weight',
              stacked = TRUE, cutoff.pvalue = 0.01, do.stat = TRUE, color.use = cond_cols[c(1, 3)]))
dev.off()


####  PER-CONDITION INTERACTION TABLES  ###########################################################

for (cond in conditions) {
  ints <- subsetCommunication(cellchat_list[[cond]])
  ints$condition <- cond
  write.csv(ints, paste0(results_dir, '/cellchat_interactions.', cond, '.csv'), row.names = FALSE)
}
