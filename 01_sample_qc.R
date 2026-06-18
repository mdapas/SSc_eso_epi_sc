# ==============================================================================
# Single-cell sample- and cell-level quality control
#
# Part of the analysis for:
#   Dapas et al. (2025) "Cellular and molecular dysregulation of the esophageal
#   epithelium in systemic sclerosis."
#
# Purpose
#   For each sequenced sample, build a Seurat object from the Cell Ranger
#   filtered matrix and apply per-cell quality control:
#     - remove cells with excessive total counts (> 100,000)
#     - remove apoptotic/lysed cells via miQC (MT% + posterior probability)
#     - remove predicted doublets via scDblFinder
#   QC diagnostic plots and pre-/post-QC summary tables are written per sample.
#   (Corresponds to the Methods "Single-cell sample and cell quality control"
#   section.)
#
# Pipeline position
#   Upstream of integration. Run once per Cell Ranger output directory; the
#   resulting per-sample Seurat objects feed 02_integrate_allcells.R.
#
# Inputs
#   <data_dir>/<sample>/filtered_feature_bc_matrix/   (Cell Ranger output)
#
# Outputs (written under each sample's analysis/sample_QC/ directory)
#   seuratObj_<sample>.RData    QC-passing Seurat object
#   QC_<sample>.pdf             scatter QC diagnostics
#   Doublets_<sample>.pdf       doublet diagnostics
#   miQC_<sample>.pdf           miQC model fit
# Outputs (written under data_dir)
#   QC_summary_all.tsv          pass/fail counts per sample
#   raw_metrics_all.tsv         pre-QC summary metrics
#   qc_metrics_all.tsv          post-QC summary metrics
#
# Usage
#   Rscript --vanilla 01_sample_qc.R -d <data_dir>
# ==============================================================================


####  LOAD LIBRARIES  ##############################################################################

suppressPackageStartupMessages({
  library(methods)
  library(Seurat)
  library(ggplot2)
  library(scater)
  library(scales)
  library(scDblFinder)
  library(cowplot)
  library(SeuratWrappers)
  library(flexmix)
  library(optparse)
  library(RColorBrewer)
})


####  SET PARAMETERS  ##############################################################################

options(scipen = 10000)
set.seed(123456789)

# Metadata column names used throughout (defined once and referenced by the
# QC helper functions below).
count    <- "nCount_RNA"
features <- "nFeature_RNA"
mt       <- "percent.mt"
rb       <- "percent.rb"
qc       <- "qc"

# QC thresholds
HI_COUNT_CUT <- 100000   # maximum total UMI counts per cell
MIQC_CUT     <- 0.75     # miQC posterior probability of being a compromised cell
MT_CUT       <- 5        # maximum mitochondrial read percentage
MIN_FEATURES <- 200      # minimum unique genes per cell (CreateSeuratObject filter)


####  FUNCTIONS  ###################################################################################

get_qc <- function(sc_obj) {
  ## Return a per-cell vector flagging which QC metric (if any) a cell failed.
  ## A cell is labelled 'Pass' only if it clears all thresholds; doublets
  ## (flagged as NA by the miQC merge below) are labelled 'Doublet'.
  qc_vector <- rep("Pass", ncol(sc_obj))
  qc_vector <- ifelse(sc_obj@meta.data[[count]] > HI_COUNT_CUT, "High_nCount", qc_vector)
  qc_vector <- ifelse(sc_obj@meta.data[["miQC.keep"]] == "discard",
                      ifelse(qc_vector == "Pass",
                             ifelse(sc_obj@meta.data[[mt]] > MT_CUT, "miQC", qc_vector),
                             paste0(qc_vector, ", miQC")),
                      qc_vector)
  qc_vector[is.na(qc_vector)] <- "Doublet"
  return(qc_vector)
}


get_qc_plots <- function(sc_obj) {
  ## Build a list of QC diagnostic plots for a sample:
  ##   [[1]] 2x2 grid of feature-scatter plots coloured by pass/fail
  ##   [[2]] doublet UMAP + nFeature violin
  ##   [[3]] UMAP coloured by pass/fail
  pass_total   <- sum(sc_obj[[qc]] == "Pass")
  pass_percent <- round(pass_total / dim(sc_obj)[2], 3) * 100
  pass_label   <- paste0("\nPass\n(", pass_total, ", ", pass_percent, "%)\n")
  fail_total   <- sum(sc_obj[[qc]] != "Pass")
  fail_percent <- round(fail_total / dim(sc_obj)[2], 3) * 100
  fail_label   <- paste0("Fail\n(", fail_total, ", ", fail_percent, "%)")
  sc_obj[["pass_n"]] <- ifelse(sc_obj[[qc]] == "Pass", pass_label, fail_label)
  sc_obj@meta.data$pass_n <- factor(sc_obj@meta.data$pass_n, levels = c(pass_label, fail_label))
  sc_obj[["pass"]] <- ifelse(sc_obj[[qc]] == "Pass", "Pass", "Fail")
  sc_obj@meta.data$pass <- factor(sc_obj@meta.data$pass, levels = c("Pass", "Fail"))

  n_singlets    <- sum(sc_obj@meta.data$scDblFinder.class == "singlet")
  n_doublets    <- sum(sc_obj@meta.data$scDblFinder.class == "doublet")
  singlet_label <- paste0("\nSinglets\n(", n_singlets, ", ", round(n_singlets / dim(sc_obj)[2], 3) * 100, "%)\n")
  doublet_label <- paste0("Doublets\n(", n_doublets, ", ", round(n_doublets / dim(sc_obj)[2], 3) * 100, "%)")
  sc_obj[["doublets"]] <- ifelse(sc_obj@meta.data$scDblFinder.class == "singlet", singlet_label, doublet_label)
  sc_obj@meta.data$doublets <- factor(sc_obj@meta.data$doublets, levels = c(singlet_label, doublet_label))
  sc_obj[["scDblFinder.class"]] <- ifelse(sc_obj@meta.data$scDblFinder.class == "singlet", "Singlet", "Doublet")
  sc_obj@meta.data$scDblFinder.class <- factor(sc_obj@meta.data$scDblFinder.class, levels = c("Singlet", "Doublet"))

  # Transient SCT/PCA/UMAP just for visualising doublets and QC pass/fail
  sc_obj <- SCTransform(sc_obj, verbose = FALSE)
  sc_obj <- RunPCA(sc_obj, npcs = 20, verbose = FALSE)
  sc_obj <- RunUMAP(sc_obj, dims = 1:20, verbose = FALSE)

  return(list(
    plot_grid(
      FeatureScatter(sc_obj, feature1 = count, feature2 = features, group.by = "pass",
                     plot.cor = FALSE, pt.size = 0.1, cols = c("#F8766D", "gray")) +
        theme(legend.position = "none", axis.title = element_text(size = 10),
              axis.text = element_text(size = 9), title = element_text(size = 11, face = "bold")),
      FeatureScatter(sc_obj, feature1 = features, feature2 = mt, group.by = "pass",
                     plot.cor = FALSE, pt.size = 0.1, cols = c("#F8766D", "gray")) +
        theme(legend.position = c(.75, .9), axis.title = element_text(size = 10),
              axis.text = element_text(size = 9), title = element_text(size = 11, face = "bold"),
              legend.title = element_blank(), legend.text = element_text(size = 8)) +
        guides(color = guide_legend(override.aes = list(size = 2))),
      FeatureScatter(sc_obj, feature1 = rb, feature2 = mt, group.by = "pass",
                     plot.cor = FALSE, pt.size = 0.1, cols = c("#F8766D", "gray")) +
        theme(legend.position = "none", axis.title = element_text(size = 10),
              axis.text = element_text(size = 9), title = element_text(size = 10, face = "bold")),
      FeatureScatter(sc_obj, feature1 = features, feature2 = rb, group.by = "pass",
                     plot.cor = FALSE, pt.size = 0.1, cols = c("#F8766D", "gray")) +
        theme(legend.position = "none", axis.title = element_text(size = 10),
              axis.text = element_text(size = 9), title = element_text(size = 11, face = "bold"))),
    plot_grid(
      DimPlot(sc_obj, reduction = "umap", group.by = "doublets", cols = c("#F8766D", "turquoise1")) +
        theme(plot.title = element_blank()),
      VlnPlot(sc_obj, features = "nFeature_RNA", group.by = "scDblFinder.class",
              pt.size = 0.1, cols = c("#F8766D", "turquoise1")) +
        theme(legend.position = "none", plot.title = element_blank(), axis.title.x = element_blank()) +
        ylab("N Features"),
      ncol = 1),
    DimPlot(sc_obj, reduction = "umap", group.by = "pass_n", cols = c("#F8766D", "gray")) +
      ggtitle(" ")))
}


####  OPTIONS  #####################################################################################

option_list <- list(
  make_option(c("-d", "--dir"), action = "store", type = "character", default = NULL,
              help = "Root directory containing per-sample Cell Ranger subdirectories")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


####  MAIN  ########################################################################################

# Root directory where per-sample subdirectories are located
root_dir <- opt$dir

# Sample names = subdirectory names under the root directory
samples <- list.dirs(root_dir, recursive = FALSE, full.names = FALSE)

# Accumulators for per-sample QC summaries
cell_qc_df <- data.frame(Filter = character())
metric_cols <- c("sample", "cells",
                 paste0("med_", count), paste0("min_", count), paste0("max_", count),
                 paste0("med_", features), paste0("min_", features), paste0("max_", features),
                 paste0("med_", mt), paste0("min_", mt), paste0("max_", mt),
                 paste0("med_", rb), paste0("min_", rb), paste0("max_", rb))
sample_preQC_df  <- setNames(data.frame(matrix(ncol = 14, nrow = 0)), metric_cols)
sample_postQC_df <- sample_preQC_df

for (i in seq_along(samples)) {
  sample <- samples[i]
  print(sample)

  # Per-sample output directory
  output_dir <- paste(root_dir, sample, "analysis/sample_QC", sep = "/")
  if (!file.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  data.10x <- Read10X(data.dir = paste(root_dir, sample, "filtered_feature_bc_matrix", sep = "/"))

  print("  creating Seurat object...")
  sc_obj <- CreateSeuratObject(counts = data.10x, project = sample, min.features = MIN_FEATURES)
  rm(data.10x)
  sc_obj[[mt]] <- PercentageFeatureSet(sc_obj, pattern = "^MT-")
  sc_obj[[rb]] <- PercentageFeatureSet(sc_obj, pattern = "^RP[SL]")

  print("  performing QC...")

  # Predict doublets (dbr.sd controls the prediction threshold; see Methods)
  sc_obj <- as.Seurat(scDblFinder(as.SingleCellExperiment(sc_obj), dbr.sd = 0.01))
  # Drop the auxiliary scDblFinder score columns added to meta.data, keeping
  # scDblFinder.class. NOTE: column positions depend on the scDblFinder version.
  sc_obj@meta.data[, c(6, 8, 9, 10)] <- NULL

  # Apply miQC to derive sample-specific MT cutoffs, fit on the non-doublet set
  sc_obj_nonDbl <- subset(sc_obj, subset = scDblFinder.class == "singlet")
  sc_obj_nonDbl <- RunMiQC(sc_obj_nonDbl, percent.mt = mt, nFeature_RNA = features,
                           posterior.cutoff = MIQC_CUT, model.slot = "flexmix_model",
                           backup.option = "percent", backup.percent = 10)

  sc_obj@meta.data[rownames(sc_obj@meta.data) %in% rownames(sc_obj_nonDbl@meta.data),
                   c("miQC.probability", "miQC.keep")] <-
    sc_obj_nonDbl@meta.data[, c("miQC.probability", "miQC.keep")]

  # miQC model fit plot (only if miQC produced non-trivial probabilities)
  if (sum(sc_obj_nonDbl@meta.data$miQC.probability) > 0) {
    pdf(file = paste0(output_dir, "/miQC_", sample, ".pdf"), width = 6, height = 5)
    myPalette <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))
    sc <- scale_colour_gradientn(colours = myPalette(100), limits = c(0, 1))
    print(PlotMiQC(sc_obj_nonDbl, color.by = "miQC.probability") + sc +
            theme(axis.title = element_text(size = 10), axis.text = element_text(size = 9),
                  title = element_text(size = 10, face = "bold"), legend.text = element_text(size = 8),
                  legend.title = element_text(size = 9), plot.title = element_text(hjust = 0.5)) +
            ggtitle(sample))
    dev.off()
  }
  rm(sc_obj_nonDbl)

  # Pre-QC summary metrics
  sample_preQC_df[i, ] <- c(sample, ncol(sc_obj),
                            median(sc_obj[[count]][, 1]),    min(sc_obj[[count]][, 1]),    max(sc_obj[[count]][, 1]),
                            median(sc_obj[[features]][, 1]), min(sc_obj[[features]][, 1]), max(sc_obj[[features]][, 1]),
                            median(sc_obj[[mt]][, 1]),       min(sc_obj[[mt]][, 1]),       max(sc_obj[[mt]][, 1]),
                            median(sc_obj[[rb]][, 1]),       min(sc_obj[[rb]][, 1]),       max(sc_obj[[rb]][, 1]))

  # Determine QC pass/fail per cell
  sc_obj[[qc]] <- get_qc(sc_obj)

  temp_df <- data.frame(table(sc_obj[[qc]]))
  names(temp_df) <- c("Filter", sample)
  cell_qc_df <- merge(cell_qc_df, temp_df, by = "Filter", all = TRUE)

  # QC diagnostic plots
  qc_plots <- get_qc_plots(sc_obj)

  pdf(file = paste0(output_dir, "/QC_", sample, ".pdf"), width = 12, height = 5.5)
  print(plot_grid(ggdraw() + draw_label(paste(sample, "- Initial QC"), fontface = "bold"),
                  plot_grid(qc_plots[[1]], qc_plots[[3]], ncol = 2),
                  ncol = 1, rel_heights = c(0.1, 1), scale = 0.98))
  dev.off()

  pdf(file = paste0(output_dir, "/Doublets_", sample, ".pdf"), width = 7, height = 7)
  print(plot_grid(ggdraw() + draw_label(paste(sample, "- Predicted Doublets"), fontface = "bold"),
                  qc_plots[[2]], ncol = 1, rel_heights = c(0.1, 1), scale = 0.98))
  dev.off()

  # Apply QC filter and record post-QC metrics
  sc_obj <- subset(sc_obj, subset = qc == "Pass")

  sample_postQC_df[i, ] <- c(sample, ncol(sc_obj),
                             median(sc_obj[[count]][, 1]),    min(sc_obj[[count]][, 1]),    max(sc_obj[[count]][, 1]),
                             median(sc_obj[[features]][, 1]), min(sc_obj[[features]][, 1]), max(sc_obj[[features]][, 1]),
                             median(sc_obj[[mt]][, 1]),       min(sc_obj[[mt]][, 1]),       max(sc_obj[[mt]][, 1]),
                             median(sc_obj[[rb]][, 1]),       min(sc_obj[[rb]][, 1]),       max(sc_obj[[rb]][, 1]))

  # Save the QC-passing object for this sample
  save(sc_obj, file = paste0(output_dir, "/seuratObj_", sample, ".RData"))
}

# Transpose the pass/fail count table so samples are rows
cell_qc_df <- setNames(data.frame(t(cell_qc_df[, -1])), cell_qc_df[, 1])

# Write cohort-level QC summaries
write.table(cell_qc_df,       paste(root_dir, "QC_summary_all.tsv", sep = "/"), quote = FALSE, row.names = TRUE,  col.names = TRUE, sep = "\t")
write.table(sample_preQC_df,  paste(root_dir, "raw_metrics_all.tsv", sep = "/"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(sample_postQC_df, paste(root_dir, "qc_metrics_all.tsv", sep = "/"),  quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
