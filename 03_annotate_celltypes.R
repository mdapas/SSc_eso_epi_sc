# ==============================================================================
# Whole-tissue cell-type annotation and composition
#
# Part of the analysis for:
#   Dapas et al. (2025) "Cellular and molecular dysregulation of the esophageal
#   epithelium in systemic sclerosis."
#
# Purpose
#   Cluster the integrated whole-tissue object, identify and annotate the major
#   esophageal mucosal cell types (and the smaller stromal/immune populations
#   that share one mixed cluster), and summarize cell-type composition. Also
#   examines FLI1 expression across cell types (negligible in epithelium, low in
#   endothelium) and re-integrates the CD45+ immune compartment to resolve
#   T-cell (CD4+/CD8+) and myeloid subtypes.
#   (Corresponds to the Methods "Sample integration and cell type annotation".)
#
# Pipeline position
#   Downstream of 02_integrate_allcells.R. The annotated object it saves
#   (All_integratedObj.RData) is the entry point for the epithelial-cell
#   analyses (05_epithelial_compartments.R, etc.); the immune-subtype labels it
#   adds feed 09_cellchat_communication.R.
#
# Inputs
#   <results_dir>/eso_integratedRefObj.RData   integrated object (PCA computed)
#   <myeloid_ref_rds>                          external myeloid reference atlas
#                                              (for myeloid subtype label transfer)
#
# Outputs (written under results_dir)
#   All_integratedObj.RData     annotated whole-tissue object (with t_cell_subtype /
#                               myeloid_subtype immune labels)
#   eso_CD45_integratedObj.RData  re-integrated CD45+ immune compartment
#   eso_ct_markers.tsv          cell-type markers
#   eso_cellType_props.tsv      per-sample cell-type proportions
#   plus diagnostic / panel PDFs (UMAPs, module plots, violins, pie chart)
#
# NOTE on ordering
#   The original working script was run interactively and out of order (e.g.
#   cell-type markers were computed before the cell types were annotated). It
#   has been reordered into a single top-to-bottom pipeline: cluster -> identify
#   -> annotate -> save -> summarize. The computational steps are unchanged.
# ==============================================================================


####  LOAD LIBRARIES  ##############################################################################

suppressPackageStartupMessages({
  library(methods)
  library(Seurat)
  library(sctransform)
  library(ggplot2)
  library(glmGamPoi)
  library(cowplot)
  library(plyr)          # mapvalues() for cluster -> cell-type relabelling
  library(RColorBrewer)
  library(future)
  library(scCustomize)   # Stacked_VlnPlot()
  library(SingleR)               # myeloid subtype label transfer
  library(SingleCellExperiment)
})


####  SET PARAMETERS  ##############################################################################

set.seed(0123456789)
plan("multisession", workers = 2)

# Memory options (adjust to your allocation)
Sys.setenv("R_MAX_VSIZE" = 192000000000)
options(future.globals.maxSize = 192000 * 1024^2)   # ~ for a 160 GB RAM request

# --- Paths (edit to match your environment) -----------------------------------
results_dir <- "./results"

# External reference for myeloid subtype label transfer (a separate annotated myeloid
# atlas; cluster labels expected in its $Clusters2 metadata column).
myeloid_ref_rds <- "./data/references/myeloid_reference.RDS"

# Metadata column names
cond <- "condition"
loc  <- "location"

# Clustering parameters (Methods: top 40 PCs, SNN k=25)
N_PCS         <- 40
SNN_K         <- 25
CLUST_RES     <- 0.25
CLUSTER_COL   <- paste0("integrated_snn_res.", CLUST_RES)

# Cluster-12 (mixed small-population) re-clustering parameters (Methods: PCs=25, k=15)
CL12_PCS      <- 25
CL12_K        <- 15
CL12_RES      <- 0.5

# Major cell-type marker panel for the whole-tissue UMAP feature plots
major_markers <- c("KRT6A", "PTPRC", "CD3D", "LYZ",       # epithelial, CD45, T cell, myeloid
                   "PECAM1", "PDGFRA", "TAGLN", "TPSAB1")  # endothelial, fibroblast, smooth muscle, mast

# Stacked-violin marker panel, in row order
celltype_violin_markers <- c("KRT6A", "KRT14", "MUC5B", "PTPRC", "CD3D", "LYZ", "HLA-DRA",
                   "TPSAB1", "PECAM1", "PDGFRA", "PDGFRB", "RGS5", "ACTA2")

# Internal cell-type codes (display labels in the figures map En/Ep/L/My/Mc/GEp/
# F/P/SMc/SMG -> Endothelial/Epithelial/Lymphocytes/Myeloid/Mast cells/Glandular
# epithelial/Fibroblasts/Pericytes/Smooth muscle/Submucosal gland).
cell_types <- c("En", "Ep", "F", "GEp", "L", "Mc", "My", "P", "SMc", "SMG")

# Cell-type colours, ordered to match the alphabetical factor levels of cell_types
type_cols <- c("#CD96CD", "#8FBC8F", "#A52A2A", "#698B69", "#7AC5CD",
               "#4682B4", "#EEA2AD", "#DAA520", "#FF7F50", "#76EE00")

# Condition colours (HC, GERD, SSc) and order
cond_levels <- c("HC", "GERD", "SSc")
cond_cols   <- c("darkseagreen", "darkseagreen4", "lightpink")


####  FUNCTIONS  ###################################################################################

loadRData <- function(fileName) {
  ## Load a single object from an .RData file and return it.
  load(fileName)
  get(ls()[ls() != "fileName"])
}

module_plot <- function(obj, name, markers) {
  ## UMAP feature plot of an AddModuleScore signature for a marker set.
  marker_list <- list(markers)
  obj <- AddModuleScore(obj, features = marker_list, name = "mod_score", assay = "RNA")
  FeaturePlot(obj, features = "mod_score1", order = TRUE) & labs(x = "UMAP 1", y = "UMAP 2") &
    scale_colour_gradientn(colours = brewer.pal(n = 11, "Spectral")[11:1], name = "Avg. Gene\nExpression\n",
                           breaks = c(min(obj@meta.data$mod_score1), max(obj@meta.data$mod_score1)),
                           labels = c("Min", "Max"), guide = guide_colorbar(frame.colour = "black", ticks = FALSE)) &
    ggtitle(name, subtitle = paste(markers, collapse = ", ")) &
    theme(plot.title = element_text(face = "bold", size = 10),
          plot.subtitle = element_text(face = "italic", size = 8, hjust = 0.5),
          axis.text = element_blank(), axis.title = element_text(size = 8), axis.ticks = element_blank(),
          legend.title = element_text(size = 7, face = "bold"), legend.text = element_text(size = 6.5))
}


####  MAIN  ########################################################################################

## --- Load integrated object --------------------------------------------------
sc_obj <- loadRData(paste0(results_dir, "/eso_integratedRefObj.RData"))


## --- Diagnostic plots: PCs and per-cell QC on the UMAP -----------------------
pdf(file = paste0(results_dir, "/PC_plot_integratedAssay.pdf"), width = 5, height = 6)
print(ElbowPlot(sc_obj, ndims = 50))
dev.off()

for (feat in c("percent.rb", "percent.mt", "nFeature_RNA")) {
  pdf(file = paste0(results_dir, "/All_featurePlot_", feat, ".pdf"), width = 7, height = 6)
  print(FeaturePlot(sc_obj, features = feat, order = TRUE) &
          scale_colour_gradientn(colours = brewer.pal(n = 9, "YlGnBu")[3:9]))
  dev.off()
}

# PCs on the integrated UMAP (10 per page)
for (i in 0:4) {
  pc_start <- 10 * i + 1
  pc_end   <- 10 * i + 10
  pcs <- paste0("PC_", pc_start:pc_end)
  pdf(file = paste0(results_dir, "/All_UMAP_pcs", pc_start, "_", pc_end, ".pdf"), width = 14, height = 6)
  print(FeaturePlot(sc_obj, features = pcs, order = TRUE, ncol = 5) &
          scale_colour_gradientn(colours = brewer.pal(n = 11, "Spectral")[1:11]) &
          theme(plot.title = element_text(face = "bold", size = 10), legend.position = "none",
                axis.text = element_blank(), axis.title = element_text(size = 8), axis.ticks = element_blank(),
                legend.title = element_text(size = 7, face = "bold"), legend.text = element_text(size = 6.5)))
  dev.off()
}


## --- Major-marker UMAPs and per-sample distributions -------------------------
DefaultAssay(sc_obj) <- "RNA"
plan("default")
sc_obj <- NormalizeData(sc_obj)

pdf(file = paste0(results_dir, "/All_UMAP_markers.pdf"), width = 14, height = 6)
print(FeaturePlot(sc_obj, features = major_markers, order = TRUE, ncol = 4) &
        scale_colour_gradientn(colours = brewer.pal(n = 9, "YlGnBu")[3:9]))
dev.off()

samples <- names(table(sc_obj@meta.data$orig.ident))
dir.create(paste0(results_dir, "/All_UMAP_SampleDistributions"), showWarnings = FALSE)
for (sample in samples) {
  sc_obj@meta.data$sampleBin <- ifelse(sc_obj@meta.data$orig.ident == sample, 1, 0)
  pdf(file = paste0(results_dir, "/All_UMAP_SampleDistributions/All_dimPlot_", sample, ".pdf"), width = 5, height = 5)
  print(DimPlot(sc_obj, reduction = "umap", label = FALSE, group.by = "sampleBin", order = TRUE,
                cols = c("grey80", "black")) + ggtitle(sample) + theme(legend.position = "none"))
  dev.off()
}


## --- Cluster the integrated assay --------------------------------------------
plan("multisession", workers = 4)
DefaultAssay(sc_obj) <- "integrated"

sc_obj <- FindNeighbors(sc_obj, reduction = "pca", dims = 1:N_PCS, k.param = SNN_K, n.trees = 50)
sc_obj <- FindClusters(sc_obj, resolution = CLUST_RES, algorithm = 1, n.iter = 10, n.start = 10)

pdf(file = paste0(results_dir, "/All_UMAP_clusters_k", SNN_K, "_res", CLUST_RES, ".pdf"), width = 7, height = 6)
print(DimPlot(sc_obj, reduction = "umap", label = TRUE, group.by = CLUSTER_COL,
              order = as.character(30:0)) + ggtitle(""))
dev.off()


## --- Cell-type module-score plots --------------------
module_plots <- list(
  module_plot(sc_obj, "Epithelial",          c("KRT13", "KRT15", "KRT19", "S100A2", "S100A8")),
  module_plot(sc_obj, "Endothelial",         c("PECAM1", "CDH5", "ADGRL4", "EMCN")),
  module_plot(sc_obj, "Myeloid",             c("LYZ", "AIF1", "HLA-DRA")),
  module_plot(sc_obj, "Lymphocytes",         c("CD3D", "CD3E", "CD2", "TRBC2")),
  module_plot(sc_obj, "Mast cells",          c("TPSAB1", "TPSB2", "HPGDS", "CPA3")),
  module_plot(sc_obj, "Submucosal glands",   c("MUC5B", "TFF3", "LYZ", "C6orf58")),
  module_plot(sc_obj, "Glandular epithelial", c("KRT19", "KRT14", "C19orf33", "SLPI")),
  module_plot(sc_obj, "Fibroblasts",         c("PDGFRA", "LUM", "DCN", "COL1A2")),
  module_plot(sc_obj, "Pericytes",           c("KCNJ8", "RGS5", "PDGFRB", "HIGD1B")),
  module_plot(sc_obj, "Smooth muscle cells", c("ACTA2", "TAGLN", "MYH11")))

pdf(file = paste0(results_dir, "/All_modulePlots.pdf"), width = 15, height = 6)
print(plot_grid(plotlist = module_plots, ncol = 5))
dev.off()


## --- Resolve the mixed small-population cluster (cluster 12) ------------------
## The major cell types separate cleanly, but several smaller stromal/immune
## populations share a single cluster (cluster 12). Subset it, recompute PCs and
## re-cluster, then annotate the resulting subclusters by marker expression.
cl12_set <- subset(sc_obj, subset = !!sym(CLUSTER_COL) == "12")

cl12_set <- RunPCA(cl12_set, assay = "integrated")
cl12_set <- RunUMAP(cl12_set, assay = "integrated", dims = 1:CL12_PCS, n.neighbors = CL12_K,
                    min.dist = 0.3, spread = 1, repulsion.strength = 1,
                    negative.sample.rate = 10L, n.epochs = 1000)

# Representative markers used to identify the small populations
cl12_markers <- c("C19orf33", "TPSAB1", "CD3D", "MUC5B",   # glandular epi, mast, T cell, submucosal gland
                  "PDGFRA", "RGS5", "PDGFRB", "ACTA2", "MYH11")  # fibroblast, pericyte, smooth muscle
DefaultAssay(cl12_set) <- "RNA"
pdf(file = paste0(results_dir, "/Cl_12_UMAP_markers.pdf"), width = 6, height = 7.5)
print(FeaturePlot(cl12_set, features = cl12_markers, order = TRUE, ncol = 2) &
        scale_colour_gradientn(colours = brewer.pal(n = 9, "YlGnBu")[3:9]) &
        theme(axis.text = element_blank(), axis.title = element_text(size = 9), axis.ticks = element_blank()))
dev.off()

DefaultAssay(cl12_set) <- "integrated"
cl12_set <- FindNeighbors(cl12_set, reduction = "pca", dims = 1:CL12_PCS, k.param = CL12_K, n.trees = 5000)
cl12_set <- FindClusters(cl12_set, resolution = CL12_RES, algorithm = 1, n.iter = 50, n.start = 50)

# Subcluster -> small-population identity (established from marker expression)
#   0 - Mast cells               4 - Pericytes
#   1 - Glandular epithelial     5 - Smooth muscle cells
#   2 - Fibroblasts              6 - Submucosal glands
#   3 - T cells (lymphoid)
cl12_type_labels <- c("Mast cells", "Glandular\nepithelial cells", "Fibroblasts", "T cells",
                      "Pericytes", "SMCs", "Submucosal glands")
cl12_set@meta.data$types <- mapvalues(cl12_set@meta.data$seurat_clusters,
                                       from = as.character(0:6), to = cl12_type_labels)

pdf(file = paste0(results_dir, "/Cl_12_UMAP_types.pdf"), width = 6, height = 6)
print(DimPlot(cl12_set, reduction = "umap", label = TRUE, group.by = "types",
              cols = c("darkseagreen3", "lightpink2", "lavenderblush3", "steelblue1",
                       "khaki", "salmon1", "honeydew3")) +
        ggtitle("") + theme(legend.position = "none"))
dev.off()

cl12_module_plots <- list(
  module_plot(cl12_set, "Mast cells",          c("TPSAB1", "TPSB2", "HPGDS", "CPA3")),
  module_plot(cl12_set, "T cells",             c("CD3D", "CD2", "IL32", "KLRB1")),
  module_plot(cl12_set, "Submucosal glands",   c("MUC5B", "TFF3", "LYZ", "C6orf58")),
  module_plot(cl12_set, "Glandular epithelium", c("KRT19", "KRT14", "C19orf33", "SLPI")),
  module_plot(cl12_set, "Fibroblasts",         c("PDGFRA", "LUM", "DCN", "COL1A2")),
  module_plot(cl12_set, "Pericytes",           c("KCNJ8", "RGS5", "PDGFRB", "HIGD1B")),
  module_plot(cl12_set, "Smooth muscle cells", c("ACTA2", "TAGLN", "MYH11")))

pdf(file = paste0(results_dir, "/Cl_12_modulePlots.pdf"), width = 14, height = 6)
print(plot_grid(plotlist = cl12_module_plots, ncol = 4))
dev.off()


## --- Annotate cell types on the whole-tissue object --------------------------
## Major clusters are assigned directly; cluster-12 cells inherit the
## small-population identities resolved above. (Final labels are determined
## entirely by the explicit assignments below.)
sc_obj@meta.data$types <- as.character(sc_obj@meta.data$seurat_clusters)

# Major clusters
sc_obj@meta.data[sc_obj@meta.data$seurat_clusters %in% c(0, 1, 2, 3, 4, 6, 7, 8, 9, 13), "types"] <- "Ep"
sc_obj@meta.data[sc_obj@meta.data$seurat_clusters == 5,  "types"] <- "L"
sc_obj@meta.data[sc_obj@meta.data$seurat_clusters == 10, "types"] <- "My"
sc_obj@meta.data[sc_obj@meta.data$seurat_clusters == 11, "types"] <- "En"

# Cluster-12 small populations (from cl12 subclustering)
sc_obj@meta.data[rownames(cl12_set@meta.data[cl12_set@meta.data$seurat_clusters == 0, ]), "types"] <- "Mc"
sc_obj@meta.data[rownames(cl12_set@meta.data[cl12_set@meta.data$seurat_clusters == 1, ]), "types"] <- "GEp"
sc_obj@meta.data[rownames(cl12_set@meta.data[cl12_set@meta.data$seurat_clusters == 2, ]), "types"] <- "F"
sc_obj@meta.data[rownames(cl12_set@meta.data[cl12_set@meta.data$seurat_clusters == 3, ]), "types"] <- "L"
sc_obj@meta.data[rownames(cl12_set@meta.data[cl12_set@meta.data$seurat_clusters == 4, ]), "types"] <- "P"
sc_obj@meta.data[rownames(cl12_set@meta.data[cl12_set@meta.data$seurat_clusters == 5, ]), "types"] <- "SMc"
sc_obj@meta.data[rownames(cl12_set@meta.data[cl12_set@meta.data$seurat_clusters == 6, ]), "types"] <- "SMG"

pdf(file = paste0(results_dir, "/All_dimPlot_annotation.pdf"), width = 7, height = 6)
print(DimPlot(sc_obj, reduction = "umap", label = TRUE, group.by = "types", shuffle = TRUE, raster = FALSE,
              cols = type_cols, label.size = 3, pt.size = 0.0001) + ggtitle("") &
        theme(axis.text = element_blank(), axis.title = element_text(size = 9),
              axis.ticks = element_blank(), legend.text = element_text(size = 9)))
dev.off()

# nFeature distribution by annotated cell type
pdf(file = paste0(results_dir, "/All_violin_nFeature.pdf"), width = 8, height = 5)
print(ggplot(sc_obj@meta.data, aes(x = types, y = nFeature_RNA, fill = types)) +
        geom_violin(scale = "width") + labs(y = "nFeature_RNA", x = "Cell type") +
        ggtitle("Esophageal Cell Clusters") + theme_classic() +
        scale_fill_manual(values = type_cols) +
        theme(axis.title = element_text(face = "bold", size = 8), axis.ticks = element_blank(),
              legend.position = "none", panel.grid.major.y = element_line(),
              panel.grid.minor.y = element_line()))
dev.off()

# Save the annotated whole-tissue object
save(sc_obj, file = paste0(results_dir, "/All_integratedObj.RData"))


## --- Cell-type markers --------------------------------
DefaultAssay(sc_obj) <- "RNA"
Idents(sc_obj) <- "types"
ct_markers <- FindAllMarkers(sc_obj, logfc.threshold = 0.1, min.pct = 0.01, only.pos = TRUE)
write.table(ct_markers, paste0(results_dir, "/eso_ct_markers.tsv"),
            quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")


## --- Per-sample cell-type proportions ----------------------------------------
eso_cluster_df <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(eso_cluster_df) <- c("Sample", "Location", "Condition", "CellType", "nCells", "totCells", "Proportion")
i <- 1
for (sample in samples) {
  dat   <- sc_obj@meta.data[sc_obj@meta.data$orig.ident == sample, ]
  s_loc <- ifelse(substr(sample, nchar(sample), nchar(sample)) == "P", "Proximal", "Distal")
  s_con <- dat[1, cond]
  total <- nrow(dat)
  for (l in names(table(dat$types))) {
    nCells <- nrow(dat[dat$types == l, ])
    eso_cluster_df[i, ] <- c(sample, s_loc, as.character(s_con), l, nCells, total, round(nCells / total, 4))
    i <- i + 1
  }
}
write.table(eso_cluster_df, paste0(results_dir, "/eso_cellType_props.tsv"),
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")


## --- Composition summaries ( B-C) ------------------------------------

# Cell-type proportions table (for ordering the pie chart and stacked violin)
type_df <- data.frame(table(sc_obj@meta.data$types) / nrow(sc_obj@meta.data))
names(type_df) <- c("type", "prop")
type_df$percent <- paste0(round(type_df$prop * 100, 2), "%")
type_df$cols <- type_cols
type_df <- type_df[order(type_df$prop, decreasing = TRUE), ]
rownames(type_df) <- 1:nrow(type_df)
type_df$type <- factor(type_df$type, level = type_df$type[1:10])

# Pie chart of overall cell-type composition
pdf(file = paste0(results_dir, "/All_types_pie.pdf"), width = 6, height = 5)
print(ggplot(type_df, aes(x = "", y = prop, fill = type)) + geom_bar(stat = "identity", width = 1) +
        coord_polar(theta = "y") +
        scale_fill_manual(name = "Cell type",
                          labels = paste0(type_df$type, " (", type_df$percent, ")"), values = type_df$cols) +
        theme_void())
dev.off()

# Stacked violin of canonical markers by cell type
type_order <- c(1, 6, 10, 2, 3, 5, 4, 7, 8, 9)   # display order of cell types
sc_obj@meta.data$types <- factor(sc_obj@meta.data$types, levels = type_df$type[type_order])

pdf(file = paste0(results_dir, "/All_types_markers_violin.pdf"), width = 5, height = 10)
print(Stacked_VlnPlot(sc_obj, celltype_violin_markers, pt.size = 0, group.by = "types", raster = FALSE,
                      plot_spacing = 0, colors_use = type_df$cols[type_order], x_lab_rotate = TRUE) &
        theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()))
dev.off()

# Collapsed cell-type proportions by condition (epithelial, endothelial,
# lymphoid, myeloid, and all remaining types grouped as "Other")
prop_df <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(prop_df) <- c("Cell type", "Condition", "Percent")
i <- 1
for (condition_i in cond_levels) {
  total <- nrow(sc_obj@meta.data[sc_obj$condition == condition_i, ])
  for (type in c("En", "Ep", "L", "My", "Other")) {
    if (type == "Other") {
      n <- nrow(sc_obj@meta.data[sc_obj$condition == condition_i &
                                   sc_obj$types %in% c("F", "GEp", "Mc", "P", "SMc", "SMG"), ])
    } else {
      n <- nrow(sc_obj@meta.data[sc_obj$condition == condition_i & sc_obj$types == type, ])
    }
    prop_df[i, ] <- c(type, condition_i, round(n / total, 5))
    i <- i + 1
  }
}
prop_df$Percent   <- as.numeric(prop_df$Percent) * 100
prop_df$Condition <- factor(prop_df$Condition, levels = cond_levels)

pdf(file = paste0(results_dir, "/All_cellType_vs_condition.pdf"), width = 7, height = 8)
print(ggplot(prop_df, aes(x = Condition, y = Percent, fill = `Cell type`)) +
        geom_bar(stat = "identity") + theme_minimal() +
        geom_text(aes(label = paste(round(Percent, 1), "%")), position = position_stack(vjust = 0.5)) +
        scale_fill_manual("Cell type", values = c("plum3", "darkseagreen", "cadetblue3", "lightpink2", "grey")))
dev.off()


## --- FLI1 expression ( A-B) -----------------------------
## FLI1 is negligible in epithelial cells and expressed at low levels in
## endothelial cells, with no significant difference between conditions.
pdf(file = paste0(results_dir, "/All_UMAP_FLI1.pdf"), width = 7, height = 6)
print(FeaturePlot(sc_obj, features = "FLI1", order = TRUE, ncol = 1, raster = TRUE) &
        scale_colour_gradientn(colours = brewer.pal(n = 9, "YlGnBu")[3:9]))
dev.off()

# Test FLI1 between conditions within endothelial cells (the only cell type with
# appreciable FLI1 expression). NOTE: the original script annotated this plot
# with a p-value mistakenly taken from an unrelated C19orf33 test; corrected here
# to test FLI1 in endothelial cells (SSc vs HC). Verify the comparison matches
# the published "p = 0.75".
en_set <- subset(sc_obj, subset = types == "En")
Idents(en_set) <- "condition"
fli1_en <- FindMarkers(en_set, ident.1 = "SSc", ident.2 = "HC",
                       features = "FLI1", min.pct = 0, logfc.threshold = 0)

pdf(file = paste0(results_dir, "/All_Vln_FLI1.pdf"), width = 6, height = 3)
print(VlnPlot(sc_obj, features = "FLI1", group.by = "types", split.by = "condition",
              pt.size = 0, cols = cond_cols) &
        theme(axis.title.x = element_blank(), legend.position = "none") &
        annotate(geom = "text", label = paste0("p=", signif(fli1_en$p_val[1], 2)),
                 x = Inf, y = Inf, hjust = 1, vjust = 1, colour = "darkred"))
dev.off()


####  CD45+ IMMUNE-CELL RE-INTEGRATION AND SUBTYPING  #############################################
## Subset the immune compartment (lymphoid, myeloid, mast), re-integrate it on its own, and resolve
## T-cell (CD4+/CD8+) and myeloid subtypes. The per-cell T-cell-subtype labels feed the cell-cell
## communication analysis (09_cellchat_communication.R); the myeloid subtypes provide a finer
## immune annotation. NOTE: the myeloid subtyping requires an external reference atlas
## (myeloid_ref_rds) and is not used by the main figures.

immune_set <- subset(sc_obj, types %in% c("L", "My", "Mc"))

## Re-integrate the immune cells (per-sample SCTransform v2, SCT/CCA anchoring).
immune_list <- SplitObject(immune_set, split.by = "orig.ident")
immune_list <- lapply(immune_list, function(x) SCTransform(x, vst.flavor = "v2", variable.features.n = 2000, verbose = FALSE))
int_features <- SelectIntegrationFeatures(immune_list, nfeatures = 2000)
immune_list  <- PrepSCTIntegration(immune_list, anchor.features = int_features)
anchors <- FindIntegrationAnchors(object.list = immune_list, normalization.method = "SCT", dims = 1:30,
                                  anchor.features = int_features, k.filter = 10, k.anchor = 10, k.score = 10)
immune_set <- IntegrateData(anchorset = anchors, normalization.method = "SCT", k.weight = 50)

immune_set <- RunPCA(immune_set, assay = "integrated")
immune_set <- RunUMAP(immune_set, assay = "integrated", dims = 1:20, n.neighbors = 25L,
                      min.dist = 0.25, n.epochs = 1000)
immune_set <- FindNeighbors(immune_set, reduction = "pca", dims = 1:15, k.param = 27, n.trees = 1000)
immune_set <- FindClusters(immune_set, graph.name = "integrated_snn", resolution = 0.5)

pdf(file = paste0(results_dir, "/Immune_UMAP_celltype.pdf"), width = 6, height = 3)
print(DimPlot(immune_set, reduction = "umap", group.by = "types", shuffle = TRUE, raster = FALSE,
              cols = c("#7AC5CD", "#4682B4", "#EEA2AD")) + ggtitle("") &
        theme(axis.text = element_blank(), axis.title = element_text(size = 9), axis.ticks = element_blank()))
dev.off()

save(immune_set, file = paste0(results_dir, "/eso_CD45_integratedObj.RData"))

## --- T cells: resolve CD4+ / CD8+ subsets --------------------------------------
## Re-cluster the lymphoid cells (ribosomal genes dropped from the variable features) and label the
## low-resolution clusters by CD8A / CD4 expression.
t_set <- subset(immune_set, types == "L")
DefaultAssay(t_set) <- "integrated"
VariableFeatures(t_set) <- grep("^RP[SL]|^MRPS|^MRPL", VariableFeatures(t_set), value = TRUE, invert = TRUE)
t_set <- RunPCA(t_set)
t_set <- RunUMAP(t_set, dims = 1:10, reduction = "pca")
t_set <- FindNeighbors(t_set, dims = 1:10, reduction = "pca")
t_set <- FindClusters(t_set, resolution = 0.2)
## Clusters 0 / 1 / 2 correspond to CD8+ T, CD4+ T, and other lymphoid (confirmed by CD8A / CD4).
t_set$t_cell_subtype <- mapvalues(as.character(t_set$seurat_clusters),
                                  from = c("0", "1", "2"),
                                  to   = c("CD8+ T", "CD4+ T", "Other Lymphoid"))

## --- Myeloid: subtype by label transfer from an external myeloid reference ------
my_set <- subset(immune_set, types == "My")
myeloid_ref <- readRDS(myeloid_ref_rds)
my_pred <- SingleR(test   = as.SingleCellExperiment(my_set),
                   ref    = GetAssayData(myeloid_ref, assay = "RNA", layer = "data"),
                   labels = myeloid_ref$Clusters2, num.threads = 2)
my_set$myeloid_subtype <- my_pred$pruned.labels

## --- Write the resolved immune subtypes back onto the whole-tissue object -------
sc_obj$t_cell_subtype <- NA
sc_obj$t_cell_subtype[colnames(t_set)]  <- t_set$t_cell_subtype
sc_obj$myeloid_subtype <- NA
sc_obj$myeloid_subtype[colnames(my_set)] <- my_set$myeloid_subtype
save(sc_obj, file = paste0(results_dir, "/All_integratedObj.RData"))
