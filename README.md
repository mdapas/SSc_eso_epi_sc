# Esophageal epithelium in systemic sclerosis — analysis code

Analysis code for:

> Dapas et al. (2025). *Cellular and molecular dysregulation of the esophageal epithelium in systemic sclerosis.* JCI Insight.

The repository reproduces the single-cell (scRNA-seq), spatial (CosMx), and clinical-association
analyses reported in the paper, starting from per-sample count matrices and ending in the figure
and supplemental-table source data.

---

## Pipeline

The scripts are numbered in run order. Each consumes the saved outputs of earlier steps, so they are
intended to be run top to bottom. Steps 01–02 and 04 are wrapped in SLURM submission scripts
(`*.sh`) because they are memory-intensive; the rest run with `Rscript <script>.R` from the
repository root.

| # | Script | Purpose |
|---|--------|---------|
| 01 | `01_sample_qc.R` | Per-sample single-cell QC (miQC mitochondrial cutoff, doublet removal) |
| 02 | `02_integrate_allcells.R` | Reference-based SCTransform/CCA integration of all samples |
| 03 | `03_annotate_celltypes.R` | Whole-tissue cell typing + composition; FLI1 expression; CD45+ immune re-integration and T-cell / myeloid subtyping |
| 04 | `04_reintegrate_epithelium.R` | Subset epithelial cells and re-integrate (cell-cycle difference regressed) |
| 05 | `05_epithelial_compartments.R` | Epithelial differentiation compartments, proportions, per-compartment DEGs |
| 06 | `06_pseudobulk_dgea.R` | edgeR pseudobulk differential expression + TF over-representation |
| 07 | `07_gsea_enrichment.R` | fgsea gene-set enrichment of the pseudobulk results |
| 08 | `08_superficial_analysis.R` | Superficial-compartment subclustering, module scores, TF targets, proportions |
| 09 | `09_cellchat_communication.R` | CellChat cell–cell communication (3-group and SSc-vs-rest) |
| 10 | `10_cosmx_spatial.R` | CosMx spatial annotation (label transfer) and cell-cell proximity |
| 11 | `11_clinical_associations.R` | Ordinal regression / PCA linking molecular features to clinical phenotypes |

### Run order and dependencies

```
01_sample_qc ─► 02_integrate_allcells ─► 03_annotate_celltypes ─┬─► 04_reintegrate_epithelium ─► 05_epithelial_compartments ─┬─► 06_pseudobulk_dgea ─► 07_gsea_enrichment
                                                                │                                                            ├─► 08_superficial_analysis
                                                                │                                                            │
                                                                └────────────────────────────────────────────────────────► 09_cellchat_communication
                                                                                                                             │
  10_cosmx_spatial  (uses the annotated references from 03 / 05 / 08)                                                        │
                                                                                                                             └─► 11_clinical_associations
                                                                                                                                 (uses 06 / 08 outputs)
```

---

## Per-script inputs and key outputs

Paths are relative to the repository root; adjust `data_dir` / `results_dir` at the top of each script.

**01_sample_qc.R** — reads Cell Ranger `filtered_feature_bc_matrix` per sample; writes a QC-passing
Seurat object and QC diagnostics per sample, plus a pass/fail summary table.

**02_integrate_allcells.R** — reads the per-sample QC objects; writes `eso_intRefFeatures.RData`,
`eso_intRefAnchors.RData`, `eso_integratedRefObj.RData`.

**03_annotate_celltypes.R** — reads `eso_integratedRefObj.RData` (and an external myeloid reference,
`myeloid_ref_rds`); writes `All_integratedObj.RData` (annotated whole-tissue object, with
`t_cell_subtype` / `myeloid_subtype` immune labels), `eso_CD45_integratedObj.RData`,
`eso_ct_markers.tsv`, `eso_cellType_props.tsv`.

**04_reintegrate_epithelium.R** — reads `All_integratedObj.RData`; writes `eso_epi_integratedObj.RData`.

**05_epithelial_compartments.R** — reads `eso_epi_integratedObj.RData`; writes per-compartment DEG
objects (`epi_degs_*`), `epi_layerRep_props.tsv`, `epi_layerRep_props_stats.tsv`,
`epi_layer_props.tsv`, `epi_cluster_props.tsv`, `diffScoreStats.tsv`, and the re-saved epithelial
object with compartment labels.

**06_pseudobulk_dgea.R** — reads `eso_epi_integratedObj.RData`; writes
`pb_deg_res_layerRep_byRegion.RData`, `pbObj_layerRep_byRegion.RData`, the combined DEG workbook,
GO/KEGG/TF over-representation tables, `proximal_distal_correlation.tsv`,
`compartment_logfc_jitter.tsv`.

**07_gsea_enrichment.R** — reads the by-location DEG lists from 06; writes `GSEAtable.tsv`.

**08_superficial_analysis.R** — reads `eso_epi_integratedObj.RData`; writes `eso_epi_superficial.RData`,
`epiSuper_cluster_props.tsv`, `epiSuper_clusters_props_stats.tsv`, and superficial figure source data
(`superficial_umap_coords.tsv`, `superficial_marker_expression.tsv`, `superficial_module_score.tsv`,
`metallothionein_umap.tsv`, `metallothionein_by_subcluster.tsv`, `superficial_resolution_scan.tsv`,
`superficial_marker_module_umap.tsv`, `tf_<tf>_targets.tsv`, `cluster5_immune_genes.tsv`,
`fli1_target_module.tsv`).

**09_cellchat_communication.R** — reads `All_integratedObj.RData` (incl. `t_cell_subtype`) and
`eso_epi_integratedObj.RData`; writes the per-condition and merged CellChat objects and interaction
tables.

**10_cosmx_spatial.R** — reads the raw CosMx per-sample `.RDS` objects and the scRNA-seq reference
objects (03 / 05 / 08); writes `Eso_CosMx_seuratObj.RData`, the SingleR label objects
(`cosMx_labels.RData`, `cosMx_epiGepiLabels.RData`, `cosMx_epiSuperLabels.RData`),
`cosMx_spatialProximityAnalysis.RData`, and spatial / proximity source data
(`adjacency_df_epi.tsv`, `all_adj_cond.tsv`, `fibroblast_neighbor_composition.tsv`,
`cosmx_celltype_counts.tsv`).

**11_clinical_associations.R** — reads the pseudobulk objects (06), the superficial object and
proportions (08), the compartment summaries (05), and the clinical tables; writes
`clinical_trait_correlations.tsv`, `compartment_props_clinical_pcs.tsv`,
`clinical_TFgene_ordinalReg.tsv`, and a PCA biplot table.

---

## Figure / table coverage

The scripts use descriptive section titles rather than manuscript figure numbers. This table maps
each manuscript panel to the script that produces its source data.

| Manuscript item | Script(s) |
|-----------------|-----------|
| Figure 1; Supplemental Figure 1 | 03 |
| Supplemental Figure 9 A–B (FLI1) | 03 |
| Supplemental Table 4 (cell-type markers) | 03 |
| Figure 2 A/B/C/E; Supplemental Figure 3 | 05 |
| Figure 2 D (spatial) | 10 |
| Figure 3 A/C/D; Supplemental Figure 5 | 06 |
| Figure 3 B (GSEA) | 07 |
| Supplemental Tables 5–12 (DEG / TF tables) | 06 |
| Figure 4 A/B/C/E/F/G/I; Supplemental Figure 6; Supplemental Figure 9 C | 08 |
| Figure 4 D/H (spatial) | 10 |
| Figure 5 A/B/C; Supplemental Figure 8 A–C | 09 |
| Figure 5 D/E/F; Supplemental Figure 8 D | 10 |
| Supplemental Figure 4 (CosMx maps + counts) | 10 |
| Figure 6; Supplemental Figure 10 | 11 |
| Supplemental Table 2 (sample QC / exclusion) | 01 / 02 |

---

## Data and conventions

**scRNA-seq cohort.** SSc (n = 10), GERD (n = 4), HC (n = 6) — 20 participants, each with a proximal
and a distal biopsy (40 samples). One sample (`SSc10-P`) is excluded for low quality, leaving 39.

**CosMx cohort.** Four samples (`SSc_11`, `SSc_12`, `HC_7`, `HC_8`), each with proximal and distal
regions. The `HC_7` proximal/distal sections were mislabeled on the slide and are swapped in code.

**Sample IDs.** Samples are read in under their de-identified IDs. Condition is taken from the ID
prefix (`HC` / `GERD` / `SSc`) and biopsy location from the `-P` (proximal) / `-D` (distal) suffix.
Original participant identifiers are not stored anywhere in this repository.

**Integration reference.** Reference-based integration anchors on three healthy-control proximal/distal
pairs: `HC3`, `HC6`, `HC5`.

**External inputs not included here.**
- An annotated myeloid reference atlas for myeloid subtyping in 03 (`myeloid_ref_rds`).
- Pathway gene-set maps for the over-representation analysis in 06.
- MSigDB C2 gene sets for the GSEA in 07.
- The raw CosMx per-sample Seurat `.RDS` objects for 10.
- Clinical tables (`clinical_data.txt`, `ssc_clin_dat.txt`) for 11, keyed on the de-identified IDs.

---

## Requirements

R (developed under 4.2.0). Core packages: Seurat (v4/SCTransform v2), sctransform, glmGamPoi, harmony,
Banksy, SingleR, SingleCellExperiment, scDblFinder, miQC, edgeR, limma, fgsea, clusterProfiler,
CellChat, rms, FactoMineR, propeller/speckle, and the usual tidyverse / plotting stack
(ggplot2, dplyr, tidyr, pheatmap, cowplot, RColorBrewer). See the `LOAD LIBRARIES` block at the top of
each script for its specific dependencies.
