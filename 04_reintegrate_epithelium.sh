#!/bin/bash
# ==============================================================================
# SLURM submission wrapper for 04_reintegrate_epithelium.R
#
# Part of the analysis for:
#   Dapas et al. (2025) "Cellular and molecular dysregulation of the esophageal
#   epithelium in systemic sclerosis."
#
# Subsets the epithelial cells from the annotated all-cell object and re-runs
# reference-based SCTransform/CCA integration (regressing the cell-cycle
# difference) on that subset. Requires a high-memory node.
# Edit the SBATCH directives and shell_dir for your cluster/layout, then submit
# with: sbatch run_reintegration.sh
# ==============================================================================
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 35:59:59
#SBATCH --mem=540G
#SBATCH --job-name="eso_epi_reintegration"
#SBATCH --output %x.%j.out
#SBATCH --partition=genhimem
#SBATCH --account=p31763
########################

module purge all
module load R/4.2.0
module load geos/3.8.1

# shell_dir: directory containing this script and 04_reintegrate_epithelium.R
shell_dir="./shell"

Rscript --vanilla "${shell_dir}/04_reintegrate_epithelium.R"
