#!/bin/bash
# ==============================================================================
# SLURM submission wrapper for 02_integrate_allcells.R
#
# Part of the analysis for:
#   Dapas et al. (2025) "Cellular and molecular dysregulation of the esophageal
#   epithelium in systemic sclerosis."
#
# Reference-based SCTransform/CCA integration of all QC-passing samples.
# Requires a high-memory node. Edit the SBATCH directives and shell_dir for your
# cluster/layout, then submit with: sbatch run_integration.sh
# ==============================================================================
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 35:59:59
#SBATCH --mem=540G
#SBATCH --job-name="eso_integration"
#SBATCH --output %x.%j.out
#SBATCH --partition=genhimem
#SBATCH --account=p31763
########################

module purge all
module load R/4.2.0
module load geos/3.8.1

# shell_dir: directory containing this script and 02_integrate_allcells.R
shell_dir="./shell"

Rscript --vanilla "${shell_dir}/02_integrate_allcells.R"
