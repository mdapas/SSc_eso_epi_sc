#!/bin/bash
# ==============================================================================
# SLURM submission wrapper for 01_sample_qc.R
#
# Part of the analysis for:
#   Dapas et al. (2025) "Cellular and molecular dysregulation of the esophageal
#   epithelium in systemic sclerosis."
#
# Runs per-sample single-cell quality control over a directory of Cell Ranger
# outputs. Edit the SBATCH directives and the two path variables below for your
# cluster/layout, then submit with: sbatch run_sc_sampleQC.sh
# ==============================================================================
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 9:59:59
#SBATCH --mem=50G
#SBATCH --job-name="sc_sampleQC"
#SBATCH --output ./logs/%x.%j.out
#SBATCH -p genomics
#SBATCH -A b1042
########################

module purge all
module load R/4.2.0
module load geos/3.8.1

# --- Paths (edit to match your environment) -----------------------------------
# shell_dir:  directory containing this script and 01_sample_qc.R
# sample_dir: directory containing per-sample Cell Ranger output subdirectories
shell_dir="./shell"
sample_dir="./data/matrices"

Rscript --vanilla "${shell_dir}/01_sample_qc.R" -d "${sample_dir}"
