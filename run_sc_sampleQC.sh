#!/bin/bash
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

shell_dir=/projects/p31763/users/mdapas/MP_SSc_scRNAseq/shell
sample_dir=/projects/p31763/users/mdapas/MP_SSc_scRNAseq/Matrices

Rscript --vanilla $shell_dir/sc_sampleQC.R -d $sample_dir
