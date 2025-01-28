#!/bin/bash
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 35:59:59
#SBATCH --mem=540G
#SBATCH --job-name="epi_reint_ccDiffRegress"
#SBATCH --output %x.%j.out 
#SBATCH --partition=genhimem
#SBATCH --account=p31763
#SBATCH --mail-user=mdapas@northwestern.edu
#SBATCH --mail-type=END
########################

module purge all
module load R/4.2.0
module load geos/3.8.1

Rscript --vanilla /projects/p31763/users/mdapas/MP_SSc_scRNAseq/shell/re_integrate_refBased.R
