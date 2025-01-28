####  LOAD LIBRARIES  ##############################################################################

.libPaths('~/local/R_libs/')    # Set path using .libPaths function

library(methods)
library(Seurat)
library(sctransform)
library(ggplot2)
library(glmGamPoi)
library(future)
library(cowplot)
#library(scater)
library(scales)
library(future)
library(future.apply)

####  SET PARAMETERS  ##############################################################################

set.seed(0123456789)
plan('multisession', workers=6)

#Memory Options
Sys.setenv('R_MAX_VSIZE'=760000000000)
options(future.globals.maxSize = 760000 * 1024^2) # NOTE Calculated for 640 Gb RAM request
maxMem <- Sys.getenv('R_MAX_VSIZE')
#print(maxMem)

scl_weak <- c('P0449-P','P0449-D','P0483-P','P0483-D','P0523-P',
              'P0523-D','P0542-P','P0542-D','P0541-P','P0541-D')

scl_absent <- c('P0491-P','P0491-D','P0628-P','P0628-D','P0630-P','P0630-D','P0656-P','P0656-D')

scl_samples <- c(scl_weak, scl_absent, 'P0080-P', 'P0080-D')

gerd <- c('IW0065-P','IW0065-D','IW0071-D','IW0071-P','IW0076-D','IW0076-P','IW0077-D','IW0077-P')

hcs <- c('HCE40-P','HCE40-D','HCE43-D', 'HCE43-P','HCE047-D','HCE048-P',
         'HCE048-D','HCE049-P','HCE049-D','HCE047-P','HCE051-D','HCE051-P')

#ref_samples <- hcs
ref_samples <- c('HCE051-P','HCE051-D','HCE047-P','HCE047-D', 'HCE048-P', 'HCE048-D')
#ref_samples <- c('HCE048-D','HCE051-P','IW0065-D','IW0076-P', 'P0491-D', 'P0523-P')
#ref_samples <- c('HCE047-P','HCE051-D','IW0065-D','IW0076-P', 'P0523-D', 'P0541-P')
#ref_samples <- c('HCE048-D','HCE051-D','HCE051-P','IW0065-D','IW0065-P',
#	         'IW0076-P','P0491-D','P0628-D','P0628-P')

#ref_samples <- c('HCE047-D','HCE048-D','HCE048-P','HCE051-P','IW0065-D','IW0076-P', 'P0491-D', 'P0523-P')

####  FUNCTIONS  ###################################################################################


#Load Options
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

####  MAIN  #######################################################################################

# Set root directory where sample sub-directories are located
root_dir <- '/projects/p31763/users/mdapas/MP_SSc_scRNAseq/Matrices/'
output_dir <- '/projects/p31763/users/mdapas/MP_SSc_scRNAseq'

condition='condition'; location='location';

# Get sample names according to directories within root directory
samples <- list.dirs(root_dir, recursive=F, full.names=F)
samples <- samples[!samples %in% c('P0080-P')]

scrna.list <- list();

for (i in 1:length(samples)) {
  sample <- samples[i]
  print(sample)
  
  scrna.list[[i]] <- loadRData(paste0(root_dir, sample, '/analysis/sample_QC/seuratObj_', sample, '.RData'))
    
  # Add case/control, sample location metadata
  if (substr(sample,1,1)=='H') {
    scrna.list[[i]][[condition]] <- 'HC'
  } else if (substr(sample,1,1)=='P') {
    scrna.list[[i]][[condition]] <- 'SSc'
  } else {
    scrna.list[[i]][[condition]] <- 'GERD'
  }
  if (substr(sample,nchar(sample),nchar(sample))=='P') {
    scrna.list[[i]][[location]] <- 'Proximal'
  } else {
    scrna.list[[i]][[location]] <- 'Distal'
  }
}
names(scrna.list) <- samples

ref_list <- which(names(scrna.list) %in% ref_samples)

scrna.list <- future_lapply(X=scrna.list, FUN = SCTransform, vst.flavor="v2")

int_features <- SelectIntegrationFeatures(scrna.list, nfeatures = 3000)
save(int_features, file=paste0(output_dir,'eso_intRefFeatures_2.RData'))
#int_features <- loadRData(paste0(output_dir,'eso_intFeatures.RData'))

scrna.list <- PrepSCTIntegration(scrna.list, anchor.features = int_features)

scrna.list <- future_lapply(X=scrna.list, FUN = RunPCA, features=int_features)

anchors <- FindIntegrationAnchors(object.list=scrna.list, normalization.method='SCT',
                                  anchor.features=int_features, reference=ref_list)
save(anchors, file=paste0(output_dir,'/eso_intRefAnchors_2.RData'))
#ref_anchors <- loadRData(paste0(output_dir,'/eso_intAnchors.RData'))

plan('multisession', workers=1)
integrated_set <- IntegrateData(anchorset=anchors, normalization.method='SCT')

integrated_set <- RunPCA(object=integrated_set, assay='integrated')
DefaultAssay(integrated_set) <- 'integrated'

save(integrated_set, file=paste0(output_dir,'/eso_integratedRefObj_2.RData'))
