####  LOAD LIBRARIES  ##############################################################################

.libPaths('~/local/R_libs/')    # Set R library path

suppressPackageStartupMessages({
  library(methods) # 
  library(Seurat)
  library(ggplot2)
  library(scater)
  library(scales)
  library(scDblFinder) # 
  library(cowplot)
  library(SeuratWrappers)
  library(flexmix)
  library(optparse)
  library(RColorBrewer)
})


####  SET PARAMETERS  ##############################################################################

options(scipen=10000)
set.seed(123456789)


####  FUNCTIONS  ###################################################################################

get_qc <- function(sc_obj) {
  ### Creates vector indicating any QC metrics each cell failed on
  hi_count <- 100000
  miqc_cut <- 0.75
  mt_cut <- 5
  qc_vector <- rep('Pass', ncol(sc_obj))
  qc_vector <- ifelse(sc_obj@meta.data[[count]] > hi_count, 'High_nCount', qc_vector)
  qc_vector <- ifelse(sc_obj@meta.data[['miQC.keep']]=='discard',
                      ifelse(qc_vector=='Pass', 
                             ifelse(sc_obj@meta.data[[mt]]>mt_cut, 'miQC', qc_vector),
                                    paste0(qc_vector, ', miQC')), qc_vector)
  qc_vector[is.na(qc_vector)] <- 'Doublet'
return(qc_vector)
}


get_qc_plots <- function(sc_obj) {
  ### Returns 2x2 grid of QC plots for given Seurat object
  pass_total <- sum(sc_obj[[qc]]=='Pass')
  pass_percent <- round(pass_total/dim(sc_obj)[2],3)*100
  pass_label <- paste0('\nPass\n(',pass_total,', ',pass_percent,'%)\n')
  fail_total <- sum(sc_obj[[qc]]!='Pass')
  fail_percent <- round(fail_total/dim(sc_obj)[2],3)*100
  fail_label <- paste0('Fail\n(',fail_total,', ',fail_percent,'%)')
  sc_obj[['pass_n']] <- ifelse(sc_obj[[qc]]=='Pass',pass_label, fail_label)
  sc_obj@meta.data$pass_n <- factor(sc_obj@meta.data$pass_n, levels=c(pass_label,fail_label))
  sc_obj[['pass']] <- ifelse(sc_obj[[qc]]=='Pass','Pass', 'Fail')
  sc_obj@meta.data$pass <- factor(sc_obj@meta.data$pass, levels=c('Pass','Fail'))
  
  n_singlets <- sum(sc_obj@meta.data$scDblFinder.class=='singlet')
  n_doublets <- sum(sc_obj@meta.data$scDblFinder.class=='doublet')
  singlet_label <- paste0('\nSinglets\n(',n_singlets,', ',round(n_singlets/dim(sc_obj)[2],3)*100,'%)\n')
  doublet_label <- paste0('Doublets\n(',n_doublets,', ',round(n_doublets/dim(sc_obj)[2],3)*100,'%)')
  sc_obj[['doublets']] <- ifelse(sc_obj@meta.data$scDblFinder.class=='singlet', singlet_label, doublet_label)
  sc_obj@meta.data$doublets <- factor(sc_obj@meta.data$doublets, levels=c(singlet_label,doublet_label))
  sc_obj[['scDblFinder.class']] <- ifelse(sc_obj@meta.data$scDblFinder.class=='singlet', 'Singlet', 'Doublet')
  sc_obj@meta.data$scDblFinder.class <- factor(sc_obj@meta.data$scDblFinder.class, levels=c('Singlet', 'Doublet'))
  
  # Pre-QC clusters and doublet identification
  sc_obj <- SCTransform(sc_obj, verbose=F)
  sc_obj <- RunPCA(sc_obj, npcs=20, verbose=F)
  sc_obj <- RunUMAP(sc_obj, dims = 1:20, verbose=F)
  
  return(list(
      plot_grid(FeatureScatter(sc_obj, feature1=count, feature2=features, group.by='pass',
                                  plot.cor=F, pt.size=0.1, cols=c('#F8766D','gray'))
            + theme(legend.position='none', axis.title=element_text(size=10), axis.text=element_text(size=9),
                    title=element_text(size=11, face='bold')),
            FeatureScatter(sc_obj, feature1=features, feature2=mt, group.by='pass',
                           plot.cor=F, pt.size=0.1, cols=c('#F8766D','gray'))
            + theme(legend.position=c(.75,.9), axis.title=element_text(size=10), axis.text=element_text(size=9),
                    title=element_text(size=11, face='bold'), legend.title=element_blank(), 
                    legend.text=element_text(size=8)) + guides(color=guide_legend(override.aes=list(size=2))),
            FeatureScatter(sc_obj, feature1=rb, feature2=mt, group.by='pass',
                           plot.cor=F, pt.size=0.1, cols=c('#F8766D','gray'))
            + theme(legend.position='none', axis.title=element_text(size=10), axis.text=element_text(size=9),
                    title=element_text(size=10, face='bold')),
            FeatureScatter(sc_obj, feature1=features, feature2=rb, group.by='pass',
                           plot.cor=F, pt.size=0.1, cols=c('#F8766D','gray'))
            + theme(legend.position='none', axis.title=element_text(size=10), axis.text=element_text(size=9),
                    title=element_text(size=11, face='bold'))),
    plot_grid(DimPlot(sc_obj,reduction='umap',group.by = 'doublets', cols=c('#F8766D','turquoise1')) + 
                theme(plot.title=element_blank()),
              VlnPlot(sc_obj, features = "nFeature_RNA", group.by = 'scDblFinder.class', 
                      pt.size = 0.1, cols=c('#F8766D','turquoise1')) + 
                theme(legend.position='none', plot.title=element_blank(), axis.title.x = element_blank()) + 
                ylab('N Features'), ncol=1),
    DimPlot(sc_obj, reduction = "umap", group.by='pass_n', cols=c('#F8766D','gray')) + 
      ggtitle(' ')))
}

get_subset <- function(sc_obj, pos_markers, neg_markers) {
  ### Identifies epithelial cells based on expression of provided positive and negative markers
  epi_clusts <- c()
  
  marker_expression <- AverageExpression(sc_obj, features=c(epi_pos_markers, epi_neg_markers), 
                                         group.by='seurat_clusters', slot='count', assays='RNA')[[1]]
  clusters <- seq(0,ncol(marker_expression)-1,1)
  
  marker_sums <- data.frame(rbind(colSums(marker_expression[epi_pos_markers,]),
                                  colSums(marker_expression[epi_neg_markers,])),
                            row.names = c('Pos','Neg')) 
  names(marker_sums) <- clusters
  
  print('     Epithelial marker expression by cluster:')
  print(marker_sums)

  for (marker in epi_pos_markers) {
    epi_clusts <- unique(c(epi_clusts, clusters[marker_expression[marker,]>10]))
  }

  for (marker in epi_neg_markers) {
    epi_clusts <- epi_clusts[!epi_clusts %in% clusters[marker_expression[marker,]>2]]
  }
  return(epi_clusts)
}


####  OPTIONS  #####################################################################################

option_list = list(
  make_option(c("-d", "--dir"), action='store', type="character", default=NULL, 
              help="root directory")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


####  MAIN  #######################################################################################


# Get root directory where sample sub-directories are located
root_dir <- opt$dir

# Variable declarations
count='nCount_RNA'; features='nFeature_RNA'; mt='percent.mt'; rb='percent.rb'; qc='qc'

# Get sample names according to directories within root directory
samples <- list.dirs(root_dir, recursive=F, full.names=F)

# Loop through samples, create Seurat objects, performing QC
cell_qc_df <- data.frame(Filter=character())
sample_preQC_df <- setNames(data.frame(matrix(ncol = 14, nrow = 0)), 
                            c('sample', 'cells', 
                              paste0('med_', count), paste0('min_', count), paste0('max_', count), 
                              paste0('med_', features), paste0('min_', features), paste0('max_', features),
                              paste0('med_', mt), paste0('min_', mt), paste0('max_', mt), 
                              paste0('med_', rb), paste0('min_', rb), paste0('max_', rb)))
sample_postQC_df <- sample_preQC_df

sample <- c('HCE43-D')

for (i in 1:length(samples)) {
  sample <- samples[i]
  print(sample)
  
  # Create output directory (if necessary)
  output_dir <- paste(root_dir, sample, 'analysis/sample_QC', sep='/')
  if (!file.exists(output_dir)) {
    dir.create(output_dir, recursive=T)
  }
  
  data.10x <- Read10X(data.dir=paste(root_dir,sample, 'filtered_feature_bc_matrix',sep='/'))
  
  print('  creating Seurat object...')
  sc_obj = CreateSeuratObject(counts=data.10x, project=sample, min.features=200)
  rm(data.10x)
  sc_obj[[mt]] <- PercentageFeatureSet(sc_obj, pattern = "^MT-")
  sc_obj[[rb]] <- PercentageFeatureSet(sc_obj, pattern = "^RP[SL]")
  
  print('  performing QC...')
  
  # Predict doublets (can adjust dbr.sd to lower prediction threshold)
  sc_obj <- as.Seurat(scDblFinder(as.SingleCellExperiment(sc_obj), dbr.sd=0.01))
  sc_obj@meta.data[,c(6,8,9,10)] <- NULL
  
  # Apply MiQC to get sample-specific MT cutoffs, run on non-doublet set
  sc_obj_nonDbl <- subset(sc_obj, subset=scDblFinder.class=='singlet')
  sc_obj_nonDbl <- RunMiQC(sc_obj_nonDbl, percent.mt=mt, nFeature_RNA=features, 
                           posterior.cutoff = 0.75, model.slot = "flexmix_model",
                           backup.option='percent', backup.percent=10)
  
  sc_obj@meta.data[rownames(sc_obj@meta.data) %in% rownames(sc_obj_nonDbl@meta.data),
                   c('miQC.probability','miQC.keep')] <- sc_obj_nonDbl@meta.data[,c('miQC.probability','miQC.keep')]

  # Print miQC probabilities
  if (sum(sc_obj_nonDbl@meta.data$miQC.probability)>0) {
  pdf(file=paste0(output_dir, '/miQC_', sample, '.pdf'), width=6, height=5)
  myPalette <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))
  sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0,1))
  print(PlotMiQC(sc_obj_nonDbl, color.by = "miQC.probability") + sc + 
    theme(axis.title=element_text(size=10), axis.text=element_text(size=9),
            title=element_text(size=10, face='bold'), legend.text=element_text(size=8),
          legend.title=element_text(size=9), plot.title=element_text(hjust = 0.5)) + ggtitle(sample))
  dev.off()
  }
  
  rm(sc_obj_nonDbl)
  
  sample_preQC_df[i,] <- c(sample, ncol(sc_obj), 
                           median(sc_obj[[count]][,1]), min(sc_obj[[count]][,1]), max(sc_obj[[count]][,1]),
                           median(sc_obj[[features]][,1]), min(sc_obj[[features]][,1]), max(sc_obj[[features]][,1]),
                           median(sc_obj[[mt]][,1]), min(sc_obj[[mt]][,1]), max(sc_obj[[mt]][,1]),
                           median(sc_obj[[rb]][,1]), min(sc_obj[[rb]][,1]), max(sc_obj[[rb]][,1]))
  
  # Determine QC pass/fail
  sc_obj[[qc]] <- get_qc(sc_obj) 
  
  temp_df <- data.frame(table(sc_obj[[qc]]))
  names(temp_df) <- c('Filter',sample)
  cell_qc_df <- merge(cell_qc_df, temp_df, by='Filter', all=T)
  
  # QC scatter plots
  qc_plots <- get_qc_plots(sc_obj)
  
  pdf(file=paste0(output_dir, '/QC_', sample, '.pdf'), width=12, height=5.5)
  print(plot_grid(ggdraw() + draw_label(paste(sample, '- Initial QC'), fontface='bold'), 
                  plot_grid(qc_plots[[1]],qc_plots[[3]], ncol=2), 
                  ncol=1, rel_heights=c(0.1, 1), scale=0.98))
  dev.off()
  
  pdf(file=paste0(output_dir, '/Doublets_', sample, '.pdf'), width=7, height=7)
  print(plot_grid(ggdraw() + draw_label(paste(sample, '- Predicted Doublets'), fontface='bold'), 
                  qc_plots[[2]], ncol=1, rel_heights=c(0.1, 1), scale=0.98))
  dev.off()
  
  # Apply QC
  sc_obj <- subset(sc_obj, subset=qc=='Pass')
  
  sample_postQC_df[i,] <- c(sample, ncol(sc_obj), 
                            median(sc_obj[[count]][,1]), min(sc_obj[[count]][,1]), max(sc_obj[[count]][,1]),
                            median(sc_obj[[features]][,1]), min(sc_obj[[features]][,1]), max(sc_obj[[features]][,1]),
                            median(sc_obj[[mt]][,1]), min(sc_obj[[mt]][,1]), max(sc_obj[[mt]][,1]),
                            median(sc_obj[[rb]][,1]), min(sc_obj[[rb]][,1]), max(sc_obj[[rb]][,1]))

  # Save individual sample data
  save(sc_obj, file=paste0(output_dir,'/seuratObj_', sample, '.RData'))
}

cell_qc_df <- setNames(data.frame(t(cell_qc_df[ , - 1])), cell_qc_df[ , 1])

# Save QC data
write.table(cell_qc_df, paste(root_dir, 'QC_summary_all.tsv',sep='/'), quote=F, row.names=T,col.names=T, sep="\t")
write.table(sample_preQC_df, paste(root_dir, 'raw_metrics_all.tsv',sep='/'), quote=F, row.names=F,col.names=T, sep="\t")
write.table(sample_postQC_df, paste(root_dir, 'qc_metrics_all.tsv',sep='/'), quote=F, row.names=F,col.names=T, sep="\t")
