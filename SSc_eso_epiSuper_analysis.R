####  LOAD LIBRARIES  ##############################################################################

.libPaths('~/local/R_libs/')    # Set path using .libPaths function

#Load Libraries
dyn.load('/home/mdn578/local/libs/usr/lib64/libfftw3.so.3')

#BiocManager::install("glmGamPoi")
library(methods)
#.libPaths("/home/mhc0155/R/x86_64-pc-linux-gnu-library/4.1")
library(metap)
library(Seurat)
library(SeuratWrappers)
library(monocle)
library(sctransform)
library(ggplot2)
library(ggpubr)
library(glmGamPoi)
library(cowplot)
library(plyr)
library(RColorBrewer)
library(future)
library(future.apply)
library(scCustomize)
library(monocle3)
library(R.utils)
library(lme4)
library(MASC)
library(speckle)
library(limma)
library(eulerr)
library(Lamian)


####  SET PARAMETERS  ##############################################################################

set.seed(0123456789)
plan('multisession', workers=6)

##Memory Options
Sys.setenv('R_MAX_VSIZE'=192000000000)
options(future.globals.maxSize = 192000 * 1024^2) # NOTE Calculated for 160 Gb RAM request
maxMem <- Sys.getenv('R_MAX_VSIZE')

cond='condition'; loc='location'; 

cond_cols <- c('#8156B3','#e26b53','#A8A39d')
layer_cols_3 <- c('#4ca5b1','#e9b85d','#c72f4c')
layer_cols_5 <- brewer.pal(n=11,'Spectral')[c(10,9,7,4,2)]

#layer_cols_5[1] <- layer_cols_3[1]
layer_cols_5[4:5] <- layer_cols_3[2:3]

root_dir <- '/projects/p31763/users/mdapas/MP_SSc_scRNAseq/Matrices/'

output_dir <- '/projects/p31763/users/mdapas/MP_SSc_scRNAseq/dapas_analysis/epi_analysis/'


####  FUNCTIONS  ##############################################################################


#Load Options
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

module_plot <- function(obj, name, markers) {
  marker_list <- list(markers)
  obj <- AddModuleScore(obj, features=marker_list, name='mod_score', assay='RNA')
  return(
    FeaturePlot(obj, features='mod_score1', order=T) & labs(x='UMAP 1', y='UMAP 2') &
      scale_colour_gradientn(colours=brewer.pal(n=11,'Spectral')[11:1], name='Avg. Gene\nExpression\n', 
                             breaks=c(min(obj@meta.data$mod_score1),max(obj@meta.data$mod_score1)),
                             labels=c('Min', 'Max'), guide=guide_colorbar(frame.colour='black', ticks=F)) & 
      ggtitle(name,subtitle=paste(markers, collapse=', ')) &
      theme(plot.title=element_text(face='bold', size=10), plot.subtitle=element_text(face='italic', size=8, hjust=0.5),
            axis.text=element_blank(), axis.title=element_text(size=8), axis.ticks = element_blank(), 
            legend.title=element_text(size=7, face='bold'), legend.text=element_text(size=6.5))
  )
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


MASC <- function(dataset, cluster, contrast, random_effects = NULL, fixed_effects = NULL, verbose = FALSE) {
  # Check inputs
  #dataset=test; cluster=test$layers; contrast='ssc_hc';random_effects='orig.ident'; fixed_effects='location'
  cluster <- cluster[!is.na(dataset[[contrast]])]
  dataset <- dataset[!is.na(dataset[[contrast]]),]
  
  if (is.factor(dataset[[contrast]]) == FALSE) {
    dataset[[contrast]] <- factor(dataset[[contrast]])
  }
  
  # Generate design matrix from cluster assignments
  cluster <- as.character(cluster)
  designmat <- model.matrix(~ cluster + 0, data.frame(cluster = cluster))
  dataset <- cbind(designmat, dataset)
  
  # Convert cluster assignments to string
  cluster <- as.character(cluster)
  # Prepend design matrix generated from cluster assignments
  designmat <- model.matrix(~ cluster + 0, data.frame(cluster = cluster))
  dataset <- cbind(designmat, dataset)
  # Create output list to hold results
  res <- vector(mode = "list", length = length(unique(cluster)))
  names(res) <- attributes(designmat)$dimnames[[2]]
  
  # Create model formulas
  if (!is.null(fixed_effects) && !is.null(random_effects)) {
    model_rhs <- paste0(c(paste0(fixed_effects, collapse="+"),
                          paste0("(1|", random_effects, ")", collapse="+")),
                        collapse = " + ")
    message(paste("Using null model:", "cluster~", model_rhs))
  } else if (is.null(fixed_effects) && !is.null(random_effects)) {
    model_rhs <- paste0("(1|", random_effects, ")", collapse="+")
    message(paste("Using null model:", "cluster~", model_rhs))
  } else {
    model_rhs <- "1" # only includes intercept
    message(paste("Using null model:", "cluster~", model_rhs))
    stop("No random or fixed effects specified")
  }
  
  # Initialize list to store model objects for each cluster
  cluster_models <- vector(mode = "list",
                           length = length(attributes(designmat)$dimnames[[2]]))
  names(cluster_models) <- attributes(designmat)$dimnames[[2]]
  
  # Run nested mixed-effects models for each cluster
  for (i in seq_along(attributes(designmat)$dimnames[[2]])) {
    test_cluster <- attributes(designmat)$dimnames[[2]][i]
    message(paste("Creating logistic mixed models for", test_cluster))
    null_fm <- as.formula(paste0(c(paste0(test_cluster, "~1+"), model_rhs), collapse = ""))
    full_fm <- as.formula(paste0(c(paste0(test_cluster, "~", contrast, "+"), model_rhs), collapse = ""))
    # Run null and full mixed-effects models
    null_model <- lme4::glmer(formula = null_fm, data = dataset,
                              family = binomial, nAGQ = 1, verbose = 0,
                              control = glmerControl(optimizer = "bobyqa"))
    full_model <- lme4::glmer(formula = full_fm, data = dataset,
                              family = binomial, nAGQ = 1, verbose = 0,
                              control = glmerControl(optimizer = "bobyqa"))
    model_lrt <- anova(null_model, full_model)
    # calculate confidence intervals for contrast term beta
    contrast_lvl2 <- paste0(contrast, levels(dataset[[contrast]])[2])
    contrast_ci <- confint.merMod(full_model, method = "Wald", parm = contrast_lvl2)
    # Save model objects to list
    cluster_models[[i]]$null_model <- null_model
    cluster_models[[i]]$full_model <- full_model
    cluster_models[[i]]$model_lrt <- model_lrt
    cluster_models[[i]]$confint <- contrast_ci
  }
  
  # Organize results into output dataframe
  output <- data.frame(cluster = attributes(designmat)$dimnames[[2]], size = colSums(designmat))
  output$model.pvalue <- sapply(cluster_models, function(x) x$model_lrt[["Pr(>Chisq)"]][2])
  output[[paste(contrast_lvl2, "OR", sep = ".")]] <- sapply(cluster_models, function(x) exp(fixef(x$full)[[contrast_lvl2]]))
  output[[paste(contrast_lvl2, "OR", "95pct.ci.lower", sep = ".")]] <- sapply(cluster_models, function(x) exp(x$confint[contrast_lvl2, "2.5 %"]))
  output[[paste(contrast_lvl2, "OR", "95pct.ci.upper", sep = ".")]] <- sapply(cluster_models, function(x) exp(x$confint[contrast_lvl2, "97.5 %"]))
  
  # Return MASC results
  return(output)
}

####  MAIN  ##############################################################################

#super_set <- loadRData(paste0(output_dir, '/MP_eso_epi_reintegratedObj.RData'))

#pdf(file=paste0(output_dir, '/epi_super_violin_nCountFilter.pdf'), width=8, height=4)
ggplot(super_set@meta.data, aes(x = condition, y=foxm1_targets1, fill=condition))+
  geom_violin(scale='width') + labs(y= "FOXM1 Targets", x= "Condition") + ggtitle('Superficial FOXM1 Target Gene Expression') + 
  theme_classic() + theme(axis.title=element_text(face='bold',size=8), axis.ticks = element_blank(), 
                          legend.position='none', panel.grid.major.y = element_line(), panel.grid.minor.y=element_line())
#dev.off()/

#epi_filtered <- subset(super_set, subset=nCount_RNA>3500)
#epi_filtered <- RunPCA(object=epi_filtered, assay='integrated')
#ElbowPlot(epi_filtered, ndims=50)
#save(epi_filtered, file=paste0(output_dir, '/MP_eso_epiFiltered_integratedObj.RData'))

#epi_set <- loadRData(paste0(output_dir, '/MP_eso_epi_reintegratedObj.RData'))
super_set <- loadRData(paste0(output_dir, '/MP_eso_epi_super.RData'))
#super_set <- RunPCA(object=super_set, assay='integrated')


pdf(file=paste0(output_dir, '/epi_super_PC_plot_integratedAssay.pdf'), width=5, height=6)
ElbowPlot(super_set, ndims=50)
dev.off()

pdf(file=paste0(output_dir, '/epi_super_UMAP_condition.pdf'), width=7, height=6)
DimPlot(super_set, reduction='umap', label=F, group.by='condition', shuffle=T, raster=T, pt.size=0.01,
        cols=c('#FF00FF','#00FFFF','#FFFF00')) + ggtitle('') &
  theme(plot.title=element_text(face='bold'  , size=10), plot.subtitle=element_text(face='italic', size=8, hjust=0.5),
        axis.text=element_blank(), axis.title=element_text(size=8), axis.ticks = element_blank(), 
        legend.title=element_text(size=8, face='bold'), legend.text=element_text(size=8))
dev.off()


pdf(file=paste0(output_dir, '/epi_super_UMAP_location.pdf'), width=7, height=6)
DimPlot(super_set, reduction='umap', label=F, group.by='location', shuffle=T, 
        cols=c('#FF00FF','#00FFFF','#FFFF00')) + ggtitle('') &
  theme(plot.title=element_text(face='bold'  , size=10), plot.subtitle=element_text(face='italic', size=8, hjust=0.5),
        axis.text=element_blank(), axis.title=element_text(size=8), axis.ticks = element_blank(), 
        legend.title=element_text(size=8, face='bold'), legend.text=element_text(size=8))
dev.off()



## Plot distribution of each individual sample within integrated UMAP
samples <- names(table(super_set@meta.data$orig.ident))
for (sample in samples) {
  super_set@meta.data$sampleBin <- ifelse(super_set@meta.data$orig.ident==sample,1,0)
  pdf(file=paste0(output_dir, '/Epi_UMAP_SampleDistributions/epi_super_dimPlot_',sample,'.pdf'), width=5, height=5)
  print(DimPlot(super_set, reduction='umap', label=F, group.by='sampleBin', order=T,
                cols=c('grey80','black')) + ggtitle(sample) + theme(legend.position = "none"))
  dev.off()
}

metallo_genes <- c("MT1A","MT1E","MT1F","MT1G","MT1H","MT1M","MT1X")

#super_set <- AddModuleScore(super_set, features=metallo_genes, name='metal_score', assay='RNA')

metal_plots <- lapply(c('Proximal','Distal'), function(loc) {
  lapply(c('HC','SSc'), function(cond) {
    print(c(loc, cond))
    return(FeaturePlot(subset(super_set, location==loc & condition==cond), features='metal_score1', order=T) & labs(x='UMAP 1', y='UMAP 2') &
             scale_colour_gradientn(colours=brewer.pal(n=11,'Spectral')[11:1], name='Avg. Gene\nExpression\n', 
                                    breaks=c(min(super_set@meta.data$metal_score1),max(super_set@meta.data$metal_score1)),
                                    labels=c('Min', 'Max'), guide=guide_colorbar(frame.colour='black', ticks=F)) & 
             ggtitle(paste(loc, ',',cond),subtitle=paste(metallo_genes, collapse=', ')) &
             theme(plot.title=element_text(face='bold', size=10), plot.subtitle=element_text(face='italic', size=8, hjust=0.5),
                   axis.text=element_blank(), axis.title=element_text(size=8), axis.ticks = element_blank(), 
                   legend.title=element_text(size=7, face='bold'), legend.text=element_text(size=6.5)))
  })
})

metal_plots <- c(metal_plots[[1]], metal_plots[[2]])

pdf(file=paste0(output_dir, '/epi_super_modulePlots_metalloGenes.pdf'), width=8, height=8)
print(plot_grid(plotlist=metal_plots, ncol=2))
dev.off()

foxm1_target_genes <- c('CDC20', 'NUF2', 'UBE2T', 'NUCKS1', 'CENPF', 'LBR', 'CKAP2L', 'CCNA2', 'CCNB1', 
                        'LMNB1', 'KIF20A', 'PTTG1', 'MDC1', 'CDCA2', 'CDK1', 'KIF20B', 'KIF11', 'CEP55', 
                        'MKI67', 'CKAP5', 'TROAP', 'RACGAP1', 'NUP37', 'HSP90B1', 'CKAP2', 'MZT1', 'BORA', 
                        'CDKN3', 'KNSTRN', 'CCNB2', 'PRC1', 'PLK1', 'CDC25B', 'TPX2', 'ASXL1', 'UBE2C', 'GTSE1', 'NEK2')

super_set <- AddModuleScore(super_set, features=list(foxm1_target_genes), name='foxm1_targets', assay='RNA')

foxm1_plots <- lapply(c('Proximal','Distal'), function(loc) {
  lapply(c('HC','SSc'), function(cond) {
    print(c(loc, cond))
    return(FeaturePlot(subset(super_set, location==loc & condition==cond), features='foxm1_targets1', order=T) & labs(x='UMAP 1', y='UMAP 2') &
             scale_colour_gradientn(colours=c('grey',brewer.pal(n=11,'Spectral')[5:1]), name='Avg. Gene\nExpression\n', 
                                    breaks=c(min(super_set@meta.data$foxm1_targets1),max(super_set@meta.data$foxm1_targets1)),
                                    labels=c('Min', 'Max'), guide=guide_colorbar(frame.colour='black', ticks=F)) & 
             ggtitle(paste(loc, ',',cond),subtitle=paste(foxm1_target_genes, collapse=', ')) &
             theme(plot.title=element_text(face='bold', size=10), plot.subtitle=element_text(face='italic', size=8, hjust=0.5),
                   axis.text=element_blank(), axis.title=element_text(size=8), axis.ticks = element_blank(), 
                   legend.title=element_text(size=7, face='bold'), legend.text=element_text(size=6.5)))
  })
})

foxm1_plots <- c(foxm1_plots[[1]], foxm1_plots[[2]])

pdf(file=paste0(output_dir, '/epi_super_modulePlots_FOXM1_UpGenes.pdf'), width=8, height=8)
print(plot_grid(plotlist=foxm1_plots, ncol=2))
dev.off()


super_set <- RunUMAP(object=super_set, assay='integrated', dims = 1:10, n.neighbors=35L, min.dist=0.3,
                     spread=1, repulsion.strength=1, negative.sample.rate = 10L, n.epochs=500)
DefaultAssay(super_set) <- 'integrated'
super_set <- FindNeighbors(super_set, reduction = "pca", dims = 1:10, k.param=35, n.trees=5000)

#library(clustree)
for (i in seq(0.05,0.5,0.05)) {
  super_set <- FindClusters(super_set, resolution=i, algorithm=1, n.iter=50, n.start=50)
}

pdf(file=paste0(output_dir, '/EpiSuper_clusterTree_n50.pdf'), width=6, height=8)
clustree(super_set, prefix='integrated_snn_res.')
dev.off()

pdf(file=paste0(output_dir, '/EpiSuper_UMAP_clusters_0.25.pdf'), width=7, height=6)
DimPlot(super_set, reduction = 'umap', label=F, group.by='integrated_snn_res.0.25',
        order=as.character(c(20:0))) + ggtitle('')
dev.off()   
#super_set@reductions$umap@cell.embeddings <- super_set@reductions$umap@cell.embeddings*-1

super_set$seurat_clusters <- super_set$integrated_snn_res.0.25

pdf(file=paste0(output_dir, '/epi_super_super_UMAP_condition.pdf'), width=7, height=6)
DimPlot(super_set, reduction='umap', label=F, group.by='condition', shuffle=T, raster=F, pt.size=0.01,
        cols=c('#FF00FF','#00FFFF','#FFFF00')) + ggtitle('') &
  theme(plot.title=element_text(face='bold'  , size=10), plot.subtitle=element_text(face='italic', size=8, hjust=0.5),
        axis.text=element_blank(), axis.title=element_text(size=8), axis.ticks = element_blank(), 
        legend.title=element_text(size=8, face='bold'), legend.text=element_text(size=8))
dev.off()

pdf(file=paste0(output_dir, '/epi_super_super_UMAP_location.pdf'), width=7, height=6)
DimPlot(super_set, reduction='umap', label=F, group.by='location', shuffle=T, 
        cols=c('#FF00FF','#00FFFF','#FFFF00')) + ggtitle('') &
  theme(plot.title=element_text(face='bold'  , size=10), plot.subtitle=element_text(face='italic', size=8, hjust=0.5),
        axis.text=element_blank(), axis.title=element_text(size=8), axis.ticks = element_blank(), 
        legend.title=element_text(size=8, face='bold'), legend.text=element_text(size=8))
dev.off()

Idents(super_set) <- 'integrated_snn_res.0.25'
DefaultAssay(super_set) <- 'RNA'
super_markers <- FindAllMarkers(super_set)
super_markers$diff.pct <- super_markers$pct.1 - super_markers$pct.2

write.table(super_markers, paste0(output_dir, '/epi_super_clusterMarkers.tsv'), 
            quote=F, row.names=F, col.names=T, sep='\t')

super_set <- AddModuleScore(super_set, features=list(foxm1_target_genes), name='foxm1_targets', assay='RNA')
super_set <- AddModuleScore(super_set, features=list(metallo_genes), name='metallothioneins', assay='RNA')

for (i in c(0:4)) {
  temp_df <- super_markers[order(super_markers$diff.pct,decreasing=T),]
  print(head(temp_df[temp_df$cluster==i,]))
}

markers <- c('FLG','RNASE7','IL36A',
             'CTSV','NAGK','GBP6',
             'metallothioneins1','foxm1_targets1','BAK1')

markers <- c('TRNP1','IL36A','PDCD4')

markers <- c('metallothioneins1')


markers <- c('metallothioneins1','FOXM1','foxm1_targets1','BAK1')


#DefaultAssay(super_set) <- 'RNA'
for (loc in c('Proximal','Distal')) {
  pdf(file=paste0(output_dir, '/epi_super_markers_split_',loc,'.pdf'), width=12, height=13)
  FeaturePlot(subset(super_set, location==loc), markers, order=T, split.by='condition') & 
    scale_colour_gradientn(colours=brewer.pal(n=9,'YlGnBu')[3:9]) & 
    theme(axis.text=element_blank(), axis.title=element_text(size=9),
          axis.ticks = element_blank())
  dev.off()
}

pdf(file=paste0(output_dir, '/epi_super_markers.pdf'), width=12, height=10)
FeaturePlot(super_set, markers, order=T) & 
  scale_colour_gradientn(colours=brewer.pal(n=9,'YlGnBu')[3:9]) & 
  theme(axis.text=element_blank(), axis.title=element_text(size=9),
        axis.ticks = element_blank())
dev.off()

foxm1_plots <- lapply(c('Proximal','Distal'), function(loc) {
  lapply(c('HC','SSc'), function(cond) {
    print(c(loc, cond))
    return(FeaturePlot(subset(super_set, location==loc & condition==cond), features='foxm1_targets1', order=T) & labs(x='UMAP 1', y='UMAP 2') &
             scale_colour_gradientn(colours=c('grey85','grey85',brewer.pal(n=9,'OrRd')[3:9],'black'), name='Avg. Gene\nExpression\n', 
                                    limits=c(min(super_set@meta.data$foxm1_targets1),max(super_set@meta.data$foxm1_targets1)),
                                    breaks=c(min(super_set@meta.data$foxm1_targets1),max(super_set@meta.data$foxm1_targets1)),
                                    labels=c('Min', 'Max'), guide=guide_colorbar(frame.colour='black', ticks=F)) & 
             ggtitle(paste(loc, ',',cond),subtitle='Upregulated FOXM1 Targets') &
             theme(plot.title=element_text(face='bold', size=10), plot.subtitle=element_text(face='italic', size=8, hjust=0.5),
                   axis.text=element_blank(), axis.title=element_text(size=8), axis.ticks = element_blank(), 
                   legend.title=element_text(size=7, face='bold'), legend.text=element_text(size=6.5)))
  })
})

legend <- get_legend(foxm1_plots[[2]][[2]])

foxm1_plots_grid <- plot_grid(foxm1_plots[[1]][[1]]+theme(legend.position='none'),
                              foxm1_plots[[1]][[2]]+theme(legend.position='none'),
                              foxm1_plots[[2]][[1]]+theme(legend.position='none'),
                              foxm1_plots[[2]][[2]]+theme(legend.position='none'), ncol=2)

pdf(file=paste0(output_dir, '/EpiSuper_modulePlots_FOXM1_UpGenes.pdf'), width=8, height=8)
print(plot_grid(foxm1_plots_grid,legend,rel_widths=c(0.9,0.1),ncol=2))
dev.off()

library(gridExtra)

metal_plots <- lapply(c('Proximal','Distal'), function(loc) {
  lapply(c('HC','SSc'), function(cond) {
    print(c(loc, cond))
    return(FeaturePlot(subset(super_set, location==loc & condition==cond), features='metallothioneins1', order=T) & labs(x='UMAP 1', y='UMAP 2') &
             scale_colour_gradientn(colours=c('grey85','grey85',brewer.pal(n=9,'OrRd')[3:9],'black'), name='Avg. Gene\nExpression\n', 
                                    limits=c(min(super_set@meta.data$metallothioneins1),max(super_set@meta.data$metallothioneins1)),
                                    breaks=c(min(super_set@meta.data$metallothioneins1),max(super_set@meta.data$metallothioneins1)),
                                    labels=c('Min', 'Max'), guide=guide_colorbar(frame.colour='black', ticks=F)) & 
             ggtitle(paste(loc, ',',cond),subtitle='Metallothioneins') &
             theme(plot.title=element_text(face='bold', size=10), plot.subtitle=element_text(face='italic', size=8, hjust=0.5),
                   axis.text=element_blank(), axis.title=element_text(size=8), axis.ticks = element_blank(), 
                   legend.title=element_text(size=7, face='bold'), legend.text=element_text(size=6.5)))
  })
})


legend <- get_legend(metal_plots[[2]][[2]])

metal_plots_grid <- plot_grid(metal_plots[[1]][[1]]+theme(legend.position='none'),
                              metal_plots[[1]][[2]]+theme(legend.position='none'),
                              metal_plots[[2]][[1]]+theme(legend.position='none'),
                              metal_plots[[2]][[2]]+theme(legend.position='none'), ncol=2)

pdf(file=paste0(output_dir, '/EpiSuper_modulePlots_metallo.pdf'), width=8, height=8)
print(plot_grid(metal_plots_grid,legend,rel_widths=c(0.9,0.1),ncol=2))
dev.off()

legend <- get_legend(metal_plots[[2]]+theme(legend.position = 'right',
                              legend.title=element_text(face='bold', size=8), 
                              legend.text=element_text(size=8)))

legend, rel_heights=c(1,0.1), scale=0.95,ncol=1

gene <- foxm1_target_genes[1]
gene <- 'CTSC'

ggarrange(
  FeaturePlot(subset(super_set,location=='Proximal'), gene,split.by='condition', order=T) &
    scale_colour_gradientn(colours=c('grey75',brewer.pal(n=11,'YlOrRd')[11:1])),
  FeaturePlot(subset(super_set,location=='Distal'), gene,split.by='condition', order=T) &
    scale_colour_gradientn(colours=c('grey75',brewer.pal(n=11,'YlOrRd')[11:1])),
  ncol=1)

pdf(file=paste0(output_dir, '/epi_super_modulePlots.pdf'), width=12, height=4)
print(plot_grid(plotlist=module_plots, ncol=3))
dev.off()


pdf(file=paste0(output_dir, '/epi_super_UMAP_nFeature_RNA.pdf'), width=7, height=6)
print(FeaturePlot(super_set, features='nFeature_RNA', order=T)  & 
        scale_colour_gradientn(colours=brewer.pal(n=9,'YlGnBu')[3:9]))
dev.off()


plan('multisession', workers=4)
DefaultAssay(super_set) <- 'integrated'
super_set <- FindNeighbors(super_set, reduction = "pca", dims = 1:35, k.param=100, n.trees=50)
super_set <- FindClusters(super_set, resolution=0.5, algorithm=1, n.iter=10, n.start=10)
#super_set <- FindClusters(super_set, resolution=0.75, algorithm=1, n.iter=10, n.start=10)


pdf(file=paste0(output_dir, '/epi_super_UMAP_clusters_k50_res0.25.pdf'), width=7, height=6)
DimPlot(super_set, reduction = 'umap', label=T, group.by='integrated_snn_res.0.25',
        order=as.character(c(30:0))) + ggtitle('') &
  theme(plot.title=element_text(face='bold', size=10), plot.subtitle=element_text(face='italic', size=8, hjust=0.5),
        axis.text=element_blank(), axis.title=element_text(size=8), axis.ticks = element_blank(), 
        legend.title=element_text(size=7, face='bold'), legend.text=element_text(size=6.5))
dev.off()

epi_repset <- subset(super_set, subset=integrated_snn_res.0.25 %in% c(1,4,5,6,7))
save(epi_repset, file=paste0(output_dir, '/MP_eso_epi_repCells.RData'))


Idents(super_set) <- 'integrated_snn_res.0.25'
x <- FindMarkers(super_set, ident.1='0',ident.2='2', min.pct=0.1, logfc.threshold=0.2, only.pos=T)


## Plot stacked violin of markers by cell type
markers <- c('PDPN','COL17A1','KRT14','PCNA','MKI67','KRT13',
             'KRT6A','SERPINB3','CNFN','KRT78','FLG')
#DefaultAssay(sc_obj) <- 'RNA'
#markers <- c('KRT6A','CLU','C19orf33','DEFB1','MMP7','KRT7','KRT14','KRT78')
#reorder <- c(1,6,10,2,3,5,4,7,8,9)

reorder <- c(6,1,8,7,3,4,5,2)
super_set@meta.data$seurat_clusters <- factor(mapvalues(super_set@meta.data$integrated_snn_res.0.25, 
                                                      from=as.character(c(0:7)), to=reorder),
                                            levels=c(1:8))

pdf(file=paste0(output_dir, '/epi_super_UMAP_clusters_k50_res0.25_ordered.pdf'), width=7, height=6)
DimPlot(super_set, reduction = 'umap', label=T, group.by='seurat_clusters',
        order=as.character(c(30:0))) + ggtitle('') &
  theme(plot.title=element_text(face='bold', size=10), plot.subtitle=element_text(face='italic', size=8, hjust=0.5),
        axis.text=element_blank(), axis.title=element_text(size=8), axis.ticks = element_blank(), 
        legend.title=element_text(size=7, face='bold'), legend.text=element_text(size=6.5))
dev.off()




## Plot stacked violin of markers by cell type
markers <- c('FLG','RNASE7',
             'CTSV','metallothioneins1',
             'GBP6','NAGK')


DefaultAssay(super_set) <- 'RNA'

pdf(file=paste0(output_dir, '/epi_super_clusters_markers_violin.pdf'), width=5, height=8)
print(Stacked_VlnPlot(super_set, markers, pt.size=0, group.by='seurat_clusters', raster=F, plot_spacing=0,
                      colors=gg_color_hue(5), x_lab_rotate=T) &
        theme(axis.line.y = element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()))
dev.off()


module_plot(super_set, 'Basal', c('PDPN, IGFBP3'))



library("RColorBrewer")

super_set@meta.data$clusters_ordered <- mapvalues(super_set@meta.data$integrated_snn_res.0.5, 
                                                from=as.character(c(0:13)),
                                                to=c('1','12','10','14','6','8','3','9','13','5','2','7','11','4'))

super_set@meta.data$clusters_ordered <- mapvalues(super_set@meta.data$integrated_snn_res.0.75, 
                                                from=as.character(c(0:20)),
                                                to=c('1','16','20','14','4','10','15','8','3','7','17',
                                                     '9','21','12','18','13','11','6','19','5','2'))

pdf(file=paste0(output_dir, '/epi_super_UMAP_clusters_k40_res0.75_ordered.pdf'), width=7, height=6)
DimPlot(super_set, reduction = 'umap', label=T, group.by='clusters_ordered',
        order=as.character(c(30:0))) + ggtitle('')
dev.off()

layers <- c(rep('Basal', 4), rep('Suprabasal', 7), rep('Superficial', 2))
layer_cols <- c('#E071A6','#5B99D0','#66B564')

super_set@meta.data$layers <- factor(mapvalues(super_set@meta.data$clusters_ordered, from=as.character(c(1:21)),
                                             to=layers), levels=c('Basal','Suprabasal','Superficial'))

pdf(file=paste0(output_dir, '/epi_super_UMAP_layers2.pdf'), width=7, height=6)
DimPlot(super_set, reduction = 'umap', label=F, group.by='layers', 
        cols=c('#E071A6','#5B99D0','#66B564')) + ggtitle('')
dev.off()

layer_clusters <- c('Basal_1','Basal_2','Basal_3','Basal_4','Basal_5','Suprabasal_1',
                    'Suprabasal_2','Superficial_1','Superficial_2')

cluster_cols <- alpha(c(brewer.pal(n=9,'Greens')[c((8:4))], brewer.pal(n=8,'GnBu')[c(5,6)],
                        brewer.pal(n=8,'BuPu')[c(5,6)]), 0.5)

super_set@meta.data$layer_clusters <- factor(mapvalues(super_set@meta.data$integrated_snn_res.0.2, 
                                                     from=as.character(c(0:7)),
                                                     to=layer_clusters), levels=layer_clusters)

super_set@meta.data$compartment <- factor(mapvalues(super_set@meta.data$integrated_snn_res.0.2, from=as.character(c(0:7)),
                                                  to=c('Basal','Basal','Basal','Basal','Basal','Suprabasal',
                                                       'Suprabasal','Superficial','Superficial')), 
                                        levels=c('Basal','Suprabasal','Superficial'))



samples <- names(table(super_set@meta.data$orig.ident))

epi_layer_df <- data.frame(matrix(ncol=4, nrow=length(samples)))
colnames(epi_layer_df) <- c('Sample','Condition','Layer','Proportion')
i=1

for (sample in samples) {
  c <- super_set@meta.data[super_set@meta.data$orig.ident==sample, cond][1]
  total <- nrow(super_set@meta.data[super_set$orig.ident==sample,])
  for (l in c(1:8)) {
    epi_layer_df[i,] <- c(sample, c, l, round(nrow(super_set@meta.data[super_set$orig.ident==sample & 
                                                                       super_set$seurat_clusters==l,])/total,2))
    i <- i+1
  }
}

epi_layer_df$Proportion <- as.numeric(epi_layer_df$Proportion)

epi_layer_df$Condition <- factor(epi_layer_df$Condition, levels=c('HC','GERD','SSc'))

epi_layer_df



pdf(file=paste0(output_dir, '/epi_super_layerProp_vs_condition.pdf'), width=8, height=6)
print(ggbarplot(epi_layer_df, x='Layer', y='Proportion', add=c('mean_sd','jitter'), fill='Condition', 
                position=position_dodge(0.8), add.params=list(color='Condition')) + 
        scale_color_manual(values=rep('black',3)))
dev.off()

non_super <- subset(super_set, layers_rep!='superficial')
props <- as.data.frame(table(non_super@meta.data[,c('condition','layers_rep','location')]))[c(1:12,16:27),]

totals_p <- rowSums(table(non_super@meta.data[non_super@meta.data$location=='Proximal',c('condition','layers_rep')]))
totals_d <- rowSums(table(non_super@meta.data[non_super@meta.data$location=='Distal',c('condition','layers_rep')]))

props$percent_loc <- ifelse(props$location=='Proximal',
                            apply(props[props$location=='Proximal',], 1, function(x) {
                              as.numeric(x['Freq'])/as.numeric(totals_p[x['condition']])})*100,
                            apply(props[props$location=='Distal',], 1, function(x) {
                              as.numeric(x['Freq'])/as.numeric(totals_d[x['condition']])})*100)
props$condition <- factor(props$condition, levels=c('HC','GERD','SSc'))

totals_p_basal <- rowSums(table(non_super@meta.data[non_super@meta.data$location=='Proximal',
                                                    c('condition','layers_rep')])[,1:2])
totals_p_suprabasal <- rowSums(table(non_super@meta.data[non_super@meta.data$location=='Proximal',
                                                         c('condition','layers_rep')])[,3:4])
totals_d_basal <- rowSums(table(non_super@meta.data[non_super@meta.data$location=='Distal',
                                                    c('condition','layers_rep')])[,1:2])
totals_d_suprabasal <- rowSums(table(non_super@meta.data[non_super@meta.data$location=='Distal',
                                                         c('condition','layers_rep')])[,3:4])

pdf(file=paste0(output_dir, '/epi_super_repLayer_vs_condition_loc.pdf'), width=4.5, height=7)
ggarrange(ggplot(props[props$location=='Proximal',], aes(x=condition, y=percent_loc, fill=layers_rep)) + 
            geom_bar(stat='identity', colour='black', linewidth=0.2) + theme_minimal() + ggtitle('Proximal') +
            geom_text(aes(label=paste(round(percent_loc,1),'%')), position=position_stack(vjust=0.5)) +
            scale_fill_manual('Cell type', values=c(layer_cols_5)) + 
            theme(axis.title.x=element_blank(), axis.text.x=element_text(face='bold', color=cond_cols[3:1]),
                  axis.title.y=element_blank()),
          ggplot(props[props$location=='Distal',], aes(x=condition, y=percent_loc, fill=layers_rep)) + 
            geom_bar(stat='identity', colour='black', linewidth=0.2) + theme_minimal() + ggtitle('Distal') +
            geom_text(aes(label=paste(round(percent_loc,1),'%')), position=position_stack(vjust=0.5)) +
            scale_fill_manual('Cell type', values=c(layer_cols_5)) + 
            theme(axis.title.x=element_blank(), axis.text.x=element_text(face='bold',color=cond_cols[3:1]),
                  axis.title.y=element_blank()), 
          ncol=1, common.legend=T, legend='right')
dev.off()


nonSuper_layerRep_df <- data.frame(matrix(ncol=5, nrow=length(samples)))
colnames(nonSuper_layerRep_df) <- c('Sample','Location','Condition','LayerRep','PropRep')
i=1
for (sample in samples) {
  dat <- non_super@meta.data[non_super@meta.data$orig.ident==sample,]
  loc <- ifelse(substr(sample,nchar(sample),nchar(sample))=='P','Proximal','Distal')
  c <- dat[1, cond]
  total <- nrow(dat)
  for (l in names(table(dat$layers_rep)[1:4])) {
    total <- ifelse(l %in% c('basal','replicating basal'),
                    sum(table(dat[dat$location==loc & dat$condition==c,'layers_rep'])[1:2]),
                    sum(table(dat[dat$location==loc & dat$condition==c,'layers_rep'])[3:4]))
    nonSuper_layerRep_df[i,] <- c(sample, loc, as.character(c), l, 
                                  round(nrow(dat[dat$layers_rep==l,])/total,4))
    i <- i+1
  }
}

nonSuper_layerRep_df$PropRep <- as.numeric(nonSuper_layerRep_df$PropRep)
nonSuper_layerRep_df$Condition <- factor(nonSuper_layerRep_df$Condition, levels=c('HC','GERD','SSc'))
nonSuper_layerRep_df$LayerRep <- gsub(" ", '\n', capitalize(nonSuper_layerRep_df$LayerRep))
nonSuper_layerRep_df$LayerRep <- gsub("Replicating", 'Proliferating', capitalize(nonSuper_layerRep_df$LayerRep))

write.table(nonSuper_layerRep_df, paste0(output_dir, '/epi_nonSuper_layerRep_RepProps.tsv'), 
            quote=F, row.names=F, col.names=T, sep='\t')


# Compare proportions
prop_bar_cond <- function(dat, pal, groupBy, splitBy) {
  ggbarplot(dat, x=groupBy, y='PropRep', add=c('mean_sd','jitter'), fill=splitBy, 
            position=position_dodge(0.8), add.params=list(color=splitBy), palette=pal) + 
    scale_color_manual(values=rep('black',3)) + labs(y='Proportion Proliferating') +
    theme(axis.title.x=element_blank(), axis.text.y=element_text(size=8),
          axis.text.x=element_text(face='bold',color=cond_cols[3:1]),
          plot.title=element_text(face='bold',hjust=0.5))
}


pdf(file=paste0(output_dir, '/epi_super_layerRepProp_vs_condition.pdf'), width=7, height=4)
prop_bar_cond(nonSuper_layerRep_df[!nonSuper_layerRep_df$LayerRep %in% c('Basal','Suprabasal'),], 
              layer_cols_5[2:3], 'Condition', 'LayerRep')
dev.off()

p1 <- prop_bar_cond(nonSuper_layerRep_df[!nonSuper_layerRep_df$LayerRep %in% c('Basal','Suprabasal') &
                                           nonSuper_layerRep_df$Location=='Proximal',], 
                    layer_cols_5[2:3], 'Condition', 'LayerRep') + ggtitle('Proximal')
p2 <- prop_bar_cond(nonSuper_layerRep_df[!nonSuper_layerRep_df$LayerRep %in% c('Basal','Suprabasal') &
                                           nonSuper_layerRep_df$Location=='Distal',], 
                    layer_cols_5[2:3], 'Condition', 'LayerRep') + ggtitle('Distal')
legend <- get_legend(p1+theme(legend.box.margin=margin(0,0,0,0), legend.position = 'bottom',
                              legend.title=element_blank()))
p_bars <- plot_grid(p1+theme(legend.position='none', plot.margin=unit(c(0,0,1,0), 'lines')),
                    p2+theme(legend.position='none', plot.margin=unit(c(1,0,0,0), 'lines')), ncol=1)

pdf(file=paste0(output_dir, '/epi_super_layerRepProp_vs_condition_Region.pdf'), width=3.8, height=10)
plot_grid(p_bars, legend, rel_heights=c(1,0.1), scale=0.95,ncol=1)
dev.off()


#### Calculate cell type proportions by condition
prop_df <- data.frame(matrix(ncol=3, nrow=length(table(super_set$seurat_clusters))))
colnames(prop_df) <- c('Cell type','Condition','Percent')
i=1
for (cond in c('HC','GERD','SSc')) {
  total <- nrow(super_set@meta.data[super_set$condition==cond,])
  for (type in c(unique(super_set$seurat_clusters))) {
    print(type)
    prop_df[i,] <- c(type, cond, 
                     round(nrow(super_set@meta.data[super_set$condition==cond & 
                                                    super_set$seurat_clusters==type,])/total,5))
    i <- i+1
  }
}


prop_df$Percent <- as.numeric(prop_df$Percent)*100
prop_df$Condition <- factor(prop_df$Condition, levels=c('HC','GERD','SSc'))

pdf(file=paste0(output_dir, '/epi_super_layer_vs_condition.pdf'), width=7, height=8)
print(ggplot(prop_df, aes(x=Condition, y=Percent, fill=`Cell type`)) +
        geom_bar(stat='identity') + theme_minimal() +
        geom_text(aes(label=paste(round(Percent,1),'%')), position=position_stack(vjust=0.5)) +
        scale_fill_manual('Layer', values=gg_color_hue(length(unique(super_set$seurat_clusters)))))
dev.off()

##### Differentation score
basal_markers <- c('PDPN', 'IGFBP3','COL17A1','KRT15','DST')
suprabasal_markers <- c('KRT4','KRT13','SERPINB3','DSC2','IVL')
superficial_markers <- c('KRT17','KRT78','FLG','CNFN','CRCT1')

super_set <- AddModuleScore(super_set, features=list(basal_markers), name='basal', assay='RNA')
super_set <- AddModuleScore(super_set, features=list(suprabasal_markers), name='suprabasal', assay='RNA')
super_set <- AddModuleScore(super_set, features=list(superficial_markers), name='superficial', assay='RNA')

#super_set$basal11 <- (mean(super_set$basal1) - super_set$basal1)/sd(super_set$basal1)
#super_set$suprabasal1 <- (mean(super_set$suprabasal1) - super_set$suprabasal1)/sd(super_set$suprabasal1)
#super_set$superficial1 <- (mean(super_set$superficial1) - super_set$superficial1)/sd(super_set$superficial1)

super_set$diff.score <- super_set$superficial1+super_set$suprabasal1-super_set$basal1
#super_set$diff.score <- super_set$suprabasal1-super_set$basal1*2


pdf(file=paste0(output_dir, '/epi_super_UMAP_diff.score.pdf'), width=7, height=7)
print(FeaturePlot(super_set, features='diff.score', order=T) & labs(x='UMAP 1', y='UMAP 2') &
        scale_colour_gradientn(colours=brewer.pal(n=11,'Spectral')[11:1], name='Differentiation\nScore\n', 
                               breaks=c(min(super_set@meta.data$diff.score),max(super_set@meta.data$diff.score)),
                               labels=c('Min', 'Max'), guide=guide_colorbar(frame.colour='black', ticks=F)) & 
        ggtitle('Differentiation Score') &
        theme(plot.title=element_text(face='bold', size=10), plot.subtitle=element_text(face='italic', size=8, hjust=0.5),
              axis.text=element_blank(), axis.title=element_text(size=8), axis.ticks = element_blank(), 
              legend.title=element_text(size=7, face='bold'), legend.text=element_text(size=6.5)))
dev.off()

DefaultAssay(super_set) <- 'RNA'
layer_km <- kmeans(FetchData(super_set, vars=c(basal_markers, 
                                             suprabasal_markers, superficial_markers)), 3,
                   iter.max=1000, nstart=10)

#layer_km <- kmeans(super_set@meta.data[,c('basal1', 'suprabasal1', 'superficial1')],3,
#                   iter.max=1000, nstart=10)

super_set@meta.data$layer_cluster <- layer_km$cluster

DimPlot(super_set, reduction = 'umap', label=T, group.by='layer_cluster',
        raster=F, pt.size = 0.01, cols=c('orange','navy','aquamarine')) + ggtitle('') &
  theme(plot.title=element_text(face='bold', size=10), plot.subtitle=element_text(face='italic', size=8, hjust=0.5),
        axis.text=element_blank(), axis.title=element_text(size=8), axis.ticks = element_blank(), 
        legend.title=element_text(size=7, face='bold'), legend.text=element_text(size=6.5))


rep_subset <- subset(super_set, seurat_clusters %in% c(3,4,5))

DefaultAssay(rep_subset) <- 'RNA'
layer_km <- kmeans(FetchData(rep_subset, vars=c(basal_markers, suprabasal_markers)), 2,
                   iter.max=1000, nstart=10)

#layer_km <- kmeans(rep_subset@meta.data[,c('basal1', 'suprabasal1')],2,
#                   iter.max=1000, nstart=10)

rep_subset@meta.data$layer_cluster <- layer_km$cluster

pdf(file=paste0(output_dir, '/epi_super_UMAP_repSubset_kmeans.pdf'), width=7, height=7)
DimPlot(rep_subset, reduction = 'umap', label=T, group.by='layer_cluster',
        raster=F, pt.size = 0.01, cols=c('orange','navy','aquamarine')) + ggtitle('') &
  theme(plot.title=element_text(face='bold', size=10), plot.subtitle=element_text(face='italic', size=8, hjust=0.5),
        axis.text=element_blank(), axis.title=element_text(size=8), axis.ticks = element_blank(), 
        legend.title=element_text(size=7, face='bold'), legend.text=element_text(size=6.5))
dev.off()

super_set$layers <- ifelse(super_set@meta.data$seurat_clusters %in% c(1,2), 'basal',
                         ifelse(super_set@meta.data$seurat_clusters %in% c(8), 'superficial',
                                ifelse(super_set@meta.data$seurat_clusters %in% c(6,7), 'suprabasal',
                                       ifelse(rownames(super_set@meta.data) %in% 
                                                rownames(rep_subset@meta.data[rep_subset@meta.data$layer_cluster==2,]),'basal',
                                              'suprabasal'))))


pdf(file=paste0(output_dir, '/epi_super_UMAP_layers2.pdf'), width=7, height=6)
DimPlot(super_set, reduction = 'umap', label=F, group.by='layers', 
        cols=c('#E071A6','#5B99D0','#66B564')) + ggtitle('')
dev.off()

# plot basal vs suprabasal module scores--difficult to distinguish
library(scales)
pdf(file=paste0(output_dir, '/epi_super_scatter_suprabasal_basal.pdf'), width=7, height=7)
plot(super_set$suprabasal1, super_set$basal1, pch=16,cex=0.25, col='navy',
     xlab='Suprabasal score', ylab='Basal score')
dev.off()

# use k-means clustering to distinguish layers
superficial_subset <- subset(super_set, seurat_clusters %in% c(6,7,8))

DefaultAssay(superficial_subset) <- 'RNA'
layer_km <- kmeans(FetchData(superficial_subset, vars=c(superficial, suprabasal_markers)), 2,
                   iter.max=1000, nstart=10)

layer_km <- kmeans(superficial_subset@meta.data[,c('superficial1', 'suprabasal1')],2,
                   iter.max=1000, nstart=10)

superficial_subset@meta.data$layer_cluster <- layer_km$cluster

pdf(file=paste0(output_dir, '/epi_super_UMAP_superficialSubset_kmeans.pdf'), width=7, height=7)
DimPlot(superficial_subset, reduction = 'umap', label=T, group.by='layer_cluster',
        raster=F, pt.size = 0.01, cols=c('orange','navy','aquamarine')) + ggtitle('') &
  theme(plot.title=element_text(face='bold', size=10), plot.subtitle=element_text(face='italic', size=8, hjust=0.5),
        axis.text=element_blank(), axis.title=element_text(size=8), axis.ticks = element_blank(), 
        legend.title=element_text(size=7, face='bold'), legend.text=element_text(size=6.5))
dev.off()


super_set$layers <- ifelse(rownames(super_set@meta.data) %in% 
                           rownames(superficial_subset@meta.data[superficial_subset@meta.data$layer_cluster==1,]),'superficial',
                         super_set$layers)

super_set$layers <- factor(super_set$layers, levels=c('basal','suprabasal','superficial'))

layer_cols_3 <- c('#4ca5b1','#e9b85d','#c72f4c')

pdf(file=paste0(output_dir, '/epi_super_UMAP_layers_kmeans.pdf'), width=7, height=6)
DimPlot(super_set, reduction = 'umap', label=F, group.by='layers', 
        cols=layer_cols_3) + ggtitle('') &
  theme(plot.title=element_text(face='bold', size=10), plot.subtitle=element_text(face='italic', size=8, hjust=0.5),
        axis.text=element_blank(), axis.title=element_text(size=8), axis.ticks = element_blank(), 
        legend.title=element_text(size=7, face='bold'), legend.text=element_text(size=8))
dev.off()


super_set$layers_rep <- factor(ifelse((super_set$layers=='basal') & (super_set$Phase %in% c('S','G2M')),
                                    'replicating basal',
                                    ifelse((super_set$layers=='suprabasal') & (super_set$Phase %in% c('S','G2M')),
                                           'replicating suprabasal',
                                           ifelse(super_set$layers=='basal','basal',
                                                  ifelse(super_set$layers=='suprabasal','suprabasal',
                                                         'superficial')))),
                             levels=c('basal','replicating basal','replicating suprabasal',
                                      'suprabasal','superficial'))

layer_cols_5 <- brewer.pal(n=11,'Spectral')[c(10,9,7,4,2)]
#layer_cols_5 <- c("#3288BD", "#66C2A5", "#E6F598","#d6bb7c","#D53E4F")
#layer_cols_5 <- c('#031D44','#43829d','#afdc7b','#d6bb7c','#b0434e')
#layer_cols_5 <- c('#80a776','#bbd0b3','#95c5da','#43829d','#e29e8d')

pdf(file=paste0(output_dir, '/epi_super_UMAP_layersRep_kmeans.pdf'), width=8, height=6)
DimPlot(super_set, reduction = 'umap', label=F, group.by='layers_rep',raster=F,
        cols=layer_cols_5) + ggtitle('') &
  theme(plot.title=element_text(face='bold', size=10), plot.subtitle=element_text(face='italic', size=8, hjust=0.5),
        axis.text=element_blank(), axis.title=element_text(size=8), axis.ticks = element_blank(), 
        legend.title=element_text(size=7, face='bold'), legend.text=element_text(size=12))
dev.off()

par(mfrow=c(4,1))
hist(super_set$basal1, xlim=c(-4,4), breaks=40)
hist(super_set$suprabasal1, xlim=c(-4,4), breaks=40)
hist(super_set$superficial1, xlim=c(-4,4), breaks=40)
hist(super_set$diff.score, xlim=c(-4,4), breaks=40)
dev.off()


pdf(file=paste0(output_dir, '/epi_super_dens_layer_moduleScores.pdf'), width=6, height=8)
par(mfrow=c(3,1))
plot(density(super_set@meta.data[super_set$layers=='basal','basal1']), xlim=c(-2,3), ylim=c(0,3),
     main='Basal Cells',col=layer_cols_3[1],lwd=3, xlab='Module Score')
legend('topright',legend=c('Basal','Suprabasal','Superficial'), title='Module Scores', title.font=2,
       title.adj=-0.02,col=layer_cols_3,lty=1, lwd=3, cex=0.8, box.lty=0)
lines(density(super_set@meta.data[super_set$layers=='suprabasal','basal1']), col=layer_cols_3[2],lwd=3)
lines(density(super_set@meta.data[super_set$layers=='superficial','basal1']), col=layer_cols_3[3],lwd=3)
plot(density(super_set@meta.data[super_set$layers=='basal','suprabasal1']), xlim=c(-2,5), ylim=c(0,2),
     main='Suprabasal Cells',col=layer_cols_3[1],lwd=3,xlab='Module Score')
legend('topright',legend=c('Basal','Suprabasal','Superficial'), title='Module Scores', title.font=2,
       title.adj=-0.02, col=layer_cols_3,lty=1, lwd=3, cex=0.8, box.lty=0)
lines(density(super_set@meta.data[super_set$layers=='suprabasal','suprabasal1']), col=layer_cols_3[2],lwd=3)
lines(density(super_set@meta.data[super_set$layers=='superficial','suprabasal1']), col=layer_cols_3[3],lwd=3)
plot(density(super_set@meta.data[super_set$layers=='basal','superficial1']), xlim=c(-2,6), ylim=c(0,2.5),
     main='Superficial Cells',col=layer_cols_3[1],lwd=3,xlab='Module Score')
legend('topright',legend=c('Basal','Suprabasal','Superficial'), title='Module Scores', title.font=2,
       title.adj=-0.02, col=layer_cols_3,lty=1, lwd=3, cex=0.8, box.lty=0)
lines(density(super_set@meta.data[super_set$layers=='suprabasal','superficial1']), col=layer_cols_3[2],lwd=3)
lines(density(super_set@meta.data[super_set$layers=='superficial','superficial1']), col=layer_cols_3[3],lwd=3)
dev.off()


pdf(file=paste0(output_dir, '/epi_super_dens_layer_diff.score.pdf'), width=6, height=5)
plot(density(super_set@meta.data[super_set$layers=='basal','diff.score']), xlim=c(-4,8), ylim=c(0,1),
     main='',col=layer_cols_3[1],lwd=4, xlab='Differentiation Score', cex.axis=0.85)
legend('topright',legend=c('Basal','Suprabasal','Superficial'), 
       col=layer_cols_3,lty=1, lwd=3, cex=0.9, box.lty=0)
lines(density(super_set@meta.data[super_set$layers=='suprabasal','diff.score']), col=layer_cols_3[2],lwd=4)
lines(density(super_set@meta.data[super_set$layers=='superficial','diff.score']), col=layer_cols_3[3],lwd=4)
dev.off()


##### ANALYZE COMPARTMENTS AND CELLS BY CONDITION
table(super_set@meta.data[,c(cond,'layers')])

#super_set$diff.col <- seq_gradient_pal(brewer.pal(n=11,'Spectral')[11:1])(rescale(super_set$diff.score))

super_set$diff.col <- mapvalues(super_set@meta.data$layers, 
                              from=c('basal','suprabasal','superficial'),
                              to=layer_cols_3)

super_set$diff.col5 <- mapvalues(super_set@meta.data$layers_rep, 
                               from=c('basal','replicating basal','replicating suprabasal','suprabasal','superficial'),
                               to=layer_cols_5)

#do_BeeSwarmPlot(super_set, 'diff.score',group.by='condition', continuous_feature = T)

swarm_dat <- do.call(rbind, lapply(split(super_set@meta.data, super_set@meta.data$condition),
                                   function(x) x[sample(nrow(x), 50000), 
                                                 c('condition','diff.score','diff.col5','layers_rep')]))
#swarm_dat <- super_set@meta.data[,c('diff.score','condition','diff.col')]
swarm_dat$condition <- factor(swarm_dat$condition,levels=c('HC','GERD','SSc'))
cond_cols <- c('#8156B3','#e26b53','#A8A39d')

library(beeswarm)

pdf(file=paste0(output_dir, '/epi_super_swarm_layersRep.pdf'), width=5, height=6)
par(mar=c(5,4,0,2))
beeswarm(diff.score ~ condition, data=swarm_dat, vertical=F, pwcol=as.character(swarm_dat$diff.col5),
         fast=T, method='compactswarm', spacing=0.18, cex=0.32, pch=20, ylab='', xlab='Differentiation Score',
         bty='n',axes=F, xlim=c(0,1), cex.lab=1)
axis(1, cex.axis=0.9)
Map(axis, side=2, las=2,at=1:3, col.axis=cond_cols[3:1], labels=c('HC','GERD','SSc'), font=2,tick=F, lwd=0,cex=1.1)
dev.off()


pdf(file=paste0(output_dir, '/epi_super_dens_layers_condition.pdf'), width=6, height=4)
par(mfrow=c(1,1), mar=c(5,4.5,1,2))
plot(density(super_set@meta.data[super_set$condition=='SSc','diff.score']), xlim=c(0,1), ylim=c(0,4.2),
     main='',col=cond_cols[1],lwd=4, xlab='Differentation Score', cex.axis=0.85)
legend('topright',legend=c('SSc','GERD','HC'), 
       col=cond_cols,lty=1, lwd=3, cex=1, box.lty=0)
lines(density(super_set@meta.data[super_set$condition=='GERD','diff.score']), col=cond_cols[2],lwd=4)
lines(density(super_set@meta.data[super_set$condition=='HC','diff.score']), col=cond_cols[3],lwd=4)
dev.off()

# Plot by location
library(beeswarm)
super_set$diff.col <- mapvalues(super_set@meta.data$layers_rep, 
                              from=c('basal','replicating basal','replicating suprabasal','suprabasal','superficial'),
                              to=layer_cols_5)

epi_subset_p <- subset(super_set, location=='Proximal')
epi_subset_d <- subset(super_set, location=='Distal')
swarm_dat_p <- do.call(rbind, lapply(split(epi_subset_p@meta.data, epi_subset_p@meta.data$condition),
                                     function(x) x[sample(nrow(x), min(table(epi_subset_p$condition))), 
                                                   c('condition','diff.score','diff.col','layers_rep')]))
swarm_dat_d <- do.call(rbind, lapply(split(epi_subset_d@meta.data, epi_subset_d@meta.data$condition),
                                     function(x) x[sample(nrow(x), min(table(epi_subset_d$condition))), 
                                                   c('condition','diff.score','diff.col','layers_rep')]))
#swarm_dat <- super_set@meta.data[,c('diff.score','condition','diff.col')]
swarm_dat_p$condition <- factor(swarm_dat_p$condition,levels=c('HC','GERD','SSc'))
swarm_dat_d$condition <- factor(swarm_dat_d$condition,levels=c('HC','GERD','SSc'))

pdf(file=paste0(output_dir, '/epi_super_swarm_layersRep_loc.pdf'), width=8, height=6)
par(mfrow=c(1,2),mar=c(5,4,2,1))
beeswarm(diff.score ~ condition, data=swarm_dat_p, vertical=F, pwcol=as.character(swarm_dat_p$diff.col),
         fast=T, method='compactswarm', spacing=0.18, cex=0.32, pch=20, ylab='', xlab='Differentiation Score',
         bty='n',axes=F, xlim=c(0,1), cex.lab=1,main='Proximal')
axis(1, cex.axis=0.9)
Map(axis, side=2, las=2,at=1:3, col.axis=cond_cols[3:1], labels=c('HC','GERD','SSc'), font=2,tick=F, lwd=0, cex=1.1)
beeswarm(diff.score ~ condition, data=swarm_dat_d, vertical=F, pwcol=as.character(swarm_dat_d$diff.col),
         fast=T, method='compactswarm', spacing=0.18, cex=0.32, pch=20, ylab='', xlab='Differentiation Score',
         bty='n',axes=F, xlim=c(0,1), cex.lab=1,main='Distal')
axis(1, cex.axis=0.9)
Map(axis, side=2, las=2,at=1:3, col.axis=cond_cols[3:1], labels=c('HC','GERD','SSc'), font=2,tick=F, lwd=0,cex=1.1)
dev.off()


pdf(file=paste0(output_dir, '/epi_super_dens_layers_condition_loc.pdf'), width=10, height=5)
par(mfrow=c(1,2), mar=c(5,4.5,3,2))
plot(density(epi_subset_p@meta.data[epi_subset_p$condition=='SSc','diff.score']), xlim=c(0,1), ylim=c(0,4.2),
     main='Proximal',col=cond_cols[1],lwd=4, xlab='Differentation Score', cex.axis=0.85)
legend('topright',legend=c('SSc','GERD','HC'), 
       col=cond_cols,lty=1, lwd=3, cex=1, box.lty=0)
lines(density(epi_subset_p@meta.data[epi_subset_p$condition=='GERD','diff.score']), col=cond_cols[2],lwd=4)
lines(density(epi_subset_p@meta.data[epi_subset_p$condition=='HC','diff.score']), col=cond_cols[3],lwd=4)
plot(density(epi_subset_d@meta.data[epi_subset_d$condition=='SSc','diff.score']), xlim=c(0,1), ylim=c(0,4.2),
     main='Distal',col=cond_cols[1],lwd=4, xlab='Differentation Score', cex.axis=0.85)
legend('topright',legend=c('SSc','GERD','HC'), 
       col=cond_cols,lty=1, lwd=3, cex=1, box.lty=0)
lines(density(epi_subset_d@meta.data[epi_subset_d$condition=='GERD','diff.score']), col=cond_cols[2],lwd=4)
lines(density(epi_subset_d@meta.data[epi_subset_d$condition=='HC','diff.score']), col=cond_cols[3],lwd=4)
dev.off()


# Get cluster counts by person 
epi_cluster_df <- data.frame(matrix(ncol=7, nrow=length(samples)))
colnames(epi_cluster_df) <- c('Sample','Location','Condition','Cluster','nCells','totCells','Proportion')
i=1
for (sample in samples) {
  dat <- super_set@meta.data[super_set@meta.data$orig.ident==sample,]
  loc <- ifelse(substr(sample,nchar(sample),nchar(sample))=='P','Proximal','Distal')
  c <- dat[1, cond]
  total <- nrow(dat)
  for (l in names(table(dat$seurat_clusters))) {
    nCells <- nrow(dat[dat$seurat_clusters==l,])
    epi_cluster_df[i,] <- c(sample, loc, as.character(c), l, nCells, total,
                            round(nCells/total,2))
    i <- i+1
  }
}

write.table(epi_cluster_df, paste0(output_dir, '/epi_cluster_props.tsv'), 
            quote=F, row.names=F, col.names=T, sep='\t')


# Compare proportions
prop_bar_cond <- function(dat, pal, groupBy, splitBy) {
  ggbarplot(dat, x=groupBy, y='Proportion', add=c('mean_sd','jitter'), fill=splitBy, 
            position=position_dodge(0.8), add.params=list(color=splitBy), palette=pal) + 
    scale_color_manual(values=rep('black',3)) + 
    theme(axis.title.x=element_blank(), axis.text.y=element_text(size=8),
          plot.title=element_text(face='bold',hjust=0.5))
}

samples <- names(table(super_set@meta.data$orig.ident))
epi_layer_df <- data.frame(matrix(ncol=5, nrow=length(samples)))
colnames(epi_layer_df) <- c('Sample','Location','Condition','Layer','Proportion')
i=1
for (sample in samples) {
  dat <- super_set@meta.data[super_set@meta.data$orig.ident==sample,]
  loc <- ifelse(substr(sample,nchar(sample),nchar(sample))=='P','Proximal','Distal')
  c <- dat[1, cond]
  total <- nrow(dat)
  for (l in c('basal','suprabasal','superficial')) {
    epi_layer_df[i,] <- c(sample, loc, as.character(c), l, round(nrow(dat[dat$layers==l,])/total,2))
    i <- i+1
  }
}
epi_layer_df$Proportion <- as.numeric(epi_layer_df$Proportion)
epi_layer_df$Condition <- factor(epi_layer_df$Condition, levels=c('HC','GERD','SSc'))
epi_layer_df$Layer <- capitalize(epi_layer_df$Layer)

write.table(epi_layer_df, paste0(output_dir, '/epi_layer_props.tsv'), 
            quote=F, row.names=F, col.names=T, sep='\t')


pdf(file=paste0(output_dir, '/epi_super_layerProp_vs_condition.pdf'), width=6, height=4)
prop_bar_cond(epi_layer_df, cond_cols[3:1], groupBy='Layer', splitBy='Condition')
dev.off()

p1 <- prop_bar_cond(epi_layer_df[epi_layer_df$Location=='Proximal',], 
                    cond_cols[3:1], 'Layer', 'Condition') + ggtitle('Proximal')
p2 <- prop_bar_cond(epi_layer_df[epi_layer_df$Location=='Distal',], 
                    cond_cols[3:1], 'Layer', 'Condition') + ggtitle('Distal')
legend <- get_legend(p1+theme(legend.box.margin=margin(0,0,0,12), legend.position = 'right',
                              legend.title=element_text(face='bold')))
p_bars <- plot_grid(p1+theme(legend.position='none', plot.margin=unit(c(0,0,1,0), 'lines')),
                    p2+theme(legend.position='none', plot.margin=unit(c(1,0,0,0), 'lines')), ncol=1)

pdf(file=paste0(output_dir, '/epi_super_layerProp_vs_condition_Region.pdf'), width=7, height=7)
plot_grid(p_bars, legend, rel_widths=c(1,0.2), scale=0.95)
dev.off()


epi_layerRep_df <- data.frame(matrix(ncol=5, nrow=length(samples)))
colnames(epi_layerRep_df) <- c('Sample','Location','Condition','LayerRep','Proportion')
i=1
for (sample in samples) {
  dat <- super_set@meta.data[super_set@meta.data$orig.ident==sample,]
  loc <- ifelse(substr(sample,nchar(sample),nchar(sample))=='P','Proximal','Distal')
  c <- dat[1, cond]
  total <- nrow(dat)
  for (l in names(table(dat$layers_rep))) {
    epi_layerRep_df[i,] <- c(sample, loc, as.character(c), l, round(nrow(dat[dat$layers_rep==l,])/total,2))
    i <- i+1
  }
}

epi_layerRep_df$Proportion <- as.numeric(epi_layerRep_df$Proportion)
epi_layerRep_df$Condition <- factor(epi_layerRep_df$Condition, levels=c('HC','GERD','SSc'))
epi_layerRep_df$LayerRep <- gsub(" ", '\n', capitalize(epi_layerRep_df$LayerRep))

write.table(epi_layerRep_df, paste0(output_dir, '/epi_layerRep_props.tsv'), 
            quote=F, row.names=F, col.names=T, sep='\t')


pdf(file=paste0(output_dir, '/epi_super_layerRepProp_vs_condition.pdf'), width=7, height=4)
prop_bar_cond(epi_layerRep_df, cond_cols[3:1], 'LayerRep', 'Condition')
dev.off()

p1 <- prop_bar_cond(epi_layerRep_df[epi_layerRep_df$Location=='Proximal',], 
                    cond_cols[3:1], 'LayerRep', 'Condition') + ggtitle('Proximal')
p2 <- prop_bar_cond(epi_layerRep_df[epi_layerRep_df$Location=='Distal',], 
                    cond_cols[3:1], 'LayerRep', 'Condition') + ggtitle('Distal')
legend <- get_legend(p1+theme(legend.box.margin=margin(0,0,0,12), legend.position = 'right',
                              legend.title=element_text(face='bold')))
p_bars <- plot_grid(p1+theme(legend.position='none', plot.margin=unit(c(0,0,1,0), 'lines')),
                    p2+theme(legend.position='none', plot.margin=unit(c(1,0,0,0), 'lines')), ncol=1)

pdf(file=paste0(output_dir, '/epi_super_layerRepProp_vs_condition_Region.pdf'), width=7, height=7)
plot_grid(p_bars, legend, rel_widths=c(1,0.2), scale=0.95)
dev.off()

#### Calculate cell type proportions by condition
prop_df <- data.frame(matrix(ncol=3, nrow=length(table(super_set$seurat_clusters))))
colnames(prop_df) <- c('Cell type','Condition','Percent')
i=1
for (cond in c('HC','GERD','SSc')) {
  total <- nrow(super_set@meta.data[super_set$condition==cond,])
  for (type in c(unique(super_set$layers))) {
    print(type)
    prop_df[i,] <- c(type, cond, 
                     round(nrow(super_set@meta.data[super_set$condition==cond & 
                                                    super_set$layers==type,])/total,5))
    i <- i+1
  }
}


prop_df$Percent <- as.numeric(prop_df$Percent)*100
prop_df$Condition <- factor(prop_df$Condition, levels=c('HC','GERD','SSc'))

pdf(file=paste0(output_dir, '/epi_super_layer_vs_condition.pdf'), width=7, height=8)
print(ggplot(prop_df, aes(x=Condition, y=Percent, fill=`Cell type`)) +
        geom_bar(stat='identity') + theme_minimal() +
        geom_text(aes(label=paste(round(Percent,1),'%')), position=position_stack(vjust=0.5)) +
        scale_fill_manual('Layer', values=gg_color_hue(length(unique(super_set$seurat_clusters)))))
dev.off()

### Compute p-values for proportional differences
test_df <- super_set@meta.data[,c('orig.ident','location','condition','layers','layers_rep')]
test_df$layers_rep <- gsub(' ','_',test_df$layers_rep)
test_df$ssc_hc <- ifelse(test_df$condition=='SSc',1,ifelse(test_df$condition=='HC',0,NA))
test_df$gerd_hc <- ifelse(test_df$condition=='GERD',1,ifelse(test_df$condition=='HC',0,NA))
test_df$ssc_gerd <- ifelse(test_df$condition=='SSc',1,ifelse(test_df$condition=='GERD',0,NA))

prop_comps <- data.frame(matrix(ncol=9, nrow=0))
colnames(prop_comps) <- c('cluster','size','masc.p','OR','OR.95.ci.lower','OR.95.ci.upper',
                          'comp','location','propeller.p')

for (comp in c('ssc_hc','gerd_hc','ssc_gerd')) {
  for (l in c('All','Proximal','Distal')) {
    if (l=='All') {
      temp_df <- test_df
    } else {
      temp_df <- test_df[test_df$location==l,]
    }
    masc <- MASC(temp_df, cluster=temp_df$layers_rep, contrast=comp, random_effects='orig.ident')
    masc$cluster <- gsub('cluster','',masc$cluster)
    masc$comp <- comp
    masc$location <- l
    propel <- propeller(clusters=temp_df$layers_rep, sample=temp_df$orig.ident, group=temp_df[,comp])
    masc <- merge(masc, propel[,c('BaselineProp.clusters','P.Value')], by.x='cluster', 
                  by.y='BaselineProp.clusters')
    prop_comps <- rbind(prop_comps, setNames(masc, names(prop_comps)))
  }
}

prop_comps <- prop_comps[,c(7,8,1,2,4:6,3,9)]
prop_comps[prop_comps$location=='All', 'masc.p.fdr'] <- p.adjust(prop_comps[prop_comps$location=='All', 
                                                                            'masc.p'], method='fdr')
prop_comps[prop_comps$location!='All', 'masc.p.fdr'] <- p.adjust(prop_comps[prop_comps$location!='All', 
                                                                            'masc.p'], method='fdr')
prop_comps[prop_comps$location=='All', 'propeller.p.fdr'] <- p.adjust(prop_comps[prop_comps$location=='All', 
                                                                                 'propeller.p'], method='fdr')
prop_comps[prop_comps$location!='All', 'propeller.p.fdr'] <- p.adjust(prop_comps[prop_comps$location!='All', 
                                                                                 'propeller.p'], method='fdr')
write.table(prop_comps, paste0(output_dir, '/epi_layerRep_props_stats.tsv'), 
            quote=F, row.names=F, col.names=T, sep='\t')


### Generate separate UMAPs for each condition

epi_hc <- subset(super_set, subset=condition=='HC')
epi_hc <- epi_hc[, sample(colnames(epi_hc), 50000)]
epi_gerd <- subset(super_set, subset=condition=='GERD')
epi_gerd <- epi_gerd[, sample(colnames(epi_gerd), 50000)]
epi_ssc <- subset(super_set, subset=condition=='SSc')
epi_ssc <- epi_ssc[, sample(colnames(epi_ssc), 50000)]


epi_hc <- RunUMAP(object=epi_hc, assay='integrated', dims = 1:35, n.neighbors=100L, min.dist=0.3,
                  spread=1, repulsion.strength=0.9, negative.sample.rate = 10L, n.epochs=500)

epi_gerd <- RunUMAP(object=epi_gerd, assay='integrated', dims = 1:35, n.neighbors=100L, min.dist=0.3,
                    spread=1, repulsion.strength=0.9, negative.sample.rate = 10L, n.epochs=500)

epi_ssc <- RunUMAP(object=epi_ssc, assay='integrated', dims = 1:35, n.neighbors=100L, min.dist=0.3,
                   spread=1, repulsion.strength=0.9, negative.sample.rate = 10L, n.epochs=500)


epi_hc@reductions$umap@cell.embeddings[,2] <- epi_hc@reductions$umap@cell.embeddings[,2]*-1
epi_gerd@reductions$umap@cell.embeddings[,1] <- epi_gerd@reductions$umap@cell.embeddings[,1]*-1
epi_ssc@reductions$umap@cell.embeddings <- epi_ssc@reductions$umap@cell.embeddings*-1



p1 <- FeaturePlot(epi_hc, features='diff.score', order=T, pt.size=0.01) & labs(x='UMAP 1', y='UMAP 2') &
  scale_colour_gradientn(colours=brewer.pal(n=11,'Spectral')[11:1], name='Differentiation\nScore\n', 
                         breaks=c(min(super_set@meta.data$diff.score),max(super_set@meta.data$diff.score)),
                         labels=c('Min', 'Max'), guide=guide_colorbar(frame.colour='black', ticks=F)) & 
  ggtitle('HC') & theme(plot.title=element_text(face='bold', size=10), axis.text=element_blank(), 
                        axis.title=element_text(size=8), axis.ticks = element_blank())
p2 <- FeaturePlot(epi_gerd, features='diff.score', order=T, pt.size=0.01) & labs(x='UMAP 1', y='UMAP 2') &
  scale_colour_gradientn(colours=brewer.pal(n=11,'Spectral')[11:1], 
                         breaks=c(min(super_set@meta.data$diff.score),max(super_set@meta.data$diff.score)),
                         guide=guide_colorbar(frame.colour='black', ticks=F)) & ggtitle('GERD') &
  theme(plot.title=element_text(face='bold', size=10), axis.text=element_blank(), 
        axis.title=element_text(size=8), axis.ticks = element_blank())
p3 <- FeaturePlot(epi_ssc, features='diff.score', order=T, pt.size=0.01) & labs(x='UMAP 1', y='UMAP 2') &
  scale_colour_gradientn(colours=brewer.pal(n=11,'Spectral')[11:1], 
                         breaks=c(min(super_set@meta.data$diff.score),max(super_set@meta.data$diff.score)),
                         guide=guide_colorbar(frame.colour='black', ticks=F)) & ggtitle('SSc') &
  theme(plot.title=element_text(face='bold', size=10), axis.text=element_blank(), 
        axis.title=element_text(size=8), axis.ticks = element_blank(), legend.position='none')

legend <- get_legend(p1+theme(legend.position = 'right',
                              legend.title=element_text(face='bold', size=8), 
                              legend.text=element_text(size=8)))

p_umaps <- plot_grid(p1+theme(legend.position='none'), p2, p3, ncol=3)

pdf(file=paste0(output_dir, '/epi_super_UMAP_diff_condition.pdf'), width=12, height=4)
plot_grid(p_umaps, legend, rel_widths=c(1,0.1), scale=0.95)
dev.off()



### Generate separate UMAPs for each condition for distal and proximal
epi_p <- subset(super_set, location=='Proximal')
epi_d <- subset(super_set, location=='Distal')

epi_p$diff.bin <- cut(epi_p$diff.score, 
                      breaks=unique(quantile(epi_p$diff.score,
                                             probs=seq.int(0,1,by=1/10))), include.lowest=T)

epi_d$diff.bin <- cut(epi_d$diff.score, 
                      breaks=unique(quantile(epi_d$diff.score,
                                             probs=seq.int(0,1,by=1/10))), include.lowest=T)

epi_p$bin_9 <- ifelse(epi_p$diff.bin=='(2.89,4.41]',1,0)
epi_d$bin_9 <- ifelse(epi_d$diff.bin=='(2.72,4.06]',1,0)


epi_hc_p <- subset(epi_p, condition=='HC')
#epi_hc_p <- epi_hc_p[, sample(colnames(epi_hc_p), 24000)]
epi_hc_d <- subset(epi_d, condition=='HC')
#epi_hc_d <- epi_hc_d[, sample(colnames(epi_hc_d), 24000)]

epi_gerd_p <- subset(epi_p, condition=='GERD')
#epi_gerd_p <- epi_gerd_p[, sample(colnames(epi_gerd_p), 24000)]
epi_gerd_d <- subset(epi_d, condition=='GERD')
#epi_gerd_d <- epi_gerd_d[, sample(colnames(epi_gerd_d), 24000)]

epi_ssc_p <- subset(epi_p, condition=='SSc')
#epi_ssc_p <- epi_ssc_p[, sample(colnames(epi_ssc_p), 24000)]
epi_ssc_d <- subset(epi_d, condition=='SSc')
#epi_ssc_d <- epi_ssc_d[, sample(colnames(epi_ssc_d), 24000)]


epi_hc_p <- RunUMAP(object=epi_hc_p, assay='integrated', dims = 1:25, n.neighbors=100L, min.dist=0.3,
                    spread=1, repulsion.strength=0.9, negative.sample.rate = 10L, n.epochs=500)
epi_hc_d <- RunUMAP(object=epi_hc_d, assay='integrated', dims = 1:35, n.neighbors=100L, min.dist=0.3,
                    spread=1, repulsion.strength=0.9, negative.sample.rate = 10L, n.epochs=500)
epi_gerd_p <- RunUMAP(object=epi_gerd_p, assay='integrated', dims = 1:35, n.neighbors=100L, min.dist=0.3,
                      spread=1, repulsion.strength=0.9, negative.sample.rate = 10L, n.epochs=500)
epi_gerd_d <- RunUMAP(object=epi_gerd_d, assay='integrated', dims = 1:35, n.neighbors=100L, min.dist=0.3,
                      spread=1, repulsion.strength=0.9, negative.sample.rate = 10L, n.epochs=500)
epi_ssc_p <- RunUMAP(object=epi_ssc_p, assay='integrated', dims = 1:35, n.neighbors=100L, min.dist=0.3,
                     spread=1, repulsion.strength=0.9, negative.sample.rate = 10L, n.epochs=500)
epi_ssc_d <- RunUMAP(object=epi_ssc_d, assay='integrated', dims = 1:35, n.neighbors=100L, min.dist=0.3,
                     spread=1, repulsion.strength=0.9, negative.sample.rate = 10L, n.epochs=500)

epi_hc_p@reductions$umap@cell.embeddings[,1] <- epi_hc_p@reductions$umap@cell.embeddings[,1]*-1
epi_gerd_p@reductions$umap@cell.embeddings[,2] <- epi_gerd_p@reductions$umap@cell.embeddings[,2]*-1
#epi_ssc_p@reductions$umap@cell.embeddings <- epi_ssc_p@reductions$umap@cell.embeddings*-1

epi_hc_d@reductions$umap@cell.embeddings[,1] <- epi_hc_d@reductions$umap@cell.embeddings[,1]*-1
epi_gerd_d@reductions$umap@cell.embeddings[,1] <- epi_gerd_d@reductions$umap@cell.embeddings[,1]*-1
epi_ssc_d@reductions$umap@cell.embeddings[,2] <- epi_ssc_d@reductions$umap@cell.embeddings[,2]*-1

plot_grid(FeaturePlot(epi_hc_p, feature='diff.score') &
            scale_colour_gradientn(colours=brewer.pal(n=11,'Spectral')[11:1]),
          FeaturePlot(epi_gerd_p, feature='diff.score') &
            scale_colour_gradientn(colours=brewer.pal(n=11,'Spectral')[11:1]),
          FeaturePlot(epi_ssc_p, feature='diff.score') &
            scale_colour_gradientn(colours=brewer.pal(n=11,'Spectral')[11:1]),
          FeaturePlot(epi_hc_d, feature='diff.score') &
            scale_colour_gradientn(colours=brewer.pal(n=11,'Spectral')[11:1]),
          FeaturePlot(epi_gerd_d, feature='diff.score') &
            scale_colour_gradientn(colours=brewer.pal(n=11,'Spectral')[11:1]),
          FeaturePlot(epi_ssc_d, feature='diff.score') &
            scale_colour_gradientn(colours=brewer.pal(n=11,'Spectral')[11:1]), ncol=3)


bin_cols <- c(brewer.pal(n=11,'Spectral')[11:4],'black',brewer.pal(n=11,'Spectral')[2])
plot_grid(DimPlot(epi_hc_p, group.by='diff.bin', cols=bin_cols),
          DimPlot(epi_gerd_p, group.by='diff.bin', cols=bin_cols),
          DimPlot(epi_ssc_p, group.by='diff.bin', cols=bin_cols),
          DimPlot(epi_hc_d, group.by='diff.bin', cols=bin_cols),
          DimPlot(epi_gerd_d, group.by='diff.bin', cols=bin_cols),
          DimPlot(epi_ssc_d, group.by='diff.bin', cols=bin_cols), ncol=3)


bin_cols <- c('grey','black')
plot_grid(DimPlot(epi_hc_p, group.by='bin_9', cols=bin_cols, order=T),
          DimPlot(epi_gerd_p, group.by='bin_9', cols=bin_cols, order=T),
          DimPlot(epi_ssc_p, group.by='bin_9', cols=bin_cols, order=T),
          DimPlot(epi_hc_d, group.by='bin_9', cols=bin_cols, order=T),
          DimPlot(epi_gerd_d, group.by='bin_9', cols=bin_cols, order=T),
          DimPlot(epi_ssc_d, group.by='bin_9', cols=bin_cols, order=T), ncol=3)


epi_hc@reductions$umap@cell.embeddings[,2] <- epi_hc@reductions$umap@cell.embeddings[,2]*-1
epi_gerd@reductions$umap@cell.embeddings[,1] <- epi_gerd@reductions$umap@cell.embeddings[,1]*-1
epi_ssc@reductions$umap@cell.embeddings <- epi_ssc@reductions$umap@cell.embeddings*-1



p1 <- FeaturePlot(epi_hc, features='diff.score', order=T, pt.size=0.01) & labs(x='UMAP 1', y='UMAP 2') &
  scale_colour_gradientn(colours=brewer.pal(n=11,'Spectral')[11:1], name='Differentiation\nScore\n', 
                         breaks=c(min(super_set@meta.data$diff.score),max(super_set@meta.data$diff.score)),
                         labels=c('Min', 'Max'), guide=guide_colorbar(frame.colour='black', ticks=F)) & 
  ggtitle('HC') & theme(plot.title=element_text(face='bold', size=10), axis.text=element_blank(), 
                        axis.title=element_text(size=8), axis.ticks = element_blank())
p2 <- FeaturePlot(epi_gerd, features='diff.score', order=T, pt.size=0.01) & labs(x='UMAP 1', y='UMAP 2') &
  scale_colour_gradientn(colours=brewer.pal(n=11,'Spectral')[11:1], 
                         breaks=c(min(super_set@meta.data$diff.score),max(super_set@meta.data$diff.score)),
                         guide=guide_colorbar(frame.colour='black', ticks=F)) & ggtitle('GERD') &
  theme(plot.title=element_text(face='bold', size=10), axis.text=element_blank(), 
        axis.title=element_text(size=8), axis.ticks = element_blank())
p3 <- FeaturePlot(epi_ssc, features='diff.score', order=T, pt.size=0.01) & labs(x='UMAP 1', y='UMAP 2') &
  scale_colour_gradientn(colours=brewer.pal(n=11,'Spectral')[11:1], 
                         breaks=c(min(super_set@meta.data$diff.score),max(super_set@meta.data$diff.score)),
                         guide=guide_colorbar(frame.colour='black', ticks=F)) & ggtitle('SSc') &
  theme(plot.title=element_text(face='bold', size=10), axis.text=element_blank(), 
        axis.title=element_text(size=8), axis.ticks = element_blank(), legend.position='none')

legend <- get_legend(p1+theme(legend.position = 'right',
                              legend.title=element_text(face='bold', size=8), 
                              legend.text=element_text(size=8)))

p_umaps <- plot_grid(p1+theme(legend.position='none'), p2, p3, ncol=3)

pdf(file=paste0(output_dir, '/epi_super_UMAP_diff_condition.pdf'), width=12, height=4)
plot_grid(p_umaps, legend, rel_widths=c(1,0.1), scale=0.95)
dev.off()



plot_grid(FeaturePlot(epi_hc, features='diff.score', order=T, raster=T) & labs(x='UMAP 1', y='UMAP 2') &
            scale_colour_gradientn(colours=brewer.pal(n=11,'Spectral')[11:1], name='Differentiation\nScore\n', 
                                   breaks=c(min(super_set@meta.data$diff.score),max(super_set@meta.data$diff.score)),
                                   labels=c('Min', 'Max'), guide=guide_colorbar(frame.colour='black', ticks=F)) & 
            ggtitle('HC') &
            theme(plot.title=element_text(face='bold', size=10), plot.subtitle=element_text(face='italic', size=8, hjust=0.5),
                  axis.text=element_blank(), axis.title=element_text(size=8), axis.ticks = element_blank(), 
                  legend.position='none'),
          FeaturePlot(epi_gerd, features='diff.score', order=T, raster=T) & labs(x='UMAP 1', y='UMAP 2') &
            scale_colour_gradientn(colours=brewer.pal(n=11,'Spectral')[11:1], name='Differentiation\nScore\n', 
                                   breaks=c(min(super_set@meta.data$diff.score),max(super_set@meta.data$diff.score)),
                                   labels=c('Min', 'Max'), guide=guide_colorbar(frame.colour='black', ticks=F)) & 
            ggtitle('GERD') &
            theme(plot.title=element_text(face='bold', size=10), plot.subtitle=element_text(face='italic', size=8, hjust=0.5),
                  axis.text=element_blank(), axis.title=element_text(size=8), axis.ticks = element_blank(), 
                  legend.title=element_text(size=7, face='bold'), legend.text=element_text(size=6.5)),
          FeaturePlot(epi_ssc, features='diff.score', order=T, raster=T) & labs(x='UMAP 1', y='UMAP 2') &
            scale_colour_gradientn(colours=brewer.pal(n=11,'Spectral')[11:1], name='Differentiation\nScore\n', 
                                   breaks=c(min(super_set@meta.data$diff.score),max(super_set@meta.data$diff.score)),
                                   labels=c('Min', 'Max'), guide=guide_colorbar(frame.colour='black', ticks=F)) & 
            ggtitle('SSc') &
            theme(plot.title=element_text(face='bold', size=10), plot.subtitle=element_text(face='italic', size=8, hjust=0.5),
                  axis.text=element_blank(), axis.title=element_text(size=8), axis.ticks = element_blank(), 
                  legend.title=element_text(size=7, face='bold'), legend.text=element_text(size=6.5)),
          ncol=3)



#### Look at differential gene expression ####

# basal, suprabasal, superficial
library(future.apply)
plan('multisession',workers=3)
DefaultAssay(super_set) <- 'RNA'
Idents(super_set) <- 'condition'

layers <- c('basal','suprabasal','superficial')
deg_list_ssc_hc <- list(); deg_list_ssc_gerd <- list(); deg_list_gerd_hc <- list()

deg_list_ssc_hc <- setNames(lapply(layers, function(l) {
  FindMarkers(subset(super_set, layers==l), ident.1='SSc', ident.2='HC', logfc.threshold=0.05, min.pct=0.01)}), layers)
deg_list_ssc_gerd <- setNames(lapply(layers, function(l) {
  FindMarkers(subset(super_set, layers==l), ident.1='SSc', ident.2='GERD', logfc.threshold=0.05, min.pct=0.01)}), layers)
deg_list_gerd_hc <- setNames(lapply(layers, function(l) {
  FindMarkers(subset(super_set, layers==l), ident.1='GERD', ident.2='HC', logfc.threshold=0.05, min.pct=0.01)}), layers)

# replicating basal, replicating suprabasal

layers <- names(table(super_set$layers_rep))
deg_list_ssc_hc <- list(); deg_list_ssc_gerd <- list(); deg_list_gerd_hc <- list()

deg_list_ssc_hc <- setNames(future_lapply(layers, function(l) {
  FindMarkers(subset(super_set, layers_rep==l), ident.1='SSc', ident.2='HC', logfc.threshold=0, min.pct=0)}, future.seed=T), layers)
deg_list_ssc_gerd <- setNames(future_lapply(layers, function(l) {
  FindMarkers(subset(super_set, layers_rep==l), ident.1='SSc', ident.2='GERD', logfc.threshold=0, min.pct=0)}, future.seed=T), layers)
deg_list_gerd_hc <- setNames(future_lapply(layers, function(l) {
  FindMarkers(subset(super_set, layers_rep==l), ident.1='GERD', ident.2='HC', logfc.threshold=0, min.pct=0)}, future.seed=T), layers)

# equal size bins by diff.score

super_set$diff.bin <- cut(super_set$diff.score, 
                        breaks=unique(quantile(super_set$diff.score,
                                               probs=seq.int(0,1,by=1/10))), include.lowest=T)

bins <- names(table(super_set$diff.bin))
deg_list_ssc_hc_bins <- list(); deg_list_ssc_hc_bins <- list(); deg_list_gerd_hc_bins <- list()

deg_list_ssc_hc_bins <- setNames(future_lapply(bins, function(l) {
  FindMarkers(subset(super_set, diff.bin==l), ident.1='SSc', ident.2='HC', logfc.threshold=0.2, min.pct=0.1)}, future.seed=T), bins)
deg_list_ssc_gerd_bins <- setNames(future_lapply(bins, function(l) {
  FindMarkers(subset(super_set, diff.bin==l), ident.1='SSc', ident.2='GERD', logfc.threshold=0.2, min.pct=0.1)}, future.seed=T), bins)
deg_list_gerd_hc_bins <- setNames(future_lapply(bins, function(l) {
  FindMarkers(subset(super_set, diff.bin==l), ident.1='GERD', ident.2='HC', logfc.threshold=0.2, min.pct=0.1)}, future.seed=T), bins)

## Venn Diagrams

deg_list_ssc_hc <- loadRData(paste0(output_dir, '/epi_degs_ssc_hc_layers.RData'))
deg_list_gerd_hc <- loadRData(paste0(output_dir, '/epi_degs_gerd_hc_layers.RData'))
deg_list_ssc_gerd <- loadRData(paste0(output_dir, '/epi_degs_ssc_gerd_layers.RData'))

#deg_list_ssc_hc_sig <- deg_list_ssc_hc
#deg_list_gerd_hc_sig <- deg_list_gerd_hc
#deg_list_ssc_gerd_sig <- deg_list_ssc_gerd

# Set standard for DEG inclusion
fc_min <- 0.1
deg_list_ssc_hc_sig <- lapply(deg_list_ssc_hc, subset, abs(avg_log2FC)>fc_min & p_val_adj<0.05)
deg_list_gerd_hc_sig <- lapply(deg_list_gerd_hc, subset, abs(avg_log2FC)>fc_min & p_val_adj<0.05)
deg_list_ssc_gerd_sig <- lapply(deg_list_ssc_gerd, subset, abs(avg_log2FC)>fc_min & p_val_adj<0.05)


rownames(deg_list_ssc_hc_sig$superficial[!rownames(deg_list_ssc_gerd_sig$superficial) %in%
                                           rownames(deg_list_gerd_hc_sig$superficial),])

deg_basal <- c('SSc'=sum(!rownames(deg_list_ssc_hc_sig$basal) %in% 
                           rownames(deg_list_gerd_hc_sig$basal)), 
               'GERD'=sum(!rownames(deg_list_gerd_hc_sig$basal) %in% 
                            rownames(deg_list_ssc_hc_sig$basal)),
               'SSc&GERD'=sum(rownames(deg_list_ssc_hc_sig$basal) %in% 
                                rownames(deg_list_gerd_hc_sig$basal)))
deg_supra <- c('SSc'=sum(!rownames(deg_list_ssc_hc_sig$suprabasal) %in% 
                           rownames(deg_list_gerd_hc_sig$suprabasal)), 
               'GERD'=sum(!rownames(deg_list_gerd_hc_sig$suprabasal) %in% 
                            rownames(deg_list_ssc_hc_sig$suprabasal)),
               'SSc&GERD'=sum(rownames(deg_list_ssc_hc_sig$suprabasal) %in% 
                                rownames(deg_list_gerd_hc_sig$suprabasal)))
deg_super <- c('SSc'=sum(!rownames(deg_list_ssc_hc_sig$superficial) %in% 
                           rownames(deg_list_gerd_hc_sig$superficial)), 
               'GERD'=sum(!rownames(deg_list_gerd_hc_sig$superficial) %in% 
                            rownames(deg_list_ssc_hc_sig$superficial)),
               'SSc&GERD'=sum(rownames(deg_list_ssc_hc_sig$superficial) %in% 
                                rownames(deg_list_gerd_hc_sig$superficial)))

fit_basal <- euler(deg_basal, shape='circle')
fit_supra <- euler(deg_supra, shape='circle')
fit_super <- euler(deg_super, shape='circle')


pdf(file=paste0(output_dir, '/epi_super_degVenn_layers_fc',fc_min,'.pdf'), width=8, height=4)
plot_grid(plot(fit_basal, labels=list(font=4), fills=list(fill=cond_cols, alpha=0.5), cex.main=0.8,
               main='Basal\nDEG vs. HC', quantities=list(type=c('counts','percent'), cex=0.7, font=3, round=3)),
          plot(fit_supra, labels=list(font=4), fills=list(fill=cond_cols, alpha=0.5), cex.main=0.8,
               main='Suprabasal\nDEG vs. HC', quantities=list(type=c('counts','percent'), cex=0.7, font=3, round=3)),
          plot(fit_super, labels=list(font=4), fills=list(fill=cond_cols, alpha=0.5), cex.main=0.8,
               main='Superficial\nDEG vs. HC', quantities=list(type=c('counts','percent'), cex=0.7, font=3, round=3)
          ), ncol=3)
dev.off()



### One plot with all layers for SSc

deg_ssc <- c('Basal'=sum(!rownames(deg_list_ssc_hc_sig$basal)
                         [!rownames(deg_list_ssc_hc_sig$basal) %in% 
                             rownames(deg_list_ssc_hc_sig$suprabasal)]
                         %in% rownames(deg_list_ssc_hc_sig$superficial)),
             'Suprabasal'=sum(!rownames(deg_list_ssc_hc_sig$suprabasal)
                              [!rownames(deg_list_ssc_hc_sig$suprabasal) %in% 
                                  rownames(deg_list_ssc_hc_sig$basal)]
                              %in% rownames(deg_list_ssc_hc_sig$superficial)),
             'Superficial'=sum(!rownames(deg_list_ssc_hc_sig$superficial)
                               [!rownames(deg_list_ssc_hc_sig$superficial) %in% 
                                   rownames(deg_list_ssc_hc_sig$basal)]
                               %in% rownames(deg_list_ssc_hc_sig$suprabasal)),
             'Basal&Suprabasal'=sum(rownames(deg_list_ssc_hc_sig$basal)
                                    [!rownames(deg_list_ssc_hc_sig$basal) %in%
                                        rownames(deg_list_ssc_hc_sig$superficial)] %in% 
                                      rownames(deg_list_ssc_hc_sig$suprabasal)),
             'Suprabasal&Superficial'=sum(rownames(deg_list_ssc_hc_sig$suprabasal)
                                          [!rownames(deg_list_ssc_hc_sig$suprabasal) %in%
                                              rownames(deg_list_ssc_hc_sig$basal)] %in% 
                                            rownames(deg_list_ssc_hc_sig$superficial)),
             'Basal&Superficial'=sum(rownames(deg_list_ssc_hc_sig$basal)
                                     [!rownames(deg_list_ssc_hc_sig$basal) %in%
                                         rownames(deg_list_ssc_hc_sig$suprabasal)] %in% 
                                       rownames(deg_list_ssc_hc_sig$superficial)),
             'Basal&Suprabasal&Superficial'=sum(rownames(deg_list_ssc_hc_sig$basal)
                                                [rownames(deg_list_ssc_hc_sig$basal) %in% 
                                                    rownames(deg_list_ssc_hc_sig$suprabasal)] %in%
                                                  rownames(deg_list_ssc_hc_sig$superficial)))

fit_ssc <- euler(deg_ssc, shape='ellipse')

pdf(file=paste0(output_dir, '/epi_super_degVenn_layers_SSc_fc',fc_min,'.pdf'), width=6, height=4)
plot(fit_ssc, labels=list(font=4), fills=list(fill=layer_cols_3, alpha=0.5), cex.main=0.8,
     main='SSc vs. HC DEGs', quantities=list(type=c('counts','percent'),
                                             cex=0.7, font=3, round=3))
dev.off()


### One plot with all layers for GERD
deg_gerd <- c('Basal'=sum(!rownames(deg_list_gerd_hc_sig$basal)
                          [!rownames(deg_list_gerd_hc_sig$basal) %in% 
                              rownames(deg_list_gerd_hc_sig$suprabasal)]
                          %in% rownames(deg_list_gerd_hc_sig$superficial)),
              'Suprabasal'=sum(!rownames(deg_list_gerd_hc_sig$suprabasal)
                               [!rownames(deg_list_gerd_hc_sig$suprabasal) %in% 
                                   rownames(deg_list_gerd_hc_sig$basal)]
                               %in% rownames(deg_list_gerd_hc_sig$superficial)),
              'Superficial'=sum(!rownames(deg_list_gerd_hc_sig$superficial)
                                [!rownames(deg_list_gerd_hc_sig$superficial) %in% 
                                    rownames(deg_list_gerd_hc_sig$basal)]
                                %in% rownames(deg_list_gerd_hc_sig$suprabasal)),
              'Basal&Suprabasal'=sum(rownames(deg_list_gerd_hc_sig$basal)
                                     [!rownames(deg_list_gerd_hc_sig$basal) %in%
                                         rownames(deg_list_gerd_hc_sig$superficial)] %in% 
                                       rownames(deg_list_gerd_hc_sig$suprabasal)),
              'Suprabasal&Superficial'=sum(rownames(deg_list_gerd_hc_sig$suprabasal)
                                           [!rownames(deg_list_gerd_hc_sig$suprabasal) %in%
                                               rownames(deg_list_gerd_hc_sig$basal)] %in% 
                                             rownames(deg_list_gerd_hc_sig$superficial)),
              'Basal&Superficial'=sum(rownames(deg_list_gerd_hc_sig$basal)
                                      [!rownames(deg_list_gerd_hc_sig$basal) %in%
                                          rownames(deg_list_gerd_hc_sig$suprabasal)] %in% 
                                        rownames(deg_list_gerd_hc_sig$superficial)),
              'Basal&Suprabasal&Superficial'=sum(rownames(deg_list_gerd_hc_sig$basal)
                                                 [rownames(deg_list_gerd_hc_sig$basal) %in% 
                                                     rownames(deg_list_gerd_hc_sig$suprabasal)] %in%
                                                   rownames(deg_list_gerd_hc_sig$superficial)))

fit_gerd <- euler(deg_gerd, shape='ellipse')

pdf(file=paste0(output_dir, '/epi_super_degVenn_layers_GERD_fc',fc_min,'.pdf'), width=6, height=4)
plot(fit_gerd, labels=list(font=4), fills=list(fill=layer_cols_3, alpha=0.5), cex.main=0.8,
     main='GERD vs. HC DEGs', quantities=list(type=c('counts','percent'),
                                              cex=0.7, font=3, round=3))
dev.off()


### One plot with all layers for SSc vs GERD

deg_ssc_gerd <- c('Basal'=sum(!rownames(deg_list_ssc_gerd_sig$basal)
                              [!rownames(deg_list_ssc_gerd_sig$basal) %in% 
                                  rownames(deg_list_ssc_gerd_sig$suprabasal)]
                              %in% rownames(deg_list_ssc_gerd_sig$superficial)),
                  'Suprabasal'=sum(!rownames(deg_list_ssc_gerd_sig$suprabasal)
                                   [!rownames(deg_list_ssc_gerd_sig$suprabasal) %in% 
                                       rownames(deg_list_ssc_gerd_sig$basal)]
                                   %in% rownames(deg_list_ssc_gerd_sig$superficial)),
                  'Superficial'=sum(!rownames(deg_list_ssc_gerd_sig$superficial)
                                    [!rownames(deg_list_ssc_gerd_sig$superficial) %in% 
                                        rownames(deg_list_ssc_gerd_sig$basal)]
                                    %in% rownames(deg_list_ssc_gerd_sig$suprabasal)),
                  'Basal&Suprabasal'=sum(rownames(deg_list_ssc_gerd_sig$basal)
                                         [!rownames(deg_list_ssc_gerd_sig$basal) %in%
                                             rownames(deg_list_ssc_gerd_sig$superficial)] %in% 
                                           rownames(deg_list_ssc_gerd_sig$suprabasal)),
                  'Suprabasal&Superficial'=sum(rownames(deg_list_ssc_gerd_sig$suprabasal)
                                               [!rownames(deg_list_ssc_gerd_sig$suprabasal) %in%
                                                   rownames(deg_list_ssc_gerd_sig$basal)] %in% 
                                                 rownames(deg_list_ssc_gerd_sig$superficial)),
                  'Basal&Superficial'=sum(rownames(deg_list_ssc_gerd_sig$basal)
                                          [!rownames(deg_list_ssc_gerd_sig$basal) %in%
                                              rownames(deg_list_ssc_gerd_sig$suprabasal)] %in% 
                                            rownames(deg_list_ssc_gerd_sig$superficial)),
                  'Basal&Suprabasal&Superficial'=sum(rownames(deg_list_ssc_gerd_sig$basal)
                                                     [rownames(deg_list_ssc_gerd_sig$basal) %in% 
                                                         rownames(deg_list_ssc_gerd_sig$suprabasal)] %in%
                                                       rownames(deg_list_ssc_gerd_sig$superficial)))

fit_ssc_gerd <- euler(deg_ssc_gerd, shape='ellipse')

pdf(file=paste0(output_dir, '/epi_super_degVenn_layers_SScGerd_fc',fc_min,'.pdf'), width=6, height=4)
plot(fit_ssc_gerd, labels=list(font=4), fills=list(fill=layer_cols_3, alpha=0.5), cex.main=0.8,
     main='SSc vs. GERD DEGs', quantities=list(type=c('counts','percent'),
                                               cex=0.7, font=3, round=3))
dev.off()



### One plot with one circle for each DEG comparison

ssc_hc_genes <- unique(c(rownames(deg_list_ssc_hc_sig$basal),
                         rownames(deg_list_ssc_hc_sig$suprabasal), 
                         rownames(deg_list_ssc_hc_sig$superficial)))
gerd_hc_genes <- unique(c(rownames(deg_list_gerd_hc_sig$basal),
                          rownames(deg_list_gerd_hc_sig$suprabasal), 
                          rownames(deg_list_gerd_hc_sig$superficial)))
ssc_gerd_genes <- unique(c(rownames(deg_list_ssc_gerd_sig$basal),
                           rownames(deg_list_ssc_gerd_sig$suprabasal), 
                           rownames(deg_list_ssc_gerd_sig$superficial)))

deg_conds <- c('SSc_HC'=sum(!ssc_hc_genes[!ssc_hc_genes %in% gerd_hc_genes] 
                            %in% ssc_gerd_genes),
               'GERD_HC'=sum(!gerd_hc_genes[!gerd_hc_genes %in% ssc_hc_genes] 
                             %in% ssc_gerd_genes),
               'SSc_GERD'=sum(!ssc_gerd_genes[!ssc_gerd_genes %in% ssc_hc_genes] 
                              %in% gerd_hc_genes),
               'SSc_HC&GERD_HC'=sum(ssc_hc_genes[!ssc_hc_genes %in% ssc_gerd_genes] 
                                    %in% gerd_hc_genes),
               'GERD_HC&SSc_GERD'=sum(gerd_hc_genes[!gerd_hc_genes %in% ssc_hc_genes] 
                                      %in% ssc_gerd_genes),
               'SSc_HC&SSc_GERD'=sum(ssc_hc_genes[!ssc_hc_genes %in% gerd_hc_genes] 
                                     %in% ssc_gerd_genes),
               'SSc_HC&GERD_HC&SSc_GERD'=sum(ssc_hc_genes[ssc_hc_genes %in% gerd_hc_genes] 
                                             %in% ssc_gerd_genes))

fit_conds <- euler(deg_conds, shape='ellipse')

pdf(file=paste0(output_dir, '/epi_super_degVenn_layers_conds_fc',fc_min,'.pdf'), width=6, height=4)
plot(fit_conds, labels=list(font=4), fills=list(fill=c(cond_cols[1:2],'khaki2'), alpha=0.5), cex.main=0.8,
     main='SSc vs. GERD DEGs', quantities=list(type=c('counts','percent'),
                                               cex=0.7, font=3, round=3))
dev.off()


#### Venn Diagrams by proximal and distal location: ####

deg_list_ssc_hc_p <- loadRData(paste0(output_dir, '/epi_degs_ssc_hc_layers_Proximal.RData'))
deg_list_gerd_hc_p <- loadRData(paste0(output_dir, '/epi_degs_gerd_hc_layers_Proximal.RData'))
deg_list_ssc_gerd_p <- loadRData(paste0(output_dir, '/epi_degs_ssc_gerd_layers_Proximal.RData'))
deg_list_ssc_hc_d <- loadRData(paste0(output_dir, '/epi_degs_ssc_hc_layers_Distal.RData'))
deg_list_gerd_hc_d <- loadRData(paste0(output_dir, '/epi_degs_gerd_hc_layers_Distal.RData'))
deg_list_ssc_gerd_d <- loadRData(paste0(output_dir, '/epi_degs_ssc_gerd_layers_Distal.RData'))

# Set standard for DEG inclusion
fc_min <- 0.1
deg_list_ssc_hc_p_sig <- lapply(deg_list_ssc_hc_p, subset, abs(avg_log2FC)>fc_min & p_val_adj<0.05)
deg_list_gerd_hc_p_sig <- lapply(deg_list_gerd_hc_p, subset, abs(avg_log2FC)>fc_min & p_val_adj<0.05)
deg_list_ssc_gerd_p_sig <- lapply(deg_list_ssc_gerd_p, subset, abs(avg_log2FC)>fc_min & p_val_adj<0.05)
deg_list_ssc_hc_d_sig <- lapply(deg_list_ssc_hc_d, subset, abs(avg_log2FC)>fc_min & p_val_adj<0.05)
deg_list_gerd_hc_d_sig <- lapply(deg_list_gerd_hc_d, subset, abs(avg_log2FC)>fc_min & p_val_adj<0.05)
deg_list_ssc_gerd_d_sig <- lapply(deg_list_ssc_gerd_d, subset, abs(avg_log2FC)>fc_min & p_val_adj<0.05)

deg_basal_p <- c('SSc'=sum(!rownames(deg_list_ssc_hc_p_sig$basal) %in% rownames(deg_list_gerd_hc_p_sig$basal)), 
                 'GERD'=sum(!rownames(deg_list_gerd_hc_p_sig$basal) %in% rownames(deg_list_ssc_hc_p_sig$basal)),
                 'SSc&GERD'=sum(rownames(deg_list_ssc_hc_p_sig$basal) %in% rownames(deg_list_gerd_hc_p_sig$basal)))
deg_supra_p <- c('SSc'=sum(!rownames(deg_list_ssc_hc_p_sig$suprabasal) %in% rownames(deg_list_gerd_hc_p_sig$suprabasal)), 
                 'GERD'=sum(!rownames(deg_list_gerd_hc_p_sig$suprabasal) %in% rownames(deg_list_ssc_hc_p_sig$suprabasal)),
                 'SSc&GERD'=sum(rownames(deg_list_ssc_hc_p_sig$suprabasal) %in% rownames(deg_list_gerd_hc_p_sig$suprabasal)))
deg_super_p <- c('SSc'=sum(!rownames(deg_list_ssc_hc_p_sig$superficial) %in% rownames(deg_list_gerd_hc_p_sig$superficial)), 
                 'GERD'=sum(!rownames(deg_list_gerd_hc_p_sig$superficial) %in% rownames(deg_list_ssc_hc_p_sig$superficial)),
                 'SSc&GERD'=sum(rownames(deg_list_ssc_hc_p_sig$superficial) %in% rownames(deg_list_gerd_hc_p_sig$superficial)))

deg_basal_d <- c('SSc'=sum(!rownames(deg_list_ssc_hc_d_sig$basal) %in% rownames(deg_list_gerd_hc_d_sig$basal)), 
                 'GERD'=sum(!rownames(deg_list_gerd_hc_d_sig$basal) %in% rownames(deg_list_ssc_hc_d_sig$basal)),
                 'SSc&GERD'=sum(rownames(deg_list_ssc_hc_d_sig$basal) %in% rownames(deg_list_gerd_hc_d_sig$basal)))
deg_supra_d <- c('SSc'=sum(!rownames(deg_list_ssc_hc_d_sig$suprabasal) %in% rownames(deg_list_gerd_hc_d_sig$suprabasal)), 
                 'GERD'=sum(!rownames(deg_list_gerd_hc_d_sig$suprabasal) %in% rownames(deg_list_ssc_hc_d_sig$suprabasal)),
                 'SSc&GERD'=sum(rownames(deg_list_ssc_hc_d_sig$suprabasal) %in% rownames(deg_list_gerd_hc_d_sig$suprabasal)))
deg_super_d <- c('SSc'=sum(!rownames(deg_list_ssc_hc_d_sig$superficial) %in% rownames(deg_list_gerd_hc_d_sig$superficial)), 
                 'GERD'=sum(!rownames(deg_list_gerd_hc_d_sig$superficial) %in% rownames(deg_list_ssc_hc_d_sig$superficial)),
                 'SSc&GERD'=sum(rownames(deg_list_ssc_hc_d_sig$superficial) %in% rownames(deg_list_gerd_hc_d_sig$superficial)))

fit_basal_p <- euler(deg_basal_p, shape='circle')
fit_supra_p <- euler(deg_supra_p, shape='circle')
fit_super_p <- euler(deg_super_p, shape='circle')

fit_basal_d <- euler(deg_basal_d, shape='circle')
fit_supra_d <- euler(deg_supra_d, shape='circle')
fit_super_d <- euler(deg_super_d, shape='circle')

pdf(file=paste0(output_dir, '/epi_super_degVenn_layers_loc_fc',fc_min,'.pdf'), width=8, height=4)
plot_grid(plot_grid(plot(fit_basal_p, labels=list(font=4), fills=list(fill=cond_cols, alpha=0.5), cex.main=0.7,
                         main='Basal\nDEG vs. HC', quantities=list(type=c('counts','percent'), cex=0.7, font=3, round=3)),
                    plot(fit_supra_p, labels=list(font=4), fills=list(fill=cond_cols, alpha=0.5), cex.main=0.7,
                         main='Suprabasal\nDEG vs. HC', quantities=list(type=c('counts','percent'), cex=0.7, font=3, round=3)),
                    plot(fit_super_p, labels=list(font=4), fills=list(fill=cond_cols, alpha=0.5), cex.main=0.7,
                         main='Superficial\nDEG vs. HC', quantities=list(type=c('counts','percent'), cex=0.7, font=3, round=3)
                    ), ncol=3),
          plot_grid(plot(fit_basal_d, labels=list(font=4), fills=list(fill=cond_cols, alpha=0.5), cex.main=0.7,
                         main='Basal\nDEG vs. HC', quantities=list(type=c('counts','percent'), cex=0.7, font=3, round=3)),
                    plot(fit_supra_d, labels=list(font=4), fills=list(fill=cond_cols, alpha=0.5), cex.main=0.7,
                         main='Suprabasal\nDEG vs. HC', quantities=list(type=c('counts','percent'), cex=0.7, font=3, round=3)),
                    plot(fit_super_d, labels=list(font=4), fills=list(fill=cond_cols, alpha=0.5), cex.main=0.7,
                         main='Superficial\nDEG vs. HC', quantities=list(type=c('counts','percent'), cex=0.7, font=3, round=3)
                    ), ncol=3), ncol=1, labels=c('Proximal','Distal'))
dev.off()



### One plot with all layers for SSc

deg_ssc_p <- c('Basal'=sum(!rownames(deg_list_ssc_hc_p_sig$basal)[!rownames(deg_list_ssc_hc_p_sig$basal) %in% 
                                                                    rownames(deg_list_ssc_hc_p_sig$suprabasal)] %in% rownames(deg_list_ssc_hc_p_sig$superficial)),
               'Suprabasal'=sum(!rownames(deg_list_ssc_hc_p_sig$suprabasal)[!rownames(deg_list_ssc_hc_p_sig$suprabasal) %in% 
                                                                              rownames(deg_list_ssc_hc_p_sig$basal)] %in% rownames(deg_list_ssc_hc_p_sig$superficial)),
               'Superficial'=sum(!rownames(deg_list_ssc_hc_p_sig$superficial)[!rownames(deg_list_ssc_hc_p_sig$superficial) %in% 
                                                                                rownames(deg_list_ssc_hc_p_sig$basal)] %in% rownames(deg_list_ssc_hc_p_sig$suprabasal)),
               'Basal&Suprabasal'=sum(rownames(deg_list_ssc_hc_p_sig$basal)[!rownames(deg_list_ssc_hc_p_sig$basal) %in%
                                                                              rownames(deg_list_ssc_hc_p_sig$superficial)] %in% rownames(deg_list_ssc_hc_p_sig$suprabasal)),
               'Suprabasal&Superficial'=sum(rownames(deg_list_ssc_hc_p_sig$suprabasal)[!rownames(deg_list_ssc_hc_p_sig$suprabasal) %in%
                                                                                         rownames(deg_list_ssc_hc_p_sig$basal)] %in% rownames(deg_list_ssc_hc_p_sig$superficial)),
               'Basal&Superficial'=sum(rownames(deg_list_ssc_hc_p_sig$basal)[!rownames(deg_list_ssc_hc_p_sig$basal) %in%
                                                                               rownames(deg_list_ssc_hc_p_sig$suprabasal)] %in% rownames(deg_list_ssc_hc_p_sig$superficial)),
               'Basal&Suprabasal&Superficial'=sum(rownames(deg_list_ssc_hc_p_sig$basal)[rownames(deg_list_ssc_hc_p_sig$basal) %in% 
                                                                                          rownames(deg_list_ssc_hc_p_sig$suprabasal)] %in% rownames(deg_list_ssc_hc_p_sig$superficial)))

deg_ssc_d <- c('Basal'=sum(!rownames(deg_list_ssc_hc_d_sig$basal)[!rownames(deg_list_ssc_hc_d_sig$basal) %in% 
                                                                    rownames(deg_list_ssc_hc_d_sig$suprabasal)] %in% rownames(deg_list_ssc_hc_d_sig$superficial)),
               'Suprabasal'=sum(!rownames(deg_list_ssc_hc_d_sig$suprabasal)[!rownames(deg_list_ssc_hc_d_sig$suprabasal) %in% 
                                                                              rownames(deg_list_ssc_hc_d_sig$basal)] %in% rownames(deg_list_ssc_hc_d_sig$superficial)),
               'Superficial'=sum(!rownames(deg_list_ssc_hc_d_sig$superficial)[!rownames(deg_list_ssc_hc_d_sig$superficial) %in% 
                                                                                rownames(deg_list_ssc_hc_d_sig$basal)] %in% rownames(deg_list_ssc_hc_d_sig$suprabasal)),
               'Basal&Suprabasal'=sum(rownames(deg_list_ssc_hc_d_sig$basal)[!rownames(deg_list_ssc_hc_d_sig$basal) %in%
                                                                              rownames(deg_list_ssc_hc_d_sig$superficial)] %in% rownames(deg_list_ssc_hc_d_sig$suprabasal)),
               'Suprabasal&Superficial'=sum(rownames(deg_list_ssc_hc_d_sig$suprabasal)[!rownames(deg_list_ssc_hc_d_sig$suprabasal) %in%
                                                                                         rownames(deg_list_ssc_hc_d_sig$basal)] %in% rownames(deg_list_ssc_hc_d_sig$superficial)),
               'Basal&Superficial'=sum(rownames(deg_list_ssc_hc_d_sig$basal)[!rownames(deg_list_ssc_hc_d_sig$basal) %in%
                                                                               rownames(deg_list_ssc_hc_d_sig$suprabasal)] %in% rownames(deg_list_ssc_hc_d_sig$superficial)),
               'Basal&Suprabasal&Superficial'=sum(rownames(deg_list_ssc_hc_d_sig$basal)[rownames(deg_list_ssc_hc_d_sig$basal) %in% 
                                                                                          rownames(deg_list_ssc_hc_d_sig$suprabasal)] %in% rownames(deg_list_ssc_hc_d_sig$superficial)))

fit_ssc_p <- euler(deg_ssc_p, shape='circle')
fit_ssc_d <- euler(deg_ssc_d, shape='circle')


pdf(file=paste0(output_dir, '/epi_super_degVenn_layers_loc_SSc_fc',fc_min,'.pdf'), width=6, height=4)
plot_grid(plot(fit_ssc_p, labels=list(font=4), fills=list(fill=layer_cols_3, alpha=0.5), cex.main=0.8,
               main='SSc vs. HC DEGs\nProximal', quantities=list(type=c('counts','percent'),
                                                                 cex=0.7, font=3, round=3)),
          plot(fit_ssc_d, labels=list(font=4), fills=list(fill=layer_cols_3, alpha=0.5), cex.main=0.8,
               main='SSc vs. HC DEGs\nDistal', quantities=list(type=c('counts','percent'),
                                                               cex=0.7, font=3, round=3)), ncol=2)
dev.off()


### One plot with all layers for GERD
deg_gerd_p <- c('Basal'=sum(!rownames(deg_list_gerd_hc_p_sig$basal)
                            [!rownames(deg_list_gerd_hc_p_sig$basal) %in% 
                                rownames(deg_list_gerd_hc_p_sig$suprabasal)]
                            %in% rownames(deg_list_gerd_hc_p_sig$superficial)),
                'Suprabasal'=sum(!rownames(deg_list_gerd_hc_p_sig$suprabasal)
                                 [!rownames(deg_list_gerd_hc_p_sig$suprabasal) %in% 
                                     rownames(deg_list_gerd_hc_p_sig$basal)]
                                 %in% rownames(deg_list_gerd_hc_p_sig$superficial)),
                'Superficial'=sum(!rownames(deg_list_gerd_hc_p_sig$superficial)
                                  [!rownames(deg_list_gerd_hc_p_sig$superficial) %in% 
                                      rownames(deg_list_gerd_hc_p_sig$basal)]
                                  %in% rownames(deg_list_gerd_hc_p_sig$suprabasal)),
                'Basal&Suprabasal'=sum(rownames(deg_list_gerd_hc_p_sig$basal)
                                       [!rownames(deg_list_gerd_hc_p_sig$basal) %in%
                                           rownames(deg_list_gerd_hc_p_sig$superficial)] %in% 
                                         rownames(deg_list_gerd_hc_p_sig$suprabasal)),
                'Suprabasal&Superficial'=sum(rownames(deg_list_gerd_hc_p_sig$suprabasal)
                                             [!rownames(deg_list_gerd_hc_p_sig$suprabasal) %in%
                                                 rownames(deg_list_gerd_hc_p_sig$basal)] %in% 
                                               rownames(deg_list_gerd_hc_p_sig$superficial)),
                'Basal&Superficial'=sum(rownames(deg_list_gerd_hc_p_sig$basal)
                                        [!rownames(deg_list_gerd_hc_p_sig$basal) %in%
                                            rownames(deg_list_gerd_hc_p_sig$suprabasal)] %in% 
                                          rownames(deg_list_gerd_hc_p_sig$superficial)),
                'Basal&Suprabasal&Superficial'=sum(rownames(deg_list_gerd_hc_p_sig$basal)
                                                   [rownames(deg_list_gerd_hc_p_sig$basal) %in% 
                                                       rownames(deg_list_gerd_hc_p_sig$suprabasal)] %in%
                                                     rownames(deg_list_gerd_hc_p_sig$superficial)))
deg_gerd_d <- c('Basal'=sum(!rownames(deg_list_gerd_hc_d_sig$basal)
                            [!rownames(deg_list_gerd_hc_d_sig$basal) %in% 
                                rownames(deg_list_gerd_hc_d_sig$suprabasal)]
                            %in% rownames(deg_list_gerd_hc_d_sig$superficial)),
                'Suprabasal'=sum(!rownames(deg_list_gerd_hc_d_sig$suprabasal)
                                 [!rownames(deg_list_gerd_hc_d_sig$suprabasal) %in% 
                                     rownames(deg_list_gerd_hc_d_sig$basal)]
                                 %in% rownames(deg_list_gerd_hc_d_sig$superficial)),
                'Superficial'=sum(!rownames(deg_list_gerd_hc_d_sig$superficial)
                                  [!rownames(deg_list_gerd_hc_d_sig$superficial) %in% 
                                      rownames(deg_list_gerd_hc_d_sig$basal)]
                                  %in% rownames(deg_list_gerd_hc_d_sig$suprabasal)),
                'Basal&Suprabasal'=sum(rownames(deg_list_gerd_hc_d_sig$basal)
                                       [!rownames(deg_list_gerd_hc_d_sig$basal) %in%
                                           rownames(deg_list_gerd_hc_d_sig$superficial)] %in% 
                                         rownames(deg_list_gerd_hc_d_sig$suprabasal)),
                'Suprabasal&Superficial'=sum(rownames(deg_list_gerd_hc_d_sig$suprabasal)
                                             [!rownames(deg_list_gerd_hc_d_sig$suprabasal) %in%
                                                 rownames(deg_list_gerd_hc_d_sig$basal)] %in% 
                                               rownames(deg_list_gerd_hc_d_sig$superficial)),
                'Basal&Superficial'=sum(rownames(deg_list_gerd_hc_d_sig$basal)
                                        [!rownames(deg_list_gerd_hc_d_sig$basal) %in%
                                            rownames(deg_list_gerd_hc_d_sig$suprabasal)] %in% 
                                          rownames(deg_list_gerd_hc_d_sig$superficial)),
                'Basal&Suprabasal&Superficial'=sum(rownames(deg_list_gerd_hc_d_sig$basal)
                                                   [rownames(deg_list_gerd_hc_d_sig$basal) %in% 
                                                       rownames(deg_list_gerd_hc_d_sig$suprabasal)] %in%
                                                     rownames(deg_list_gerd_hc_d_sig$superficial)))

fit_gerd_p <- euler(deg_gerd_p, shape='circle')
fit_gerd_d <- euler(deg_gerd_d, shape='circle')

pdf(file=paste0(output_dir, '/epi_super_degVenn_layers_loc_GERD_fc',fc_min,'.pdf'), width=6, height=4)
plot_grid(plot(fit_gerd_p, labels=list(font=4), fills=list(fill=layer_cols_3, alpha=0.5), cex.main=0.8,
               main='GERD vs. HC DEGs\nProximal', quantities=list(type=c('counts','percent'),
                                                                  cex=0.7, font=3, round=3)),
          plot(fit_gerd_d, labels=list(font=4), fills=list(fill=layer_cols_3, alpha=0.5), cex.main=0.8,
               main='GERD vs. HC DEGs\nDistal', quantities=list(type=c('counts','percent'),
                                                                cex=0.7, font=3, round=3)), ncol=2)
dev.off()


### One plot with all layers for SSc vs GERD

deg_ssc_gerd_p <- c('Basal'=sum(!rownames(deg_list_ssc_gerd_p_sig$basal)
                                [!rownames(deg_list_ssc_gerd_p_sig$basal) %in% 
                                    rownames(deg_list_ssc_gerd_p_sig$suprabasal)]
                                %in% rownames(deg_list_ssc_gerd_p_sig$superficial)),
                    'Suprabasal'=sum(!rownames(deg_list_ssc_gerd_p_sig$suprabasal)
                                     [!rownames(deg_list_ssc_gerd_p_sig$suprabasal) %in% 
                                         rownames(deg_list_ssc_gerd_p_sig$basal)]
                                     %in% rownames(deg_list_ssc_gerd_p_sig$superficial)),
                    'Superficial'=sum(!rownames(deg_list_ssc_gerd_p_sig$superficial)
                                      [!rownames(deg_list_ssc_gerd_p_sig$superficial) %in% 
                                          rownames(deg_list_ssc_gerd_p_sig$basal)]
                                      %in% rownames(deg_list_ssc_gerd_p_sig$suprabasal)),
                    'Basal&Suprabasal'=sum(rownames(deg_list_ssc_gerd_p_sig$basal)
                                           [!rownames(deg_list_ssc_gerd_p_sig$basal) %in%
                                               rownames(deg_list_ssc_gerd_p_sig$superficial)] %in% 
                                             rownames(deg_list_ssc_gerd_p_sig$suprabasal)),
                    'Suprabasal&Superficial'=sum(rownames(deg_list_ssc_gerd_p_sig$suprabasal)
                                                 [!rownames(deg_list_ssc_gerd_p_sig$suprabasal) %in%
                                                     rownames(deg_list_ssc_gerd_p_sig$basal)] %in% 
                                                   rownames(deg_list_ssc_gerd_p_sig$superficial)),
                    'Basal&Superficial'=sum(rownames(deg_list_ssc_gerd_p_sig$basal)
                                            [!rownames(deg_list_ssc_gerd_p_sig$basal) %in%
                                                rownames(deg_list_ssc_gerd_p_sig$suprabasal)] %in% 
                                              rownames(deg_list_ssc_gerd_p_sig$superficial)),
                    'Basal&Suprabasal&Superficial'=sum(rownames(deg_list_ssc_gerd_p_sig$basal)
                                                       [rownames(deg_list_ssc_gerd_p_sig$basal) %in% 
                                                           rownames(deg_list_ssc_gerd_p_sig$suprabasal)] %in%
                                                         rownames(deg_list_ssc_gerd_p_sig$superficial)))

deg_ssc_gerd_d <- c('Basal'=sum(!rownames(deg_list_ssc_gerd_d_sig$basal)
                                [!rownames(deg_list_ssc_gerd_d_sig$basal) %in% 
                                    rownames(deg_list_ssc_gerd_d_sig$suprabasal)]
                                %in% rownames(deg_list_ssc_gerd_d_sig$superficial)),
                    'Suprabasal'=sum(!rownames(deg_list_ssc_gerd_d_sig$suprabasal)
                                     [!rownames(deg_list_ssc_gerd_d_sig$suprabasal) %in% 
                                         rownames(deg_list_ssc_gerd_d_sig$basal)]
                                     %in% rownames(deg_list_ssc_gerd_d_sig$superficial)),
                    'Superficial'=sum(!rownames(deg_list_ssc_gerd_d_sig$superficial)
                                      [!rownames(deg_list_ssc_gerd_d_sig$superficial) %in% 
                                          rownames(deg_list_ssc_gerd_d_sig$basal)]
                                      %in% rownames(deg_list_ssc_gerd_d_sig$suprabasal)),
                    'Basal&Suprabasal'=sum(rownames(deg_list_ssc_gerd_d_sig$basal)
                                           [!rownames(deg_list_ssc_gerd_d_sig$basal) %in%
                                               rownames(deg_list_ssc_gerd_d_sig$superficial)] %in% 
                                             rownames(deg_list_ssc_gerd_d_sig$suprabasal)),
                    'Suprabasal&Superficial'=sum(rownames(deg_list_ssc_gerd_d_sig$suprabasal)
                                                 [!rownames(deg_list_ssc_gerd_d_sig$suprabasal) %in%
                                                     rownames(deg_list_ssc_gerd_d_sig$basal)] %in% 
                                                   rownames(deg_list_ssc_gerd_d_sig$superficial)),
                    'Basal&Superficial'=sum(rownames(deg_list_ssc_gerd_d_sig$basal)
                                            [!rownames(deg_list_ssc_gerd_d_sig$basal) %in%
                                                rownames(deg_list_ssc_gerd_d_sig$suprabasal)] %in% 
                                              rownames(deg_list_ssc_gerd_d_sig$superficial)),
                    'Basal&Suprabasal&Superficial'=sum(rownames(deg_list_ssc_gerd_d_sig$basal)
                                                       [rownames(deg_list_ssc_gerd_d_sig$basal) %in% 
                                                           rownames(deg_list_ssc_gerd_d_sig$suprabasal)] %in%
                                                         rownames(deg_list_ssc_gerd_d_sig$superficial)))

fit_ssc_gerd_p <- euler(deg_ssc_gerd_p, shape='circle')
fit_ssc_gerd_d <- euler(deg_ssc_gerd_d, shape='circle')

pdf(file=paste0(output_dir, '/epi_super_degVenn_layers_loc_SScGerd_fc',fc_min,'.pdf'), width=6, height=4)
plot_grid(plot(fit_ssc_gerd_p, labels=list(font=4), fills=list(fill=layer_cols_3, alpha=0.5), cex.main=0.8,
               main='SSc vs. GERD DEGs\nProximal', quantities=list(type=c('counts','percent'),
                                                                   cex=0.7, font=3, round=3)),
          plot(fit_ssc_gerd_d, labels=list(font=4), fills=list(fill=layer_cols_3, alpha=0.5), cex.main=0.8,
               main='SSc vs. GERD DEGs\nDistal', quantities=list(type=c('counts','percent'),
                                                                 cex=0.7, font=3, round=3)), ncol=2)
dev.off()



### One plot with one circle for each DEG comparison

ssc_hc_genes_p <- unique(c(rownames(deg_list_ssc_hc_p_sig$basal),
                           rownames(deg_list_ssc_hc_p_sig$suprabasal), 
                           rownames(deg_list_ssc_hc_p_sig$superficial)))
gerd_hc_genes_p <- unique(c(rownames(deg_list_gerd_hc_p_sig$basal),
                            rownames(deg_list_gerd_hc_p_sig$suprabasal), 
                            rownames(deg_list_gerd_hc_p_sig$superficial)))
ssc_gerd_genes_p <- unique(c(rownames(deg_list_ssc_gerd_p_sig$basal),
                             rownames(deg_list_ssc_gerd_p_sig$suprabasal), 
                             rownames(deg_list_ssc_gerd_p_sig$superficial)))
ssc_hc_genes_d <- unique(c(rownames(deg_list_ssc_hc_d_sig$basal),
                           rownames(deg_list_ssc_hc_d_sig$suprabasal), 
                           rownames(deg_list_ssc_hc_d_sig$superficial)))
gerd_hc_genes_d <- unique(c(rownames(deg_list_gerd_hc_d_sig$basal),
                            rownames(deg_list_gerd_hc_d_sig$suprabasal), 
                            rownames(deg_list_gerd_hc_d_sig$superficial)))
ssc_gerd_genes_d <- unique(c(rownames(deg_list_ssc_gerd_d_sig$basal),
                             rownames(deg_list_ssc_gerd_d_sig$suprabasal), 
                             rownames(deg_list_ssc_gerd_d_sig$superficial)))

deg_conds_p <- c('SSc_HC'=sum(!ssc_hc_genes_p[!ssc_hc_genes_p %in% gerd_hc_genes_p] %in% ssc_gerd_genes_p),
                 'GERD_HC'=sum(!gerd_hc_genes_p[!gerd_hc_genes_p %in% ssc_hc_genes_p] %in% ssc_gerd_genes_p),
                 'SSc_GERD'=sum(!ssc_gerd_genes_p[!ssc_gerd_genes_p %in% ssc_hc_genes_p] %in% gerd_hc_genes_p),
                 'SSc_HC&GERD_HC'=sum(ssc_hc_genes_p[!ssc_hc_genes_p %in% ssc_gerd_genes_p] %in% gerd_hc_genes_p),
                 'GERD_HC&SSc_GERD'=sum(gerd_hc_genes_p[!gerd_hc_genes_p %in% ssc_hc_genes_p] %in% ssc_gerd_genes_p),
                 'SSc_HC&SSc_GERD'=sum(ssc_hc_genes_p[!ssc_hc_genes_p %in% gerd_hc_genes_p] %in% ssc_gerd_genes_p),
                 'SSc_HC&GERD_HC&SSc_GERD'=sum(ssc_hc_genes_p[ssc_hc_genes_p %in% gerd_hc_genes_p] %in% ssc_gerd_genes_p))

deg_conds_d <- c('SSc_HC'=sum(!ssc_hc_genes_d[!ssc_hc_genes_d %in% gerd_hc_genes_d] %in% ssc_gerd_genes_d),
                 'GERD_HC'=sum(!gerd_hc_genes_d[!gerd_hc_genes_d %in% ssc_hc_genes_d] %in% ssc_gerd_genes_d),
                 'SSc_GERD'=sum(!ssc_gerd_genes_d[!ssc_gerd_genes_d %in% ssc_hc_genes_d] %in% gerd_hc_genes_d),
                 'SSc_HC&GERD_HC'=sum(ssc_hc_genes_d[!ssc_hc_genes_d %in% ssc_gerd_genes_d] %in% gerd_hc_genes_d),
                 'GERD_HC&SSc_GERD'=sum(gerd_hc_genes_d[!gerd_hc_genes_d %in% ssc_hc_genes_d] %in% ssc_gerd_genes_d),
                 'SSc_HC&SSc_GERD'=sum(ssc_hc_genes_d[!ssc_hc_genes_d %in% gerd_hc_genes_d] %in% ssc_gerd_genes_d),
                 'SSc_HC&GERD_HC&SSc_GERD'=sum(ssc_hc_genes_d[ssc_hc_genes_d %in% gerd_hc_genes_d] %in% ssc_gerd_genes_d))

fit_conds_p <- euler(deg_conds_p, shape='circle')
fit_conds_d <- euler(deg_conds_d, shape='circle')

pdf(file=paste0(output_dir, '/epi_super_degVenn_layers_conds_loc_fc',fc_min,'.pdf'), width=6, height=4)
plot_grid(plot(fit_conds_p, labels=list(font=4), fills=list(fill=c(cond_cols[1:2],'khaki2'), alpha=0.5), cex.main=0.8,
               main='DEGs\nProximal', quantities=list(type=c('counts','percent'),
                                                      cex=0.7, font=3, round=3)),
          plot(fit_conds_d, labels=list(font=4), fills=list(fill=c(cond_cols[1:2],'khaki2'), alpha=0.5), cex.main=0.8,
               main='DEGs\nDistal', quantities=list(type=c('counts','percent'),
                                                    cex=0.7, font=3, round=3)), ncol=2)
dev.off()


#### Plot by diff score bins ####

deg_list_ssc_hc_diffbins <- loadRData(paste0(output_dir, '/epi_degs_ssc_hc_diffbins.RData'))
deg_list_gerd_hc_diffbins <- loadRData(paste0(output_dir, '/epi_degs_gerd_hc_diffbins.RData'))
deg_list_ssc_gerd_diffbins <- loadRData(paste0(output_dir, '/epi_degs_ssc_gerd_diffbins.RData'))

stripchart(lapply(lapply(deg_list_ssc_hc_diffbins, subset, p_val_adj>0.05), '[',,'avg_log2FC'), method='jitter',
           vertical=T, pch=16,col='grey25', cex=0.2, jitter=0.35, ylim=c(-.9,.9))
stripchart(lapply(lapply(deg_list_ssc_hc_diffbins, subset, p_val_adj<0.05), '[',,'avg_log2FC'), method='jitter',
           vertical=T, pch=16,col='red3', cex=0.2, jitter=0.35, add=T)

pdf(file=paste0(output_dir, '/epi_super_diffBins_nDEGs.pdf'), width=9, height=4)
plot(x=1:10, y=c(lapply(deg_list_ssc_hc_diffbins, 
                        function(x) sum(x$p_val_adj<0.05 & abs(x$avg_log2FC)>0)/nrow(x))), xaxt='n',
     xlab='Differentiation Score Decile', type='o', col=cond_cols[1], pch=16, ylab='N DEGs',
     ylim=c(0,0.3), lwd=2, cex.axis=0.8, cex.lab=0.9)
lines(x=1:10, y=c(lapply(deg_list_gerd_hc_diffbins, 
                         function(x) sum(x$p_val_adj<0.05 & abs(x$avg_log2FC)>0)/nrow(x))), xaxt='n',
      xlab='', type='o', col=cond_cols[2], pch=16, ylab='', lwd=2)
lines(x=1:10, y=c(lapply(deg_list_ssc_gerd_diffbins, 
                         function(x) sum(x$p_val_adj<0.05 & abs(x$avg_log2FC)>0)/nrow(x))), xaxt='n',
      xlab='', type='o', col='khaki2', pch=16, ylab='', lwd=2)
axis(1, at=c(1:10), labels=c('[0-10%]','(10-20%]','(20-30%]','(30-40%]','(40-50%]','(50-60%]',
                             '(60-70%]','(70-80%]','(80-90%]','(90-100%]'), cex.axis=0.8)
legend('topleft',legend=c('SSc vs. HC','GERD vs. HC','SSc vs. GERD'), col=c(cond_cols[1:2],'khaki2'),
       cex=0.8, bty='n', pch=16,pt.cex=1.1)
dev.off()

pdf(file=paste0(output_dir, '/epi_super_UMAP_diffBins_80-90.pdf'), width=6, height=5)
DimPlot(super_set, group.by='diff.bin', cols=c(rep('khaki',8),'navy','khaki'))
dev.off()

# By proximal and distal

deg_list_ssc_hc_diffbins_p <- loadRData(paste0(output_dir, '/epi_degs_ssc_hc_diffbins_Proximal.RData'))
deg_list_gerd_hc_diffbins_p <- loadRData(paste0(output_dir, '/epi_degs_gerd_hc_diffbins_Proximal.RData'))
deg_list_ssc_gerd_diffbins_p <- loadRData(paste0(output_dir, '/epi_degs_ssc_gerd_diffbins_Proximal.RData'))
deg_list_ssc_hc_diffbins_d <- loadRData(paste0(output_dir, '/epi_degs_ssc_hc_diffbins_Distal.RData'))
deg_list_gerd_hc_diffbins_d <- loadRData(paste0(output_dir, '/epi_degs_gerd_hc_diffbins_Distal.RData'))
deg_list_ssc_gerd_diffbins_d <- loadRData(paste0(output_dir, '/epi_degs_ssc_gerd_diffbins_Distal.RData'))

deg_list_gerd_hc_diffbins_sig_p <- lapply(deg_list_gerd_hc_diffbins_p, subset, abs(avg_log2FC)>fc_min & p_val_adj<0.05)

Idents(super_set) <- 'condition'
FeatureScatter(super_set, feature1='diff.score',feature2='S100A9', cols=cond_cols)

super_set$diff.bin <- cut(super_set$diff.score, 
                        breaks=unique(quantile(super_set$diff.score,
                                               probs=seq.int(0,1,by=1/10))), include.lowest=T)

pdf(file=paste0(output_dir, '/epi_super_vlnPlot_cond_diffbin_s100A9.pdf'), width=9, height=4)
VlnPlot(super_set, features=c('S100A8','S100A9'),split.by='condition',group.by='diff.bin', cols=cond_cols,
        stack=T)
dev.off()

FeatureScatter(super_set, feature1='diff.score',feature2='S100A9', shuffle=T, cols=cond_cols, smooth=T)

min_fc <- 0.25

get_prop_deg <- function(deg_df, min_fc) {
  deg_df <- deg_df[(deg_df$pct.1+deg_df$pct.2)>0.01,]
  deg_df$p_val_adj <- deg_df$p_val * nrow(deg_df)
  prop <- sum(deg_df$p_val_adj<0.05 & abs(deg_df$avg_log2FC)>min_fc)/nrow(deg_df)*100
}

pdf(file=paste0(output_dir, '/epi_super_diffBins_pctDEGs_loc_fc0.25.pdf'), width=9, height=9)
par(mfrow=c(2,1), mar=c(3,4,3,2))
plot(x=1:10, y=c(lapply(deg_list_ssc_hc_diffbins_p, 
                        function(x) get_prop_deg(x, min_fc))), xaxt='n',
     xlab='Differentiation Score Decile', type='o', col=cond_cols[1], pch=16, ylab='% DEGs',
     ylim=c(0,5), lwd=2, cex.axis=0.8, cex.lab=0.9, main='Proximal', 
     panel.first={grid(nx=0, ny=NULL, lty=2, col='gray', lwd=1)})
lines(x=1:10, y=c(lapply(deg_list_gerd_hc_diffbins_p, 
                         function(x) get_prop_deg(x, min_fc))), xaxt='n',
      xlab='', type='o', col=cond_cols[2], pch=16, ylab='', lwd=2)
lines(x=1:10, y=c(lapply(deg_list_ssc_gerd_diffbins_p, 
                         function(x) get_prop_deg(x, min_fc))), xaxt='n',
      xlab='', type='o', col='khaki2', pch=16, ylab='', lwd=2)
axis(1, at=c(1:10), labels=c('[0-10%]','(10-20%]','(20-30%]','(30-40%]','(40-50%]','(50-60%]',
                             '(60-70%]','(70-80%]','(80-90%]','(90-100%]'), cex.axis=0.8)
legend('topleft',legend=c('SSc vs. HC','GERD vs. HC','SSc vs. GERD'), col=c(cond_cols[1:2],'khaki2'),
       cex=0.8, bty='o', pch=16,pt.cex=1.1, box.lwd=1, box.col='white',bg='white')
box()
plot(x=1:10, y=c(lapply(deg_list_ssc_hc_diffbins_d, 
                        function(x) get_prop_deg(x, min_fc))), xaxt='n',
     xlab='Differentiation Score Decile', type='o', col=cond_cols[1], pch=16, ylab='% DEGs',
     ylim=c(0,5), lwd=2, cex.axis=0.8, cex.lab=0.9 ,main='Distal', 
     panel.first={grid(nx=0, ny=NULL, lty=2, col='gray', lwd=1)})
lines(x=1:10, y=c(lapply(deg_list_gerd_hc_diffbins_d, 
                         function(x) get_prop_deg(x, min_fc))), xaxt='n',
      xlab='', type='o', col=cond_cols[2], pch=16, ylab='', lwd=2)
lines(x=1:10, y=c(lapply(deg_list_ssc_gerd_diffbins_d, 
                         function(x) get_prop_deg(x, min_fc))), xaxt='n',
      xlab='', type='o', col='khaki2', pch=16, ylab='', lwd=2)
axis(1, at=c(1:10), labels=c('[0-10%]','(10-20%]','(20-30%]','(30-40%]','(40-50%]','(50-60%]',
                             '(60-70%]','(70-80%]','(80-90%]','(90-100%]'), cex.axis=0.8)
legend('topleft',legend=c('SSc vs. HC','GERD vs. HC','SSc vs. GERD'), col=c(cond_cols[1:2],'khaki2'),
       cex=0.8, bty='o', pch=16,pt.cex=1.1, box.lwd=1, box.col='white',bg='white')
box()
dev.off()

#### DEG by layer volcano plots ####




# Pseudobulk

epi_bulk <- AggregateExpression(super_set, assays='')

#### DE by pseudotime analysis ####


# create gene x cell expression matrix
lamian_data <- list()
loc <- 'Distal'

for (loc in c('Proximal', 'Distal')) {
  
  sc_obj <- subset(super_set, location==loc)
  
  lamian_data$expr <- GetAssayData(sc_obj, assay='RNA', slot='data')
  
  # remove genes with low expression (will remove genes with counts>0 in at least X cells)
  cell_counts <- rowSums(lamian_data$expr>0)
  min_cells <- ncol(lamian_data$expr)/100 # 1% of cells
  expr_genes <- names(cell_counts[which(cell_counts>=min_cells)])
  lamian_data$expr <- as.matrix(lamian_data$expr[expr_genes,])
  
  lamian_data$cellanno <- data.frame(Cell=colnames(sc_obj), Sample=sc_obj$orig.ident,
                                     row.names=NULL)
  lamian_data$pseudotime <- (sc_obj$diff.score-min(sc_obj$diff.score)) / 
    (max(sc_obj$diff.score)-min(sc_obj$diff.score))
  
  lamian_data$design <- data.frame(row.names=unique(sc_obj$orig.ident))
  lamian_data$design$intercept <- 1
  lamian_data$design$ssc <- ifelse(substr(rownames(lamian_data$design),1,1)=='P',1,0)
  lamian_data$design$gerd <- ifelse(substr(rownames(lamian_data$design),1,1)=='I',1,0)
  lamian_data$design$location <- ifelse(stringr::str_sub(rownames(lamian_data$design),-1)=='P',
                                        0,1)
  saveh5(expr=lamian_data$expr, pseudotime=lamian_data$pseudotime, 
         cellanno=lamian_data$cellanno, path=paste0(output_dir, 'lamian_data_',loc,'.h5'))
  save(lamian_data, file=paste0(output_dir, '/lamian_data_',loc,'.RData'))
  
  
  #lamian_data <- loadRData(paste0(output_dir, '/lamian_data_',loc,'.RData'))
  
  #lamian_data$expr <- lamian_data$expr[1:500,]
  
  
  pseudo_limits <- quantile(lamian_data$pseudotime,c(0.01,0.99))
  
  lamian_data$pseudotime <- lamian_data$pseudotime[lamian_data$pseudotime>=pseudo_limits[1] & lamian_data$pseudotime<=pseudo_limits[2]]
  
  lamian_data$expr <- lamian_data$expr[,names(lamian_data$pseudotime)]
  lamian_data$cellanno <- lamian_data$cellanno[lamian_data$cellanno$Cell %in% names(lamian_data$pseudotime),]
  
  save(lamian_data, file=paste0(output_dir, '/lamian_data_ns_',loc,'.RData'))
  
  saveh5(expr=lamian_data$expr, pseudotime=lamian_data$pseudotime, 
         cellanno=lamian_data$cellanno, path=paste0(output_dir, 'lamian_data_ns_',loc,'.h5'))
  
  for (n in 1:134) {
    print(n)
    temp_data <- lamian_data
    temp_data$expr <- temp_data$expr[(n*100-99):(n*100),]
    saveh5(expr=temp_data$expr, pseudotime=temp_data$pseudotime, 
           cellanno=temp_data$cellanno, path=paste0(output_dir, 'lamian_data_',loc,'_',n,'.h5'))
    save(temp_data, file=paste0(output_dir, '/lamian_data_',loc,'_',n,'.RData'))
  }
  temp_data <- lamian_data
  temp_data$expr <- temp_data$expr[13301:13334,]
  saveh5(expr=temp_data$expr, pseudotime=temp_data$pseudotime, 
         cellanno=temp_data$cellanno, path=paste0(output_dir, 'lamian_data_',loc,'_134.h5'))
  save(temp_data, file=paste0(output_dir, '/lamian_data_',loc,'_134.RData'))
  
  
  #saveh5(expr=temp_data$expr, pseudotime=temp_data$pseudotime, 
  #       cellanno=temp_data$cellanno, path=paste0(output_dir, 'lamian_data_',loc,'_',n,'.h5'))
  
}

###### Create set of objects for each pairwise comparison

# create gene x cell expression matrix
conds <- c('SSc','GERD','HC')

for (loc in c('Proximal','Distal')) {
  for (cond in conds) {
    lamian_data <- list()
    comp_pair <- conds[!conds %in% cond]
    print(comp_pair)
    epi_subset <- subset(super_set, subset=location==loc)
    epi_subset <- subset(epi_subset, subset=condition %in% comp_pair)
    lamian_data$expr <- GetAssayData(epi_subset, assay='RNA', slot='data')
    
    # remove genes with low expression (will remove genes with counts>0 in at least X cells)
    cell_counts <- rowSums(lamian_data$expr>0)
    min_cells <- ncol(lamian_data$expr)/100 # 1% of cells
    expr_genes <- names(cell_counts[which(cell_counts>=min_cells)])
    lamian_data$expr <- as.matrix(lamian_data$expr[expr_genes,])
    
    lamian_data$cellanno <- data.frame(Cell=colnames(epi_subset), Sample=epi_subset$orig.ident,
                                       row.names=NULL)
    lamian_data$pseudotime <- (epi_subset$diff.score-min(epi_subset$diff.score)) / 
      (max(epi_subset$diff.score)-min(epi_subset$diff.score))
    
    lamian_data$design <- data.frame(row.names=unique(epi_subset$orig.ident))
    lamian_data$design$intercept <- 1
    lamian_data$design$cond <- ifelse(rownames(lamian_data$design) %in% 
                                        names(which(table(epi_subset$condition, epi_subset$orig.ident)[comp_pair[1],]>0)),1,0)
    #lamian_data$design$location <- ifelse(stringr::str_sub(rownames(lamian_data$design),-1)=='P',
    #                                      0,1)
    saveh5(expr=lamian_data$expr, pseudotime=lamian_data$pseudotime, 
           cellanno=lamian_data$cellanno, path=paste0(output_dir, 'lamian_data.',loc,'.',paste(comp_pair,collapse='_'),'.h5'))
    save(lamian_data, file=paste0(output_dir, '/lamian_data.',loc,'.',paste(comp_pair,collapse='_'),'.RData'))
    
    #lamian_data <- loadRData(paste0(output_dir, '/lamian_data.',paste(comp_pair,collapse='_'),'.RData'))
    
    pseudo_limits <- quantile(lamian_data$pseudotime,c(0.01,0.99))
    lamian_data$pseudotime <- lamian_data$pseudotime[lamian_data$pseudotime>=pseudo_limits[1] & lamian_data$pseudotime<=pseudo_limits[2]]
    lamian_data$expr <- lamian_data$expr[,names(lamian_data$pseudotime)]
    lamian_data$cellanno <- lamian_data$cellanno[lamian_data$cellanno$Cell %in% names(lamian_data$pseudotime),]
    
    save(lamian_data, file=paste0(output_dir, '/lamian_data_ns.',loc,'.',paste(comp_pair,collapse='_'),'.RData'))
    
    saveh5(expr=lamian_data$expr, pseudotime=lamian_data$pseudotime, 
           cellanno=lamian_data$cellanno, path=paste0(output_dir, 'lamian_data_ns.',loc,'.',paste(comp_pair,collapse='_'),'.h5'))
    
    for (n in 1:(length(expr_genes)%/%100)) {
      print(n)
      temp_data <- lamian_data
      temp_data$expr <- temp_data$expr[(n*100-99):(n*100),]
      saveh5(expr=temp_data$expr, pseudotime=temp_data$pseudotime, 
             cellanno=temp_data$cellanno, path=paste0(output_dir, 'lamian_data_',n,'.',loc,'.',paste(comp_pair,collapse='_'),'.h5'))
      save(temp_data, file=paste0(output_dir, '/lamian_data_',n,'.',loc,'.',paste(comp_pair,collapse='_'),'.RData'))
    }
    temp_data <- lamian_data
    temp_data$expr <- temp_data$expr[(((length(expr_genes)%/%100)*100)+1):length(expr_genes),]
    saveh5(expr=temp_data$expr, pseudotime=temp_data$pseudotime, 
           cellanno=temp_data$cellanno, path=paste0(output_dir, 'lamian_data_',(length(expr_genes)%/%100)+1,'.',loc,'.',paste(comp_pair,collapse='_'),'.h5'))
    save(temp_data, file=paste0(output_dir, '/lamian_data_',(length(expr_genes)%/%100)+1,'.',loc,'.',paste(comp_pair,collapse='_'),'RData'))
    
    
    #saveh5(expr=temp_data$expr, pseudotime=temp_data$pseudotime, 
    #       cellanno=temp_data$cellanno, path=paste0(output_dir, 'lamian_data_',n,'.h5'))
  }
}


lamian_results <- lamian_test(expr=lamian_data$expr[gene,], cellanno=lamian_data$cellanno,
                              pseudotime=lamian_data$pseudotime, design=lamian_data$design, 
                              verbose.output=T, overall.only=T,
                              test.type='variable',testvar=2, ncores=4, ncores.fit = 4)


lamian_results_chisq <- loadRData(paste0(output_dir, '/lamian_results_chisq.RData'))


stat <- lamian_results_chisq$statistics
stat <- stat[order(stat[,1],-stat[,3]),]

data(expdata)
Res <- lamian_test(expr=expdata$expr, cellanno=expdata$cellanno,
                   pseudotime=expdata$pseudotime, design=expdata$design, 
                   verbose.output=T, overall.only=T,
                   test.type='variable',testvar=2, ncores=4)
stat_exp <- Res$statistics
stat_exp <- stat_exp[order(stat_exp[,1],-stat_exp[,3]),]
diffgene_exp <- rownames(stat_exp)[stat_exp[,1] < 0.05]
Res$populationFit <- getPopulationFit(Res, gene=diffgene_exp, type='variable')

Res$covariateGroupDiff <- getCovariateGroupDiff(Res, gene=diffgene)

Res$cluster <- clusterGene(Res, gene=diffgene, type='variable', k=5)

Res$expr <- expdata$expr

diffgene <- rownames(stat[stat[,grep('^fdr.*overall$', colnames(stat))] < 0.05, ])

lamian_results_chisq$expr <- lamian_data$expr

# convert pseudotime to rank
lamian_results_chisq$pseudotime <- rank(lamian_results_chisq$pseudotime)
lamian_results_chisq$populationFit <- getPopulationFit(lamian_results_chisq, gene=diffgene, 
                                                       type='variable',num.timepoint=1000)

lamian_results_chisq$covariateGroupDiff <- getCovariateGroupDiff(lamian_results_chisq,
                                                                 gene=diffgene,
                                                                 num.timepoint=1000)

lamian_results_chisq$cluster <- clusterGene(lamian_results_chisq, gene=diffgene, 
                                            type='variable')

colnames(lamian_results_chisq$populationFit[[1]]) <- 
  colnames(lamian_results_chisq$populationFit[[2]]) <- colnames(lamian_results_chisq$expr)

plotXDEHm(lamian_results_chisq, cellWidthTotal=180, cellHeightTotal=350, subsampleCell=F,
          sep=':.*')

plotClusterMeanAndDiff(lamian_results_chisq)

saveh5(expr=lamian_data$expr, pseudotime=lamian_data$pseudotime, 
       cellanno=lamian_data$cellanno, path=paste0(output_dir, 'lamian_data.h5'))


expdata_1 <- expdata
expdata_2 <- expdata

expdata_1$expr <- expdata_1$expr[1:60,]
expdata_2$expr <- expdata_2$expr[1:2,]
saveh5(expr=expdata_1$expr, pseudotime=expdata_1$pseudotime, 
       cellanno=expdata_1$cellanno, path=paste0(output_dir, 'expdata_1.h5'))
saveh5(expr=expdata_2$expr, pseudotime=expdata_2$pseudotime, 
       cellanno=expdata_2$cellanno, path=paste0(output_dir, 'expdata_2.h5'))

save(expdata, file=paste0(output_dir, '/lamian_data.RData'))


expdata <- loadRData(paste0(output_dir, '/expdata.RData'))

s <- Sys.time()
res_2 <- lamian_test(expr=expdata_2$expr, cellanno=expdata_2$cellanno,
                     pseudotime=expdata_2$pseudotime, 
                     design=expdata_2$design, verbose.output=T, overall.only=T,
                     test.type='variable',testvar=2,permuiter=100, ncores=4)
f <- Sys.time()
t <- f-s
t


s <- Sys.time()
res_chi <- lamian_test(expr=expdata_1$expr, cellanno=expdata_1$cellanno,
                       pseudotime=expdata_1$pseudotime, ncores.fit=4,
                       design=expdata_1$design, verbose.output=T, test.method='chisq',
                       test.type='variable',testvar=2,ncores=4)
f <- Sys.time()
t <- f-s
t



s <- Sys.time()
res_1 <- lamian_test_h5(expr=paste0(output_dir, 'expdata_1.h5'), cellanno=expdata_1$cellanno,
                        pseudotime=expdata_1$pseudotime, 
                        design=expdata_1$design, verbose.output=T,
                        test.type='variable',testvar=2,permuiter=100, ncores=4)
#res_2 <- lamian_test_h5(expr=paste0(output_dir, 'expdata_2.h5'), cellanno=expdata_2$cellanno,
#                      pseudotime=expdata_2$pseudotime, 
#                      design=expdata_2$design, verbose.output=T,
#                      test.type='variable',testvar=2,permuiter=100, ncores=4)
f <- Sys.time()
t <- f-s
t



##### PseudotimeDE analysis


diff_score_all <- cbind(cell=rownames(super_set@meta.data), differentation.score=super_set$diff.score)

cells <- 1:ncol(super_set)
n_sample <- length(cells)*.2
sub_indices <- future_lapply(seq_len(200), function(x) {
  sample(x=cells, size=0.2*n_sample, replace=F)
}, future.seed=T)

basal_markers <- c('PDPN', 'IGFBP3','COL17A1','KRT15','DST')
suprabasal_markers <- c('KRT4','KRT13','SERPINB3','DSC2','IVL')
superficial_markers <- c('KRT17','KRT78','FLG','CNFN','CRCT1')

sub_scores <- list()

for (i in 1:length(sub_indices)) {
  print(i)
  sc_obj <- super_set[, sub_indices[[i]]]
  sc_obj <- AddModuleScore(sc_obj, features=list(basal_markers), name='basal', assay='RNA')
  sc_obj <- AddModuleScore(sc_obj, features=list(suprabasal_markers), name='suprabasal', assay='RNA')
  sc_obj <- AddModuleScore(sc_obj, features=list(superficial_markers), name='superficial', assay='RNA')
  sc_obj$diff.score <- sc_obj$superficial1+sc_obj$suprabasal1-sc_obj$basal1
  sub_scores[[i]] <- cbind(cell=rownames(sc_obj@meta.data), pseudotime=sc_obj$diff.score)
}

save(sub_scores, file=paste0(output_dir, '/MP_eso_epi_diffScore_subsets.RData'))
sub_scores <- loadRData(paste0(output_dir, '/MP_eso_epi_diffScore_subsets.RData'))

library(tibble)
sub_scores <- lapply(sub_scores, function (x){
  x <- as.data.frame(x)
  names(x) <- c('cell','pseudotime')
  x[,2] <- as.double(x[,2])
  x$pseudotime <- (x$pseudotime-min(x$pseudotime))/
    (max(x$pseudotime)-min(x$pseudotime))
  x <- tibble(x)
  return(x)
})

diff_scores <- as.data.frame(cbind(cell=rownames(super_set@meta.data), pseudotime=super_set$diff.score))
diff_scores$pseudotime <- as.double(diff_scores$pseudotime)
diff_scores$pseudotime <- (diff_scores$pseudotime-min(diff_scores$pseudotime))/
  (max(diff_scores$pseudotime)-min(diff_scores$pseudotime))
diff_scores <- tibble(diff_scores)

expv <- GetAssayData(super_set, slot='counts')

res <- PseudotimeDE::runPseudotimeDE(gene.vec=c('S100A9','STAT4'), ori.tbl=diff_scores, sub.tbl=sub_scores[1:100],
                                     mat=expv, model='nb',mc.cores=4)
print(res)



PseudotimeDE::plotCurve(gene.vec=res$gene, ori.tbl=diff_scores, mat=expv.mat, model.fit=res$gam.fit)

#PseudotimeDE::plotUncertainty(diff_scores, sub_scores[1:3])

res <- PseudotimeDE::runPseudotimeDE(gene.vec=c('CXCL10'), ori.tbl=LPS_ori_tbl, sub.tbl=LPS_sub_tbl[1:100],
                                     mat=LPS_sce, model='nb',mc.cores=4)
print(res)

### Hierarchical clustering



##### MONOCLE ANALYSIS
epi_cds <- as.cell_data_set(super_set)

seurat.partitions <- c(rep(1,length(epi_cds@colData@rownames)))
names(seurat.partitions) <- epi_cds@colData@rownames
epi_cds@clusters@listData[['UMAP']][['partitions']] <- as.factor(seurat.partitions)
epi_cds@clusters@listData[['UMAP']][['clusters']] <- super_set@active.ident

epi_cds@int_colData@listData[['reducedDims']]@listData[['UMAP']] <- super_set@reductions$umap@cell.embeddings

epi_cds <- learn_graph(epi_cds, use_partition=F, verbose=T,
                       learn_graph_control=list(
                         minimal_branch_len=100,
                         nn.k=150,
                         nn.cores=4
                       ))

jpeg(file=paste0(output_dir, '/epi_super_monocleTrajectory.jpg'), width=7, height=8, units='in',res=150)
plot_cells(epi_cds, color_cells_by='cluster', label_groups_by_cluster=F, label_branch_points=T,
           label_roots=T, label_leaves=F,group_label_size=4)
dev.off()
