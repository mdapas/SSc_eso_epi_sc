####  LOAD LIBRARIES  ##############################################################################

.libPaths('~/local/R_libs/')    # Set path using .libPaths function

#Load Libraries
dyn.load('/home/mdn578/local/libs/usr/lib64/libfftw3.so.3')

#BiocManager::install("glmGamPoi")
library(methods)
#.libPaths("/home/mhc0155/R/x86_64-pc-linux-gnu-library/4.1")
library(metap)
library(Seurat)
library(sctransform)
library(ggplot2)
library(glmGamPoi)
library(cowplot)
library(plyr)
library(RColorBrewer)
library(future)
library(scCustomize)

####  SET PARAMETERS  ##############################################################################

set.seed(0123456789)
plan('multisession', workers=2)

##Memory Options
Sys.setenv('R_MAX_VSIZE'=192000000000)
options(future.globals.maxSize = 192000 * 1024^2) # NOTE Calculated for 160 Gb RAM request
maxMem <- Sys.getenv('R_MAX_VSIZE')


scl_weak <- c('P0449-P','P0449-D','P0483-P','P0483-D','P0523-P',
              'P0523-D','P0542-P','P0542-D','P0541-P','P0541-D')

scl_absent <- c('P0491-P','P0491-D','P0628-P','P0628-D','P0630-P','P0630-D','P0656-P','P0656-D')

scl_noStage <- 'P0080-D'

hcs <- c('HCE40-D','HCE40-P','HCE043-D', 'HCE043-P','HCE047-D','HCE047-P',
         'HCE048-P','HCE048-D','HCE049-P','HCE049-D','HCE051-D','HCE051-P')

gerds <- c('IW0065-D','IW0065-P','IW0071-D','IW0071-P','IW0076-D','IW0076-P','IW0077-D','IW0077-P')

#ref_samples <- c('HCE051-P','HCE051-D','HCE047-P','HCE047-D', 'HCE048-P', 'HCE048-D')

samples <- c(scl_weak, scl_absent, scl_noStage, hcs, gerds)

cond='condition'; loc='location'; 

root_dir <- '/projects/p31763/users/mdapas/MP_SSc_scRNAseq/Matrices/'

output_dir <- '/projects/p31763/users/mdapas/MP_SSc_scRNAseq/dapas_analysis'


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

####  MAIN  ##############################################################################


sc_obj <- loadRData(paste0(output_dir, '/All_integratedObj.RData'))
#sc_obj <- loadRData(paste0(output_dir,'/../eso_integratedRefObj_2.RData'))

DefaultAssay(sc_obj) <- 'RNA'
Idents(sc_obj) <- 'types'
ct_markers <- FindAllMarkers(sc_obj, logfc.threshold=0.1, min.pct=0.01, only.pos=T)

write.table(ct_markers, paste0(output_dir, '/eso_ct_markers.tsv'), 
            quote=F, row.names=T, col.names=T, sep='\t')


# Get cluster counts by person 
eso_cluster_df <- data.frame(matrix(ncol=7, nrow=length(samples)))
colnames(eso_cluster_df) <- c('Sample','Location','Condition','CellType','nCells','totCells','Proportion')
i=1
for (sample in samples) {
  dat <- sc_obj@meta.data[sc_obj@meta.data$orig.ident==sample,]
  loc <- ifelse(substr(sample,nchar(sample),nchar(sample))=='P','Proximal','Distal')
  c <- dat[1, cond]
  total <- nrow(dat)
  for (l in names(table(dat$types))) {
    nCells <- nrow(dat[dat$types==l,])
    eso_cluster_df[i,] <- c(sample, loc, as.character(c), l, nCells, total,
                            round(nCells/total,4))
    i <- i+1
  }
}

write.table(eso_cluster_df, paste0(output_dir, '/eso_cellType_props.tsv'), 
            quote=F, row.names=F, col.names=T, sep='\t')



pdf(file=paste0(output_dir, '/PC_plot_integratedAssay.pdf'), width=5, height=6)
ElbowPlot(sc_obj, ndims=50)
dev.off()

#sc_obj <- RunUMAP(object=sc_obj, assay='integrated', dims = 1:40, n.neighbors=200L, min.dist=0.3,
#                  spread=1, repulsion.strength=0.9, negative.sample.rate = 10L, n.epochs=500)


pdf(file=paste0(output_dir, '/All_featurePlot_percent.rb.pdf'), width=7, height=6)
print(FeaturePlot(sc_obj, features='percent.rb', order=T) & 
        scale_colour_gradientn(colours=brewer.pal(n=9,'YlGnBu')[3:9]))
dev.off()

pdf(file=paste0(output_dir, '/All_featurePlot_percent.mt.pdf'), width=7, height=6)
print(FeaturePlot(sc_obj, features='percent.mt', order=T) & 
        scale_colour_gradientn(colours=brewer.pal(n=9,'YlGnBu')[3:9]))
dev.off()

pdf(file=paste0(output_dir, '/All_featurePlot_nFeature_RNA.pdf'), width=7, height=6)
print(FeaturePlot(sc_obj, features='nFeature_RNA', order=T) & 
        scale_colour_gradientn(colours=brewer.pal(n=9,'YlGnBu')[3:9]))
dev.off()


## Plot major cell type markers
markers <- c('KRT6A','PTPRC','CD3D','LYZ',  # epi, CD45, T cells, lymphoid,
             'PECAM1','PDGFRA','TAGLN','TPSAB1') # endo, fibro, macrophage, mast
DefaultAssay(sc_obj) <- 'RNA'
plan('default')
sc_obj <- NormalizeData(sc_obj)
pdf(file=paste0(output_dir, '/All_UMAP_markers.pdf'), width=14, height=6)
print(FeaturePlot(sc_obj, features=markers, order=T, ncol=4) & 
        scale_colour_gradientn(colours=brewer.pal(n=9,'YlGnBu')[3:9]))
dev.off()

## Plot distribution of each individual sample within integrated UMAP
samples <- names(table(sc_obj@meta.data$orig.ident))
for (sample in samples) {
  sc_obj@meta.data$sampleBin <- ifelse(sc_obj@meta.data$orig.ident==sample,1,0)
  pdf(file=paste0(output_dir, '/All_UMAP_SampleDistributions/All_dimPlot_',sample,'.pdf'), width=5, height=5)
  print(DimPlot(sc_obj, reduction='umap', label=F, group.by='sampleBin', order=T,
          cols=c('grey80','black')) + ggtitle(sample) + theme(legend.position = "none"))
  dev.off()
}

## Plot PCs on integrated UMAP
for (i in c(0,1,2,3,4)) {
  pc_start <- 10*i+1
  pc_end <- 10*i+10
  pcs <- paste0('PC_',c(pc_start:pc_end))
  pdf(file=paste0(output_dir, '/All_UMAP_pcs',as.character(pc_start),'_',as.character(pc_end),'.pdf'), width=14, height=6)
  print(FeaturePlot(sc_obj, features=pcs, order=T, ncol=5) & 
          scale_colour_gradientn(colours=brewer.pal(n=11,'Spectral')[1:11])
        & theme(plot.title=element_text(face='bold', size=10), legend.position='none',
                axis.text=element_blank(), axis.title=element_text(size=8), axis.ticks = element_blank(), 
                legend.title=element_text(size=7, face='bold'), legend.text=element_text(size=6.5)))
  dev.off()
}

plan('multisession', workers=4)
DefaultAssay(sc_obj) <- 'integrated'

sc_obj <- FindNeighbors(sc_obj, reduction = "pca", dims = 1:40, k.param=25, n.trees=50)
sc_obj <- FindClusters(sc_obj, resolution=0.25, algorithm=1, n.iter=10, n.start=10)

pdf(file=paste0(output_dir, '/All_UMAP_clusters_k25_res0.25.pdf'), width=7, height=6)
DimPlot(sc_obj, reduction = 'umap', label=T, group.by='integrated_snn_res.0.25',
        order=as.character(c(30:0))) + ggtitle('')
dev.off()


## Plot cell type module heatmaps
module_plots <- list(
  module_plot(sc_obj, 'Epithelial', c('KRT13','KRT15','KRT19','S100A2','S100A8')),
  module_plot(sc_obj, 'Endothelial', c('PECAM1','CDH5','ADGRL4', 'EMCN')),
  module_plot(sc_obj, 'Myeloid', c('LYZ','AIF1','HLA-DRA')),
  module_plot(sc_obj, 'Lymphocytes', c('CD3D','CD3E','CD2','TRBC2')),
  module_plot(sc_obj, 'Mast cells', c('TPSAB1','TPSB2','HPGDS', 'CPA3')),
  module_plot(sc_obj, 'Submucosal glands', c('MUC5B','TFF3','LYZ', 'C6orf58')),
  module_plot(sc_obj, 'Glandular epithelial', c('KRT19', 'KRT14','C19orf33','SLPI')),
  module_plot(sc_obj, 'Fibroblasts', c('PDGFRA','LUM','DCN','COL1A2')),
  module_plot(sc_obj, 'Pericytes', c('KCNJ8','RGS5','PDGFRB', 'HIGD1B')),
  module_plot(sc_obj, 'Smooth muscle cells', c('ACTA2','TAGLN','MYH11')))

pdf(file=paste0(output_dir, '/All_modulePlots.pdf'), width=15, height=6)
print(plot_grid(plotlist=module_plots, ncol=5))
dev.off()


save(sc_obj, file=paste0(output_dir, '/All_integratedObj.RData'))


###################################################################################################################
#### Investigate cluster 12 (with smaller cell type populations)

pdf(file=paste0(output_dir, '/All_UMAP_cl12.pdf'), width=6.2, height=6)
print(DimPlot(sc_obj, reduction = 'umap', label=F, group.by='integrated_snn_res.0.25',
              order=as.character(c(30:0)), cols=c(rep('grey',12),'firebrick1','grey')) + ggtitle('')
      & theme(plot.title=element_text(face='bold', size=10), legend.position='none',
              axis.text=element_blank(), axis.title=element_text(size=8), axis.ticks = element_blank(), 
              legend.title=element_text(size=7, face='bold'), legend.text=element_text(size=6.5)))
dev.off()

#cl12_set <- subset(sc_obj, subset=integrated_snn_res.0.25==12)
cl12_set <- loadRData(paste0(output_dir, '/MP_eso_cl12_integratedObj.RData'))

cl12_set <- RunPCA(object=cl12_set, assay='integrated')
ElbowPlot(cl12_set, ndims=50)

cl12_set <- RunUMAP(object=cl12_set, assay='integrated', dims = 1:25, n.neighbors=15L, min.dist=0.3,
                       spread=1, repulsion.strength=1, negative.sample.rate = 10L, n.epochs=1000)
#DimPlot(myeloid_set, group.by='Phase', label=F)


pdf(file=paste0(output_dir, '/Cl12_dimPlot_condition.pdf'), width=7, height=6)
DimPlot(cl12_set, reduction='umap', label=F, group.by='condition', shuffle=T, 
        cols=c('#FF00FF','#00FFFF','#FFFF00')) + ggtitle('')
dev.off()

pdf(file=paste0(output_dir, '/Cl12_dimPlot_location.pdf'), width=7, height=6)
DimPlot(cl12_set, reduction='umap', label=F, group.by='location', shuffle=T, 
        cols=c('#FF00FF','#00FFFF','#FFFF00')) + ggtitle('')
dev.off()

DimPlot(cl12_set, reduction='umap', label=T, group.by='seurat_clusters')

markers <- c('KRT19','PTPRC','CD3D','LYZ',  # epi, CD45, T cells, lymphoid,
             'PECAM1','PDGFRA','TAGLN','TPSAB1') # endo, fibro, macorphage, mast
markers <- c('C19orf33','TPSAB1','IL32','KLRB1','MUC5B','PTPRC','PDGFRA','PDGFRB','RGS5','ACTA2','MYH11')
markers <- c('LYZ','ZG16B','MUC7','MUC5B','CRISP3','PIGR','C6orf58','DMBT1') # salivary gland cells
markers <- c('MUC5B','TFF3','KRT19','C19orf33','ELF3','EHF')
markers <- c('RGS5','CSPG4','PDGFRB','NOTCH3','HIGD1B','KCNJ8','NDUFA4L2','APOLD1')
markers <- c('PRSS23','THBS1','SFRP4','FNDC1','ACTA2','MYH11','PECAM1','POSTN')
markers <- c('ACTA2','MYH11','MUSTN1', 'TAGLN')
cl12_set@meta.data$hce051d <- ifelse(cl12_set@meta.data$orig.ident=='HCE051-D','HCE051-D','Other')

DefaultAssay(cl12_set) <- 'RNA'
pdf(file=paste0(output_dir, '/Cl_12_UMAP_SMGMarkers.pdf'), width=6, height=7.5)
print(FeaturePlot(cl12_set, features=markers, order=T, ncol=2) & 
        scale_colour_gradientn(colours=brewer.pal(n=9,'YlGnBu')[3:9]) & 
        theme(axis.text=element_blank(), axis.title=element_text(size=9),
              axis.ticks = element_blank()))
dev.off()

#0 - Mast cells
#1 - Glandular epithelial cells
#2 - Fibroblasts
#3 - T cells
#4 - Pericytes
#5 - SMCs
#6 - Submucosal glands

DefaultAssay(cl12_set) <- 'integrated'
cl12_set <- FindNeighbors(cl12_set, reduction = "pca", dims = 1:25, k.param=15, n.trees=5000)
cl12_set <- FindClusters(cl12_set, resolution=0.5, algorithm=1, n.iter=50, n.start=50)


pdf(file=paste0(output_dir, '/Cl_12_UMAP_clusters.pdf'), width=7, height=6)
DimPlot(cl12_set, reduction = 'umap', label=F, group.by='seurat_clusters',
        order=as.character(c(20:0))) + ggtitle('')
dev.off()        

types <- c('Mast cells','Glandular\nepithelial cells','Fibroblasts','T cells',
           'Pericytes', 'SMCs', 'Submucosal glands')

cl12_set@meta.data$types <- mapvalues(cl12_set@meta.data$seurat_clusters, 
                                  from=as.character(c(0:6)), to=types)

pdf(file=paste0(output_dir, '/Cl_12_UMAP_types.pdf'), width=6, height=6)
DimPlot(cl12_set, reduction = 'umap', label=T, group.by='types', 
        cols=c('darkseagreen3','lightpink2','lavenderblush3','steelblue1','khaki','salmon1','honeydew3')) + 
  ggtitle('') + theme(legend.position = "none")
dev.off()

# Cell type module heatmaps
module_plots <- list(
module_plot(cl12_set, 'Mast cells', c('TPSAB1','TPSB2','HPGDS', 'CPA3')),
module_plot(cl12_set, 'T cells', c('CD3D','CD2','IL32','KLRB1')),
module_plot(cl12_set, 'Submucosal glands', c('MUC5B','TFF3','LYZ', 'C6orf58')),
module_plot(cl12_set, 'Glandular epithelium', c('KRT19', 'KRT14','C19orf33','SLPI')),
module_plot(cl12_set, 'Fibroblasts', c('PDGFRA','LUM','DCN','COL1A2')),
module_plot(cl12_set, 'Pericytes', c('KCNJ8','RGS5','PDGFRB', 'HIGD1B')),
module_plot(cl12_set, 'Smooth muscle cells', c('ACTA2','TAGLN','MYH11')))

pdf(file=paste0(output_dir, '/Cl_12_modulePlots.pdf'), width=14, height=6)
print(plot_grid(plotlist=module_plots, ncol=4))
dev.off()


# Examine pericytes vs myofibroblasts
fps_set <- subset(cl12_set, subset=seurat_clusters %in% c(2,4))

fps_set <- RunPCA(object=fps_set, assay='integrated')
ElbowPlot(fps_set, ndims=50)

fps_set <- RunUMAP(object=fps_set, assay='integrated', dims = 1:15, n.neighbors=15L, min.dist=0.1,
                   spread=1, repulsion.strength=1, negative.sample.rate = 10L, n.epochs=5000)
#DimPlot(myeloid_set, group.by='Phase', label=F)

DimPlot(fps_set, reduction='umap', label=T, group.by='condition')

markers <- c('ACTA2','TGFB1','PDGFA','PDGFB', 'FN1', 'COL1A1','PDGFRB','CSPG4')

DefaultAssay(fps_set) <- 'RNA'
pdf(file=paste0(output_dir, '/FP_UMAP_PericyteMarkers.pdf'), width=12, height=6)
print(FeaturePlot(fps_set, features=markers, order=T, ncol=4) & 
        scale_colour_gradientn(colours=brewer.pal(n=9,'YlGnBu')[3:9]) & 
        theme(axis.text=element_blank(), axis.title=element_text(size=9),
              axis.ticks = element_blank()))
dev.off()

pdf(file=paste0(output_dir, '/FP_UMAP_MyoFibroMarkers.pdf'), width=4, height=4)
module_plot(fps_set, 'Myofibroblasts',c('ACTA2','CSPG4','NOTCH1','TGFB1','FN1'))
dev.off()

pdf(file=paste0(output_dir, '/FP_UMAP_PericyteMarkers.pdf'), width=4, height=4)
module_plot(fps_set, 'Pericytes', c('KCNJ8','RGS5','PDGFRB', 'HIGD1B'))
dev.off()


pdf(file=paste0(output_dir, '/All_Vln_ACTA2.pdf'), width=4, height=5)
print(VlnPlot(fps_set, features='ACTA2', group.by='condition', 
              pt.size=0, cols=c('darkseagreen', 'darkseagreen4','lightpink')) & geom_boxplot(width=0.1, fill='white',
                                                                                             outlier.shape=NA) & 
        theme(axis.title.x=element_blank(), legend.position='none')) &
  annotate(geom='text', label=paste0('p=',signif(x[5],2)), x=Inf, y=Inf, hjust=1, vjust=1, colour='darkred')
dev.off()

pdf(file=paste0(output_dir, '/FP_vln_ACTA2.pdf'), width=4, height=4)
VlnPlot(fps_set, features='ACTA2', group.by='condition', 
        pt.size=1) + theme(legend.position='none')
dev.off()

Idents(fps_set) <- 'condition'
FindMarkers(fps_set, ident.1='SSc',ident.2='HC',features='ACTA2', min.pct=0, logfc.threshold=0)



#######################################################################################################
## Annotate sc_obj based on cl_12 results

types <- c('En','Ep','F','GEp','L','Mc','My','P','SMc','SMG')

sc_obj@meta.data$types <- as.character(sc_obj@meta.data$seurat_clusters)

sc_obj@meta.data[sc_obj@meta.data %in% cl12_set@meta.data[cl12_set@meta.data$types=='Mast cells',],'types']

sc_obj@meta.data$types <- mapvalues(sc_obj@meta.data$seurat_clusters, 
                                    from=as.character(c(0:13)), to=types)

sc_obj@meta.data[rownames(cl12_set@meta.data[cl12_set@meta.data$seurat_clusters==0,]),'types'] <- 'Mc'
sc_obj@meta.data[rownames(cl12_set@meta.data[cl12_set@meta.data$seurat_clusters==1,]),'types'] <- 'GEp'
sc_obj@meta.data[rownames(cl12_set@meta.data[cl12_set@meta.data$seurat_clusters==2,]),'types'] <- 'F'
sc_obj@meta.data[rownames(cl12_set@meta.data[cl12_set@meta.data$seurat_clusters==3,]),'types'] <- 'L'
sc_obj@meta.data[rownames(cl12_set@meta.data[cl12_set@meta.data$seurat_clusters==4,]),'types'] <- 'P'
sc_obj@meta.data[rownames(cl12_set@meta.data[cl12_set@meta.data$seurat_clusters==5,]),'types'] <- 'SMCs'
sc_obj@meta.data[rownames(cl12_set@meta.data[cl12_set@meta.data$seurat_clusters==6,]),'types'] <- 'SMG'
sc_obj@meta.data[sc_obj@meta.data$seurat_clusters %in% c(0,1,2,3,4,6,7,8,9,13), 'types'] <- 'Ep'
sc_obj@meta.data[sc_obj@meta.data$seurat_clusters==5, 'types'] <- 'L'
sc_obj@meta.data[sc_obj@meta.data$seurat_clusters==10, 'types'] <- 'My'
sc_obj@meta.data[sc_obj@meta.data$seurat_clusters==11, 'types'] <- 'En'

type_cols <- c('#CD96CD','#8FBC8F','#A52A2A','#698B69','#7AC5CD',
                    '#4682B4','#EEA2AD','#DAA520','#FF7F50','#76EE00')

pdf(file=paste0(output_dir, '/All_dimPlot_annotation_0001.pdf'), width=7, height=6)
print(DimPlot(sc_obj, reduction='umap', label=T, group.by='types', shuffle=T, raster=F,
        cols=type_cols, label.size=3, pt.size=0.0001) + ggtitle('') & 
  theme(axis.text=element_blank(), axis.title=element_text(size=9),
        axis.ticks = element_blank(), legend.text=element_text(size=9)))
dev.off()


pdf(file=paste0(output_dir, '/All_violin_nFeature.pdf'), width=8, height=5)
ggplot(sc_obj@meta.data, aes(x = types, y = nFeature_RNA, fill=types))+
  geom_violin(scale='width') + labs(y= "nFeature_RNA", x= "Cluster") + ggtitle('Esophageal Cell Clusters') + 
  theme_classic() + scale_fill_manual(values=type_cols) +
  theme(axis.title=element_text(face='bold',size=8), axis.ticks = element_blank(), 
        legend.position='none', panel.grid.major.y = element_line(), panel.grid.minor.y=element_line())
dev.off()


pdf(file=paste0(output_dir, '/All_violin_nCount.pdf'), width=8, height=5)
ggplot(sc_obj@meta.data, aes(x = types, y = FLI1, fill=types))+
  geom_violin(scale='width') + labs(y= "nCount_RNA", x= "Cluster") + ggtitle('Esophageal Cell Clusters') + 
  theme_classic() + scale_fill_manual(values=type_cols) +
  theme(axis.title=element_text(face='bold',size=8), axis.ticks = element_blank(), 
        legend.position='none', panel.grid.major.y = element_line(), panel.grid.minor.y=element_line())
dev.off()


# Compare C19orf33 between GEp and Ep
pdf(file=paste0(output_dir, '/All_UMAP_AIRE.pdf'), width=7, height=6)
print(FeaturePlot(sc_obj, features='AIRE', order=T, ncol=1, raster=T) & 
        scale_colour_gradientn(colours=brewer.pal(n=9,'YlGnBu')[3:9]))
dev.off()

DefaultAssay(sc_obj) <- 'RNA'
Idents(sc_obj) <- 'types'
x <- FindMarkers(sc_obj, ident.1='GEp', ident.2='Ep',logfc.threshold=0, min.pct=0, features=c('C19orf33'))

pdf(file=paste0(output_dir, '/All_Vln_FLI1.pdf'), width=6, height=3)
print(VlnPlot(sc_obj, features='FLI1', group.by='types', split.by = 'condition',
              pt.size=0, cols=cond_cols) &
        theme(axis.title.x=element_blank(), legend.position='none')) &
  annotate(geom='text', label=paste0('p=',signif(x[1],2)), x=Inf, y=Inf, hjust=1, vjust=1, colour='darkred')
dev.off()


#### Calculate cell type proportions by condition
prop_df <- data.frame(matrix(ncol=3, nrow=length(table(sc_obj$types))))
colnames(prop_df) <- c('Cell type','Condition','Percent')
i=1
for (cond in c('HC','GERD','SSc')) {
  total <- nrow(sc_obj@meta.data[sc_obj$condition==cond,])
  for (type in c('En','Ep','L','My','Other')) {
    if (type=='Other') {
      prop_df[i,] <- c(type, cond, 
                       round(nrow(sc_obj@meta.data[sc_obj$condition==cond & 
                                                     sc_obj$types %in% c('F','GEp','Mc','P','SMc','SMG'),])/total,5))
    } else {
      prop_df[i,] <- c(type, cond, round(nrow(sc_obj@meta.data[sc_obj$condition==cond & sc_obj$types==type,])/total,5))
    }
    i <- i+1
  }
}

prop_df$Percent <- as.numeric(prop_df$Percent)*100
prop_df$Condition <- factor(prop_df$Condition, levels=c('HC','GERD','SSc'))

pdf(file=paste0(output_dir, '/All_cellType_vs_condition.pdf'), width=7, height=8)
print(ggplot(prop_df, aes(x=Condition, y=Percent, fill=`Cell type`)) +
        geom_bar(stat='identity') + theme_minimal() +
        geom_text(aes(label=paste(round(Percent,1),'%')), position=position_stack(vjust=0.5)) +
        scale_fill_manual('Cell type', values=c('plum3','darkseagreen','cadetblue3','lightpink2','grey')))
dev.off()



type_df <- data.frame(table(sc_obj@meta.data$types)/nrow(sc_obj@meta.data))
names(type_df) <- c('type','prop')
type_df$percent <- paste0(round(type_df$prop*100,2),'%')
type_df$cols <- type_cols
type_df <- type_df[order(type_df$prop,decreasing=T),]
rownames(type_df) <- 1:nrow(type_df)
type_df$type <- factor(type_df$type, level=type_df$type[1:10])

pdf(file=paste0(output_dir, '/All_types_pie.pdf'), width=6, height=5)
ggplot(type_df, aes(x="", y=prop, fill=type)) + geom_bar(stat='identity',width=1) +
  coord_polar(theta='y') + scale_fill_manual(name='Cell type',
                    labels=paste0(type_df$type,' (',type_df$percent,')'),values=type_df$cols) +
  theme_void()
dev.off()

## Plot stacked violin of markers by cell type
markers <- c('KRT6A','KRT14','MUC5B','PTPRC','CD3D','LYZ','HLA-DRA','TPSAB1','PECAM1',
             'PDGFRA','PDGFRB','RGS5','TAGLN','ACTA2') # endo, fibro, macrophage, mast
#DefaultAssay(sc_obj) <- 'RNA'
markers <- c('KRT6A','CLU','C19orf33','DEFB1','MMP7','KRT7','KRT14','KRT78')
reorder <- c(1,6,10,2,3,5,4,7,8,9)

sc_obj@meta.data$types <- factor(sc_obj@meta.data$types, levels=type_df$type[reorder])

pdf(file=paste0(output_dir, '/All_types_markers_violin.pdf'), width=5, height=10)
print(Stacked_VlnPlot(sc_obj, markers, pt.size=0, group.by='types', raster=F, plot_spacing=0,
                      colors_use=type_df$cols[reorder], x_lab_rotate=T) &
  theme(axis.line.y = element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()))
dev.off()
  
Idents(sc_obj) <- sc_obj@meta.data$types
GEp <- FindMarkers(sc_obj,ident.1='GEp',assay='RNA')

# explore markers
sc_obj@meta.data$DEG_group <- ifelse(sc_obj@meta.data$types=='Ep',1,
                                     ifelse(sc_obj@meta.data$types=='GEp',2,NA))
Idents(sc_obj) <- sc_obj@meta.data$DEG_group
ep <- FindMarkers(sc_obj, ident.1=1, ident.2=2)

#############################################################################################
# filter for cd45 cells

pdf(file=paste0(output_dir, '/All_UMAP_cd45.pdf'), width=6.2, height=6)
print(DimPlot(sc_obj, reduction = 'umap', label=F, group.by='types',
              order=as.character(c(30:0)), cols=c(rep('grey',2),'firebrick1',rep('grey',3),'firebrick1','grey','firebrick1','grey')) + ggtitle('')
      & theme(plot.title=element_text(face='bold', size=10), legend.position='none',
              axis.text=element_blank(), axis.title=element_text(size=8), axis.ticks = element_blank(), 
              legend.title=element_text(size=7, face='bold'), legend.text=element_text(size=6.5)))
dev.off()


#cd45_set <- subset(sc_obj, subset=types %in% c('L','My','Mc'))
cd45_set <- loadRData(paste0(output_dir, '/MP_eso_CD45_integratedObj.RData'))
#

#########  split by sample for integration ####
cd45_list <- SplitObject(cd45_set, split.by='orig.ident')

#Set default assay
for (i in 1:length(cd45_list)) { 
  cd45_list[[i]] <- subset(cd45_list[[i]])
  cd45_list[[i]] <- SCTransform(cd45_list[[i]], vst.flavor="v2", verbose=T)
  DefaultAssay(cd45_list[[i]]) <- 'SCT'
}


#Select integration features + Prep integration
int_features <- SelectIntegrationFeatures(cd45_list, nfeatures = 2000)
cd45_list <- PrepSCTIntegration(cd45_list, anchor.features=int_features)

table(cd45_set$orig.ident)

#Perform integration
ref_anchors <- FindIntegrationAnchors(object.list=cd45_list, normalization.method='SCT', dims=10,
                                      anchor.features=int_features, k.filter=10, k.anchor=10, k.score=10)

#save(ref_anchors, file=paste0(output_dir,'RefObj_anchors_myeloid.RData'))
cd45_set <- IntegrateData(anchorset=ref_anchors, normalization.method = 'SCT', k.weight=30)

#########  Plot CD45+ cells ####


cd45_set <- RunPCA(object=cd45_set, assay='integrated')
ElbowPlot(cd45_set, ndims=50)

cd45_set <- RunUMAP(object=cd45_set, assay='integrated', dims = 1:20, n.neighbors=25L, min.dist=0.25,
                    spread=1, repulsion.strength=1, negative.sample.rate = 10L, n.epochs=1000)
#DimPlot(myeloid_set, group.by='Phase', label=F)

DimPlot(cd45_set, reduction='umap', label=T, group.by='types')

pdf(file=paste0(output_dir, '/CD45_dimPlot_annotation.pdf'), width=6, height=3)
print(DimPlot(cd45_set, reduction='umap', label=T, group.by='types', shuffle=T, raster=F,
              cols=c('#7AC5CD','#4682B4','#EEA2AD'),
              label.size=3) + ggtitle('') & 
        theme(axis.text=element_blank(), axis.title=element_text(size=9),
              axis.ticks = element_blank(), legend.text=element_text(size=9)))
dev.off()


markers <- c('KRT19','PTPRC','CD3D','LYZ',  # epi, CD45, T cells, lymphoid,
             'PECAM1','PDGFRA','TAGLN','TPSAB1') # endo, fibro, macorphage, mast
# b cells, macrophage, mast cells, epi, tregs

markers <- c('LYZ','CD1A','CD14','percent.rb','G2M.Score',
             'LAMP3','CLEC9A','FCN1','CD19','KRT19') 

markers <- c('TREM2','CD14','TIMD4','MERTK','CD68',
             'CD163','CCR2','ITGAX','CD19','CD1A') 

markers <- c('FOXP3','IL17A','CD8A','NKG7','PDGFRA','TAGLN','TPSAB1','FABP5','TXN','percent.rb') 
markers <- c('PTPRC','TPSAB1','CD3D','LYZ')
  
DefaultAssay(cd45_set) <- 'RNA'
pdf(file=paste0(output_dir, '/CD45_UMAP_markers.pdf'), width=6, height=6)
print(FeaturePlot(cd45_set, features=markers, order=T, ncol=2) & 
        scale_colour_gradientn(colours=brewer.pal(n=9,'YlGnBu')[3:9]) & 
        theme(axis.text=element_blank(), axis.title=element_text(size=9),
              axis.ticks = element_blank(), legend.text=element_text(size=9)))
dev.off()


DefaultAssay(cd45_set) <- 'integrated'
cd45_set <- FindNeighbors(cd45_set, reduction = "pca", dims = 1:30, k.param=27, n.trees=1000)
cd45_set <- FindClusters(cd45_set, resolution=0.8, algorithm=1, n.iter=30, n.start=30)

pdf(file=paste0(output_dir, '/CD45_UMAP_condition.pdf'), width=7, height=6)
DimPlot(cd45_set, reduction = 'umap', label=F, group.by='condition',
        order=as.character(c(20:0))) + ggtitle('')
dev.off()        

save(cd45_set, file=paste0(output_dir, '/MP_eso_CD45_integratedObj.RData'))



