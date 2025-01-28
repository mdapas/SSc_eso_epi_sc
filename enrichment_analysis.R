
.libPaths('~/local/R_libs/')    # Set path using .libPaths function

library(fgsea)
library(dplyr)
library(ggplot2)
library(msigdbr)
library(clusterProfiler)
library(stringr)
library(RColorBrewer)
library(VIM)


#Load Libraries
dyn.load('/home/mdn578/local/libs/usr/lib64/libfftw3.so.3')

####  SET PARAMETERS  ##############################################################################

set.seed(0123456789)
plan('multisession', workers=4)

##Memory Options
Sys.setenv('R_MAX_VSIZE'=192000000000)
options(future.globals.maxSize = 192000 * 1024^2) # NOTE Calculated for 160 Gb RAM request
maxMem <- Sys.getenv('R_MAX_VSIZE')

cond='condition'; loc='location'; 

cond_cols <- c('#8156B3','#e26b53','#A8A39d')
layer_cols_3 <- c('#4ca5b1','#e9b85d','#c72f4c')
layer_cols_5 <- brewer.pal(n=11,'Spectral')[c(10,9,7,4,2)]

root_dir <- '/projects/p31763/users/mdapas/MP_SSc_scRNAseq/Matrices/'

output_dir <- '/projects/p31763/users/mdapas/MP_SSc_scRNAseq/dapas_analysis/epi_analysis/'


####  FUNCTIONS  ##############################################################################


#Load Options
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}



deg_list_ssc_hc <- loadRData(paste0(output_dir, '/epi_degs_ssc_hc_layers.RData'))
deg_list_gerd_hc <- loadRData(paste0(output_dir, '/epi_degs_gerd_hc_layers.RData'))
deg_list_ssc_gerd <- loadRData(paste0(output_dir, '/epi_degs_ssc_gerd_layers.RData'))

deg_list_ssc_hc_p <- loadRData(paste0(output_dir, '/epi_degs_ssc_hc_layersRep_Proximal.RData'))
deg_list_gerd_hc_p <- loadRData(paste0(output_dir, '/epi_degs_gerd_hc_layersRep_Proximal.RData'))
deg_list_ssc_gerd_p <- loadRData(paste0(output_dir, '/epi_degs_ssc_gerd_layersRep_Proximal.RData'))
deg_list_ssc_hc_d <- loadRData(paste0(output_dir, '/epi_degs_ssc_hc_layersRep_Distal.RData'))
deg_list_gerd_hc_d <- loadRData(paste0(output_dir, '/epi_degs_gerd_hc_layersRep_Distal.RData'))
deg_list_ssc_gerd_d <- loadRData(paste0(output_dir, '/epi_degs_ssc_gerd_layersRep_Distal.RData'))



deg_list_ssc_hc_diffbins <- loadRData(paste0(output_dir, '/epi_degs_ssc_hc_diffbins.RData'))
deg_list_gerd_hc_diffbins <- loadRData(paste0(output_dir, '/epi_degs_gerd_hc_diffbins.RData'))
deg_list_ssc_gerd_diffbins <- loadRData(paste0(output_dir, '/epi_degs_ssc_gerd_diffbins.RData'))

deg_list_ssc_hc_diffbins_p <- loadRData(paste0(output_dir, '/epi_degs_ssc_hc_diffbins_Proximal.RData'))
deg_list_gerd_hc_diffbins_p <- loadRData(paste0(output_dir, '/epi_degs_gerd_hc_diffbins_Proximal.RData'))
deg_list_ssc_gerd_diffbins_p <- loadRData(paste0(output_dir, '/epi_degs_ssc_gerd_diffbins_Proximal.RData'))
deg_list_ssc_hc_diffbins_d <- loadRData(paste0(output_dir, '/epi_degs_ssc_hc_diffbins_Distal.RData'))
deg_list_gerd_hc_diffbins_d <- loadRData(paste0(output_dir, '/epi_degs_gerd_hc_diffbins_Distal.RData'))
deg_list_ssc_gerd_diffbins_d <- loadRData(paste0(output_dir, '/epi_degs_ssc_gerd_diffbins_Distal.RData'))


c2_df <- msigdbr(species='Homo sapiens', category='C2')
c2_df <- c2_df[grepl('^CP', c2_df$gs_subcat),]
c7_df <- msigdbr(species='Homo sapiens', category='C7')

fgsea_c2_sets<- c2_df %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_c7_sets<- c7_df %>% split(x = .$gene_symbol, f = .$gs_name)

get_genes <- function(df) {
  #df <- df[abs(df$avg_log2FC)<0.1,]
  #df <- df[df$pct.1+df$pct.2>0,]
  df <- df[df$pct.1>0.01 & df$pct.2>0.01,]
  df <- df[!is.na(df$p_val),]
  df$p.diff <- df$pct.1-df$pct.2
  df[df$p_val==0,'p_val'] <- 5e-324
  df$weight <- -log10(df$p_val) * abs(df$avg_log2FC) #(df$p.diff)^2 * 
  df <- df[order(df$weight, decreasing=T),]
  v <- df[,'weight']
  names(v) <- rownames(df)
  v <- v[!is.na(v)]
  v <- v[!grepl('^MT-', names(v))]
  v <- v[!grepl('^RP', names(v))]
  v <- v[!grepl('^MRP', names(v))]
  return(v)
}

enrich_ssc_hc <- list(); enrich_gerd_hc <- list()
enrich_ssc_hc_p <- list(); enrich_gerd_hc_p <- list()
enrich_ssc_hc_d <- list(); enrich_gerd_hc_d <- list()

layers <- names(deg_list_ssc_hc)

# not the most efficient code
for (layer in layers_rep) {
  print(layer)
  ssc_hc_genes <- get_genes(deg_list_ssc_hc_p[[layer]])
  enrich_ssc_hc_p[[layer]] <- fgseaMultilevel(fgsea_c2_sets, stats=ssc_hc_genes,
                                            minSize=10, scoreType='pos', nPermSimple=10000)
  enrich_ssc_hc_p[[layer]]$cond <- 'SSc'
  enrich_ssc_hc_p[[layer]]$loc <- 'Proximal'
  enrich_ssc_hc_p[[layer]]$layer <- layer
  
  gerd_hc_genes <- get_genes(deg_list_gerd_hc_p[[layer]])
  enrich_gerd_hc_p[[layer]] <- fgseaMultilevel(fgsea_c2_sets, stats=gerd_hc_genes, 
                                             minSize=10, scoreType='pos', nPermSimple=10000)
  enrich_gerd_hc_p[[layer]]$cond <- 'GERD'
  enrich_gerd_hc_p[[layer]]$loc <- 'Proximal'
  enrich_gerd_hc_p[[layer]]$layer <- layer
  
  ssc_hc_genes <- get_genes(deg_list_ssc_hc_d[[layer]])
  enrich_ssc_hc_d[[layer]] <- fgseaMultilevel(fgsea_c2_sets, stats=ssc_hc_genes,
                                              minSize=10, scoreType='pos', nPermSimple=10000)
  enrich_ssc_hc_d[[layer]]$cond <- 'SSc'
  enrich_ssc_hc_d[[layer]]$loc <- 'Distal'
  enrich_ssc_hc_d[[layer]]$layer <- layer
  
  gerd_hc_genes <- get_genes(deg_list_gerd_hc_d[[layer]])
  enrich_gerd_hc_d[[layer]] <- fgseaMultilevel(fgsea_c2_sets, stats=gerd_hc_genes, 
                                               minSize=10, scoreType='pos', nPermSimple=10000)
  enrich_gerd_hc_d[[layer]]$cond <- 'GERD'
  enrich_gerd_hc_d[[layer]]$loc <- 'Distal'
  enrich_gerd_hc_d[[layer]]$layer <- layer
}

enrich_plot_df <- data.frame()
top_pathways <- c()
for (i in 1:3) {
  for (df in list(enrich_ssc_hc_p[[i]],enrich_ssc_hc_d[[i]],enrich_gerd_hc_p[[i]],enrich_gerd_hc_d[[i]])) {
    enrich_plot_df <- rbind(enrich_plot_df, head(df[order(df$pval),], n=10))
    top_pathways <- c(top_pathways, head(df[order(df$pval),], n=10)$pathway)
  }
}


enrich_plot_df$db <- sub("\\_.*", "",enrich_plot_df$pathway)
enrich_plot_df$pathway <- str_extract(enrich_plot_df$pathway, '_.*')
enrich_plot_df$pathway <- gsub("_"," ", substr(enrich_plot_df$pathway,2,nchar(enrich_plot_df$pathway)))
enrich_plot_df$layer <- factor(enrich_plot_df$layer, levels=c('basal','suprabasal','superficial'))
enrich_plot_df$cond <- factor(enrich_plot_df$cond, levels=c('SSc','GERD'))
enrich_plot_df$loc <- factor(enrich_plot_df$loc, levels=c('Proximal','Distal'))
enrich_plot_df$pathway <- factor(enrich_plot_df$pathway,
                                 levels=names(sort(table(enrich_plot_df$pathway))))

pdf(file=paste0(output_dir, '/epi_DEG.SSc_HC.GERD_HC.GSEA_layer.pdf'), width=7.5, height=10)
ggplot(enrich_plot_df, aes(x=cond, y=pathway, color=NES, size=-log10(pval))) + 
  geom_point() + facet_grid2(layer~loc, space='free', scale='free_y',
                             strip=strip_themed(background_y=elem_list_rect(fill=layer_cols_3),
                             text_x=element_text(face='bold', size=7.5), 
                             text_y=element_text(face='bold', size=8))) + 
  theme(text=element_text(size=7), axis.text.x=element_text(face='bold', size=10, color=cond_cols[1:2]),
        axis.title = element_blank()) + scale_color_gradient(low='skyblue',high='navy')
dev.off()

#n enrich_dfs <- list()
for (layer in names(deg_list)) {
  enrich_dfs[['ssc_hc']] <- data.frame(row.names=names(fgsea_c2_sets))
  enrich_dfs[['gerd_hc']] <- data.frame(row.names=names(fgsea_c2_sets))
  enrich_dfs[['ssc_gerd']] <- data.frame(row.names=names(fgsea_c2_sets))
  
  enrich_dfs[['ssc_hc']] <- merge(enrich_dfs[['ssc_hc']], enrich_ssc_hc[[layer]])
  enrich_dfs[['ssc_hc']][,layer] <- enrich_ssc_hc[[layer]]
}

##### GENE EXPRESSION PLOTS FOR ENRICHED PATHWAY GENES

s100s <- c('S100A7','S100A8','S100A9')
serpins <- c('SERPINB1','SERPINB3','SERPINB4')
sprrs <- c('SPRR1A','SPRR1B','SPRR2E')
keratins <- c('KRT1','KRT6A','KRT14','KRT16','KRT17')

epi_set <- AddModuleScore(epi_set, features=list(s100s), name='S100s', assay='RNA')
epi_set <- AddModuleScore(epi_set, features=list(serpins), name='Serpins', assay='RNA')
epi_set <- AddModuleScore(epi_set, features=list(sprrs), name='SPRRs', assay='RNA')
epi_set <- AddModuleScore(epi_set, features=list(keratins), name='Keratins', assay='RNA')

layers_rep <- c('basal','replicating basal','replicating suprabasal','suprabasal','superficial')
require(vioplot)
for (set in c('S100s1','Serpins1','SPRRs1','Keratins1')) {
    data <- epi_set@meta.data
    data$condition <- factor(data$condition, c('HC','GERD','SSc'))
    range <- range(data[,set])
    ymin <- min(range)-0.1*sum(abs(range)); ymax <- max(range)+0.1*sum(abs(range))
    
    pdf(file=paste0(output_dir, '/epi_',set,'.',loc,'.pdf'), width=5, height=2.5)
    par(mfrow=c(1,5),mar=c(2,0.4,2,0.4))
    for (l in layers_rep) {
      vioplot(as.formula(paste(set,'~condition')), data=data[data$layers_rep==l & data$location=='Proximal',], ylab='', col=rev(cond_cols),
              xlab='', cex.main=1.5, yaxt='n',colMed='black',colMed2=rev(cond_cols),cex=1.5, xaxt='n', rectCol='white', side='left',
              ylim=c(ymin,ymax), panel.first={axis(2,tck=1,col.ticks='gray80', lty=3,labels=T)}, pchMed=23,lwd=1.05,plotCentre='line')
      vioplot(as.formula(paste(set,'~condition')), data=data[data$layers_rep==l & data$location=='Distal',], col=rev(cond_cols),
              cex.main=1.5, yaxt='n',colMed='black',colMed2=rev(cond_cols),cex=1.5, xaxt='n', rectCol='white', side='right',add=T,
              pchMed=23,lwd=1.05,plotCentre='line')
      mtext(str_to_title(l), side=3, padj=-1, las=1 ,cex=0.8, col='grey30')
    }
    dev.off()
}         

leadingEdgeGenes <- c('S100A7','PI3','SERPINB4','MUC21','CLECL2B','CTSC','SLPI','SPRR1B','SPRR2E','IL18',
                      'LCE3D','SPRR2B','KLK12','LCE3E','KRT1','SPRR2A','KLK8','KRT6A','PLCG2','FOS','PLCG2',
                      'SAA1','UBA52','SLPI','LRG1','CTSH','SERPINB3','MUC4','FABP5','S100A8','S100A9','AKR1C3',
                      'CXCL14','SPON2','LGALS7','LGALS3','KLK13','SERPINB1','SPRR2F','ANXA3','AKR1B10')
leadingEdgeGenes <- leadingEdgeGenes[leadingEdgeGenes %in% rownames(epi_set)]

keratins<- c('KRT1', 'KRT4', 'KRT5', 'KRT6A', 'KRT6B', 'KRT6C', 'KRT7', 'KRT8', 'KRT10', 'KRT14', 'KRT15',
             'KRT16','KRT17','KRT19','KRT23','KRT78','KRT80')

keratins_ord <- c('KRT4','KRT1','KRT6A','KRT6B','KRT6C','KRT16','KRT78','KRT23','KRT80','KRT17','KRT7',
                  'KRT5','KRT10','KRT19','KRT14','KRT8','KRT15')
s100s_ord <- c('S100A2','S100A10','S100A9','S100A12','S100A11','S100A14','S100A16','S100A7','S100A8','S100A6')
serpins_ord <- c('SERPINB8','SERPINB1','SERPINB3','SERPINB4','SERPINB11','SERPINB2')
sprrs_ord <- c('SPRR3','SPRR1A','SPRR1B','SPRR2F','SPRR2E','SPRR2A','SPRR2B','SPRR2D')

s100s <- c('S100A2','S100A6','S100A7','S100A8','S100A9','S100A10','S100A11','S100A12','S100A14','S100A16')
serpins <- c('SERPINB1','SERPINB2','SERPINB3','SERPINB4','SERPINB8','SERPINB11')
sprrs <- c('SPRR1A','SPRR1B','SPRR2A','SPRR2B','SPRR2D','SPRR2E','SPRR2F','SPRR3')

featurePlot <- function(feature, layers) {
  suppressWarnings({
    lapply(c(layers_rep), function(l) {
      return(VlnPlot(subset(epi_set, layers_rep==l), feature, group.by='location', split.by='condition',
                     cols=cond_cols[c(3,1,2)], pt.size=0) &
               theme(plot.title=element_text(face='bold', size=10), plot.subtitle=element_text(face='italic', size=8, hjust=0.5),
                     axis.text=element_blank(), axis.title=element_text.y(size=8), axis.ticks = element_blank(), 
                     legend.title=element_text(size=8, face='bold'), legend.text=element_text(size=8),
                     axis.title.x=element_blank(), axis.text.x=element_text(size=8)))
    })
  })
}

for (loc in c('Proximal','Distal')) {
  for (gene in leadingEdgeGenes) {
    x <- subset(epi_set, location==loc)
    pdf(file=paste0(output_dir, '/epi_violins',gene,'.',loc,'.pdf'), width=5, height=4)
    print(VlnPlot(x, gene, group.by='layers_rep', split.by='condition',
                  cols=cond_cols[c(3,1,2)], pt.size=0) &
            theme(plot.title=element_text(face='bold', size=10), plot.subtitle=element_text(face='italic', size=8, hjust=0.5),
                  axis.text.y=element_text(size=7), axis.title.y=element_text(size=8), axis.ticks.x = element_blank(), 
                  legend.title=element_text(size=8, face='bold'), legend.text=element_text(size=7),
                  axis.title.x=element_blank(), axis.text.x=element_text(size=7)))
    dev.off()
  }
}

pdf(file=paste0(output_dir, '/epi.dotplot.pdf'), width=6.5, height=9)
DotPlot(epi_set, assay='RNA', features=c(rev(keratins),rev(s100s),rev(sprrs),rev(serpins)), group.by='layers_rep',split.by='condition', 
        cols='Reds') + coord_flip() + scale_color_distiller(palette='OrRd',direction = 0) &
  theme(plot.title=element_text(face='bold', size=10), plot.subtitle=element_text(face='italic', size=8, hjust=0.5),
        axis.text.y=element_text(size=7), axis.ticks.x = element_blank(), 
        legend.title=element_text(size=8, face='bold'), legend.text=element_text(size=7),
        axis.title=element_blank(), axis.text.x=element_text(size=7))
dev.off()


pdf(file=paste0(output_dir, '/epi.dotplot.ord.pdf'), width=6.5, height=9)
DotPlot(epi_set, assay='RNA', features=c(rev(keratins_ord),rev(s100s_ord),rev(sprrs_ord),rev(serpins_ord)), group.by='layers_rep',split.by='condition', 
        cols='Reds') + coord_flip() + scale_color_distiller(palette='OrRd',direction = 0) &
  theme(plot.title=element_text(face='bold', size=10), plot.subtitle=element_text(face='italic', size=8, hjust=0.5),
        axis.text.y=element_text(size=7), axis.ticks.x = element_blank(), 
        legend.title=element_text(size=8, face='bold'), legend.text=element_text(size=7),
        axis.title=element_blank(), axis.text.x=element_text(size=7))
dev.off()

scCustomize::Clustered_DotPlot(epi_set, assay='RNA', features=c(rev(sprrs)), group.by='layers_rep') + coord_flip() + scale_color_distiller(palette='OrRd',direction = 0) &
  theme(plot.title=element_text(face='bold', size=10), plot.subtitle=element_text(face='italic', size=8, hjust=0.5),
        axis.text.y=element_text(size=7), axis.ticks.x = element_blank(), 
        legend.title=element_text(size=8, face='bold'), legend.text=element_text(size=7),
        axis.title=element_blank(), axis.text.x=element_text(size=7))

for (gene in leadingEdgeGenes) {
  pdf(file=paste0(output_dir, '/epi_',gene,'.pdf'), width=5, height=2.5)
  par(mfrow=c(1,5),mar=c(2,0.4,2,0.4))
  for (l in layers_rep) {
    #data <- cbind.data.frame(gene=as.numeric(t(count_norm[gene,])), condition, layer,location)
    print(VlnPlot(subset(epi_set, layers_rep==l), gene, cols=cond_cols[c(3,1,2)], pt.size=0))
    mtext(str_to_title(l), side=3, padj=-1, las=1 ,cex=0.8, col='grey30')
  }
  dev.off()
}








deg_list <- deg_list_ssc_hc_diffbins_p

for (layer in rev(names(deg_list))) {
  print(layer)
  print(' - SSc vs. HC')
  ssc_hc_genes <- get_genes(deg_list_ssc_hc_diffbins_p[[layer]])
  enrich_ssc_hc[[layer]] <- fgseaMultilevel(fgsea_c2_sets, stats=ssc_hc_genes, minSize=10, scoreType='pos')
  print(' - GERD vs. HC')
  #gerd_hc_genes <- get_genes(deg_list_gerd_hc_diffbins[[layer]])
  #enrich_gerd_hc[[layer]] <- fgseaMultilevel(fgsea_c2_sets, stats=gerd_hc_genes, minSize=10, scoreType='pos')
  print(' - SSc vs. GERD')
  #ssc_gerd_genes <- get_genes(deg_list_ssc_gerd_diffbins[[layer]])
  #enrich_ssc_gerd[[layer]] <- fgseaMultilevel(fgsea_c2_sets, stats=ssc_gerd_genes, minSize=10, scoreType='pos')
  print('  - saving...')
  #save(enrich_ssc_hc, enrich_gerd_hc, enrich_ssc_gerd, file=paste0(output_dir, '/enrichment_diffbins.RData'))
}

enrich_dfs <- list()
for (layer in names(deg_list)) {
  enrich_dfs[['ssc_hc']] <- data.frame(row.names=names(fgsea_c2_sets))
  enrich_dfs[['gerd_hc']] <- data.frame(row.names=names(fgsea_c2_sets))
  enrich_dfs[['ssc_gerd']] <- data.frame(row.names=names(fgsea_c2_sets))
  
  enrich_dfs[['ssc_hc']] <- merge(enrich_dfs[['ssc_hc']], )
  enrich_dfs[['ssc_hc']][,layer] <- enrich_ssc_hc[[layer]]
}


# Look to see if genes are same across layers
# SSc
fc_min <- 0
deg_list_ssc_hc_diffbins_sig <- lapply(deg_list_ssc_hc_diffbins, subset, 
                                       abs(avg_log2FC)>fc_min & p_val_adj<0.05)

# get expressed genes
expressed_genes <- c()
for (i in 1:10) {
  expressed_genes <- unique(c(expressed_genes,
                      rownames(deg_list_ssc_hc_diffbins[[i]])
                               [rowSums(deg_list_ssc_hc_diffbins[[i]]
                               [,c('pct.1','pct.2')])>0]))
}

deg_sig_df_ssc_hc <- data.frame(row.names=rownames(deg_list_ssc_hc_diffbins[[1]][]))
for (i in 1:10) {
  print(i)
  
  deg_sig_df_ssc_hc[,paste0('bin_',i)] <- rownames(deg_sig_df_ssc_hc) %in% 
    rownames(deg_list_ssc_hc_diffbins_sig[[i]])
}
#deg_sig_df_ssc_hc <- deg_sig_df_ssc_hc[rowSums(deg_sig_df_ssc_hc)>0,]
x <- ifelse(deg_sig_df_ssc_hc==1,NA,1)
pdf(file=paste0(output_dir, '/Epi_DEGs_x_diffBins_SSc.pdf'), width=12, height=6)
aggr(x, numbers=T, ylabs=c('Proportion of Genes Differentially Expressed','DEGs by Differentiation Bin'), 
     col=c('grey90',cond_cols[1]),combined=F, varheight=T)
dev.off()

# GERD
fc_min <- 0
deg_list_gerd_hc_diffbins_sig <- lapply(deg_list_gerd_hc_diffbins, subset, abs(avg_log2FC)>fc_min & p_val_adj<0.05)

# get expressed genes
expressed_genes <- c()
for (i in 1:10) {
  expressed_genes <- unique(c(expressed_genes,
                              rownames(deg_list_gerd_hc_diffbins[[i]])
                              [rowSums(deg_list_gerd_hc_diffbins[[i]]
                                       [,c('pct.1','pct.2')])>0]))
}

deg_sig_df_gerd_hc <- data.frame(row.names=rownames(deg_list_gerd_hc_diffbins[[1]][]))
for (i in 1:10) {
  print(i)
  
  deg_sig_df_gerd_hc[,paste0('bin_',i)] <- rownames(deg_sig_df_gerd_hc) %in% 
    rownames(deg_list_gerd_hc_diffbins_sig[[i]])
}
#deg_sig_df_gerd_hc <- deg_sig_df_gerd_hc[rowSums(deg_sig_df_gerd_hc)>0,]
x <- ifelse(deg_sig_df_gerd_hc==1,NA,1)
pdf(file=paste0(output_dir, '/Epi_DEGs_x_diffBins_GERD.pdf'), width=12, height=6)
aggr(x, numbers=T, ylabs=c('Proportion of Genes Differentially Expressed','DEGs'), 
     col=c('grey90',cond_cols[2]),combined=F, varheight=T)
dev.off()






#### DEclust analysis ####

# Do pseudobulk

# Make input files
