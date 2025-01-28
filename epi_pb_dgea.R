
.libPaths('~/local/R_libs/')    # Set path using .libPaths function

library(Seurat)
library(RColorBrewer)
library(edgeR)
library(pals)
library(stringr)
library(beeswarm)
library(matrixStats)
library(scales)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)

####  SET PARAMETERS  ##############################################################################

set.seed(0123456789)

cond_cols <- c('#8156B3','#e26b53','#A8A39d')
layer_cols_3 <- c('#4ca5b1','#e9b85d','#c72f4c')
layer_cols_5 <- c('#4B86B8', '#7DC0A6','#E9F5A3','#E1BA6C','#B73D4F')

cond='condition'; loc='location'; 

comps <- list(c('SSc','HC'), c('GERD','HC'), c('SSc', 'GERD'))

layers <- c('basal','suprabasal','superficial')
layers_rep <- c('Basal','ProliferatingBasal','ProliferatingSuprabasal','Suprabasal','Superficial')


s_h_1 <- 'SSc_HC_basal'; s_h_2 <- 'SSc_HC_suprabasal'; s_h_3 <- 'SSc_HC_superficial';
g_h_1 <- 'GERD_HC_basal'; g_h_2 <- 'GERD_HC_suprabasal'; g_h_3 <- 'GERD_HC_superficial';
s_g_1 <- 'SSc_GERD_basal'; s_g_2 <- 'SSc_GERD_suprabasal'; s_g_3 <- 'SSc_GERD_superficial';

output_dir <- '/projects/p31763/users/mdapas/MP_SSc_scRNAseq/dapas_analysis/epi_analysis/'
#output_dir <- '/Users/mdn578/Documents/Research/Winter/Tetrault_SSc/Analysis/pseudobulk_analysis/'


####  FUNCTIONS  ##############################################################################


#Load Options
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}


rntransform <- function(x) {
  var <- ztransform(x)
  out <- rank(var) - 0.5
  out[is.na(var)] <- NA
  out <- out/(max(out,na.rm=T)+.5)
  out <- qnorm(out)
  return(out)
}

ztransform <- function(x) {
  x <- x[!is.na(x)]
  return((x-mean(x))/sd(x))
}

#-------------------------------------------
# Median Ratio Normalization function (MRN)
#-------------------------------------------

mrnFactors=function(rawCounts,conditions) {
  rawCounts <- as.matrix(rawCounts)
  totalCounts <- colSums(rawCounts)
  normFactors <- totalCounts
  medianRatios <- rep(1,length(conditions))
  names(medianRatios) <- names(normFactors)
  if (sum(conditions==1)>1)
    meanA <- apply(rawCounts[,conditions==1]%*%diag(1/totalCounts[conditions==1]),1,mean)
  else
    meanA <- rawCounts[,conditions==1]/totalCounts[conditions==1]
  for (i in 2:max(conditions)) {
    if (sum(conditions==i)>1)
      meanB <- apply(rawCounts[,conditions==i]%*%diag(1/totalCounts[conditions==i]),1,mean)
    else
      meanB <- rawCounts[,conditions==i]/totalCounts[conditions==i]
    meanANot0 <- meanA[meanA>0&meanB>0]
    meanBNot0 <- meanB[meanA>0&meanB>0]
    ratios <- meanBNot0/meanANot0
    medianRatios[conditions==i] <- median(ratios)
    normFactors[conditions==i] <- medianRatios[conditions==i]*totalCounts[conditions==i]
  }
  medianRatios <- medianRatios/exp(mean(log(medianRatios)))
  normFactors <- normFactors/exp(mean(log(normFactors)))
  return(list(medianRatios=medianRatios,normFactors=normFactors))
}

#---------
# The End
#---------


####  MAIN  ##############################################################################

epi_set <- loadRData(paste0(output_dir, '/MP_eso_epi_reintegratedObj.RData'))

epi_set@meta.data$layers_rep <- str_to_title(gsub('replicating', 'proliferating', epi_set@meta.data$layers_rep))
epi_set@meta.data$layers_rep <- gsub(' ','',epi_set@meta.data$layers_rep)
epi_set@meta.data$layers_rep <- factor(epi_set@meta.data$layers_rep, 
                                       levels=names(table(epi_set@meta.data$layers_rep))[c(1,2,3,5,4)])

#epi_super <- subset(epi_set, layers=='superficial')

# Aggregate expression
epi_aggro <- AggregateExpression(epi_set, group.by = c('condition','location', 'layers_rep'), 
                                 return.seurat=TRUE, assays='RNA', slot='counts')
#load(paste0(output_dir,'epi_aggro.RData'))
save(epi_aggro, file=paste0(output_dir,'epi_aggro_condLocLayersRep.RData'))


pb_obj <- DGEList(epi_aggro@assays$RNA@counts)  # create DGEList
#pb_obj <- calcNormFactors(pb_obj)

# add metadata to psuedobulk object
pb_obj$samples$sample <- rownames(pb_obj$samples)
pb_obj$samples$condition <- ifelse(grepl('^I', pb_obj$samples$sample), 'GERD',
                                   ifelse(grepl('^P', pb_obj$samples$sample), 'SSc', 'HC'))
pb_obj$samples$layer <- word(rownames(pb_obj$samples),2,sep='\\_')
pb_obj$samples$location <- ifelse(grepl('D_', pb_obj$samples$sample), 'Distal','Proximal')
pb_obj$samples$group <- paste(pb_obj$samples$condition, pb_obj$samples$layer, sep='_')

samples <- factor(sub("\\-.*$","", pb_obj$samples$sample))
sample_cols <- sample(alphabet(26),20)

layer <- factor(pb_obj$samples$layer, levels=names(table(pb_obj$samples$layer))[c(1,2,3,5,4)])
condition <- factor(pb_obj$samples$condition, levels=c('SSc','GERD','HC'))
location <- factor(pb_obj$samples$location, levels=c('Distal','Proximal'))
group <- factor(pb_obj$samples$group, levels=unique(pb_obj$samples$group))

# filter non-expressed genes by subtype (logCPM >0 in >50%)
count_cpm <- cpm(pb_obj)
count_lcpm <- cpm(pb_obj, log=T)
min_cpm <- 1
ct.genes <- lapply(setNames(levels(layer), levels(layer)),
                   function(ct){
                     ct.df <- as.data.frame(count_cpm[,layer==ct])
                     ct.df.l <- as.data.frame(count_lcpm[,layer==ct])
                     ct.ssc.df <- as.data.frame(count_cpm[,layer==ct & condition=='SSc'])
                     ct.gerd.df <- as.data.frame(count_cpm[,layer==ct & condition=='GERD'])
                     ct.hc.df <- as.data.frame(count_cpm[,layer==ct & condition=='HC'])
                     ct.ssc.genes <- rownames(ct.ssc.df[rowSums(ct.ssc.df>min_cpm)>ncol(ct.ssc.df)*0.75,])
                     ct.gerd.genes <- rownames(ct.gerd.df[rowSums(ct.gerd.df>min_cpm)>ncol(ct.gerd.df)*0.75,])
                     ct.hc.genes <- rownames(ct.hc.df[rowSums(ct.hc.df>min_cpm)>=ncol(ct.hc.df)*0.75,])
                     ct.outlier.genes <- rownames(ct.df.l[rowSums(ct.df.l>(rowMeans(ct.df.l)+3*rowSds(as.matrix(ct.df.l))))==1,])
                     ct.outlier.genes2 <- rownames(ct.df.l[rowSums(ct.df.l>(rowMeans(ct.df.l)+3*rowSds(as.matrix(ct.df.l))))==2,])
                     for (gene in ct.outlier.genes2) {
                       top_two <- sub("\\-.*$","", names(ct.df.l[gene,order(as.numeric(ct.df.l[gene,]),decreasing=T)][1:2]))
                       if (top_two[1]==top_two[2]) {
                         ct.outlier.genes <- c(ct.outlier.genes, gene)}}
                     ct.genes <- unique(c(ct.ssc.genes,ct.gerd.genes,ct.hc.genes))
                     ct.genes <- ct.genes[!ct.genes %in% ct.outlier.genes]
                     return(ct.genes)
                   })

pb_obj <- pb_obj[unique(unlist(ct.genes)), , keep=FALSE] 


# Panel of MDS plots by different variables
x_range=c(-3,6.2); y_range=c(-2.1,3.5)
pdf(file=paste0(output_dir, '/PseudoDEG_epi_layerRep_MDS.pdf'), width=11, height=11)
par(mfrow=c(2,2))
plotMDS(pb_obj, pch=21, bg=sample_cols[samples], main="MDS by Sample", cex=1.3, xlim=x_range)
legend("bottomright", legend=levels(samples),
       pch=21, pt.bg=sample_cols, cex=0.65)

plotMDS(pb_obj, pch=21, bg=layer_cols_5[layer], main="MDS by Epithelial Compartment", cex=1.3, xlim=x_range)
legend("bottomright", legend=levels(layer),
       pch=21, pt.bg=layer_cols_5, cex=0.8)

plotMDS(pb_obj, pch=21, bg=cond_cols[condition], main="MDS by Condition", cex=1.3, xlim=x_range)
legend("bottomright", legend=levels(condition),
       pch=21, pt.bg=cond_cols, cex=0.8)

plotMDS(pb_obj, pch=21, bg=c(6,13)[location], main="MDS by Biopsy Location", cex=1.3, xlim=x_range)
legend("bottomright", legend=levels(location),
       pch=21, pt.bg=c(6,13), cex=0.8)
dev.off()

# create design matrix
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))

# fit EdgeR model
pb_obj <- estimateDisp(pb_obj, design, robust=TRUE)
pb_obj$common.dispersion
plotBCV(pb_obj)

fit <- glmQLFit(pb_obj, design, robust=T)
plotQLDisp(fit)

fit_lrt <- glmFit(pb_obj, design, robust=TRUE)

res <- glmQLFTest(fit, coef=ncol(design)) #?
summary(decideTests(res)) #?

# Make contrasts matrix (for DEG comparisons)
contrast_comps <- c(); contrast_names <- c()
for (c in comps) {
  for (l in layers) {
    contrast_names <- c(contrast_names, paste(c[1],c[2],l,sep='_'))
    contrast_comps <- c(contrast_comps, paste0(paste(c[1],l,sep='_'), '-', paste(c[2],l, sep='_')))
  }
}

contrasts <- makeContrasts(contrasts=contrast_comps, levels=design)
colnames(contrasts) <- contrast_names

# Loop through contrasts and compute DEGs, save to list
qlf.res <- lapply(setNames(1:ncol(contrasts), names(contrasts[1,])),
                  function(i){
                    pb_obj.ct <- pb_obj[ct.genes[[strsplit(colnames(contrasts)[i],'_')[[1]][3]]], , keep=FALSE] 
                    fit <- glmQLFit(pb_obj.ct, design, robust=T)
                    qlf <- glmQLFTest(fit, contrast=contrasts[,i])
                    qlf$comparison <- colnames(contrasts)[i]
                    qlf$table$p_adj <- p.adjust(qlf$table$PValue, method='fdr')
                    return(qlf$table)
                  })

lrt.res <- lapply(setNames(1:ncol(contrasts), names(contrasts[1,])),
                  function(i){
                    pb_obj.ct <- pb_obj[ct.genes[[strsplit(colnames(contrasts)[i],'_')[[1]][3]]], , keep=FALSE] 
                    fit_lrt <- glmFit(pb_obj.ct, design, robust=TRUE)
                    lrt <- glmLRT(fit_lrt, contrast=contrasts[,i])
                    lrt$comparison <- colnames(contrasts)[i]
                    lrt$table$p_adj <- p.adjust(lrt$table$PValue, method='fdr')
                    return(lrt$table)
                  })


# Alternative approach, using limma voom
limma.res <- lapply(setNames(1:ncol(contrasts), names(contrasts[1,])),
                    function(i) {
                      pb_obj.ct <- pb_obj[ct.genes[[strsplit(colnames(contrasts)[i],'_')[[1]][3]]], , keep=FALSE] 
                      voom_weights <- voom(pb_obj.ct, design, plot=T)
                      limma_fit <- lmFit(voom_weights, design)
                      limma_results <- contrasts.fit(limma_fit, contrasts)
                      limma_results <- eBayes(limma_results)
                      r <- cbind(as.data.frame(limma_results$coefficients)[colnames(contrasts)[i]],
                                 as.data.frame(limma_results$p.value)[colnames(contrasts)[i]],
                                 p.adjust(as.data.frame(limma_results$p.value)[colnames(contrasts)[i]][[1]], method='fdr'))
                      names(r) <- c('logFC','p','p_adj')
                      return(r)})


# compute manually using Wilcoxon test and evaluate correlation

## Perform TMM normalization and convert to CPM (Counts Per Million)
#d <- calcNormFactors(d, method = "TMM")  # already did this
count_norm <- cpm(pb_obj, log=T)
count_norm <- as.data.frame(count_norm)

# Run the Wilcoxon rank-sum test for each gene
wilcox.res <- list()
for (c in comps) {
  print(c)
  comp_results <- list()
  for (l in layers) {
    print(l)
    pb_obj.ct <- pb_obj[ct.genes[[l]], , keep=FALSE] 
    count_norm <- cpm(pb_obj.ct, log=T)
    count_norm <- as.data.frame(count_norm)
    gene.results <- lapply(1:nrow(pb_obj.ct$counts), 
                           function(i){
                             gene.dat <- cbind.data.frame(gene=as.numeric(t(count_norm[i,])), condition, layer)
                             gene.test <- wilcox.test(jitter(gene, amount=0.0000001)~condition, 
                                                      gene.dat[(gene.dat$condition %in% c) & (gene.dat$layer==l),])
                             gene.stats <- data.frame(gene=rownames(count_norm)[i],
                                                      fc=mean(gene.dat[gene.dat$condition==c[1] & gene.dat$layer==l,'gene']) -
                                                        mean(gene.dat[gene.dat$condition==c[2] & gene.dat$layer==l,'gene']),
                                                      p=gene.test$p.value, stringsAsFactors = F)
                             return(gene.stats)
                           })
    layer_results <- as.data.frame(do.call(rbind, gene.results))
    comp_results[[paste(c[1],c[2],l,sep='_')]] <- layer_results
  }
  wilcox.res <- c(wilcox.res, comp_results)
}

cor_fc_list <- lapply(setNames(names(contrasts[1,]), names(contrasts[1,])),
                      function(c){
                        x <- cbind(qlf.res[[c]]$logFC, lrt.res[[c]]$logFC, limma.res[[c]]$logFC, wilcox.res[[c]]$fc)
                        cor_fc <- cor(x, method='spearman')
                        colnames(cor_fc) <- c('edgeR_QLF', 'edgeR_LRT','limma_voom','wilcox')
                        rownames(cor_fc) <- c('edgeR_QLF', 'edgeR_LRT','limma_voom','wilcox')
                        return(cor_fc)
                      })
cor_p_list <- lapply(setNames(names(contrasts[1,]), names(contrasts[1,])),
                     function(c){
                       x <- cbind(qlf.res[[c]]$PValue, lrt.res[[c]]$PValue, limma.res[[c]]$p, wilcox.res[[c]]$p)
                       cor_fc <- cor(x, method='spearman')
                       colnames(cor_fc) <- c('edgeR_QLF', 'edgeR_LRT','limma_voom','wilcox')
                       rownames(cor_fc) <- c('edgeR_QLF', 'edgeR_LRT','limma_voom','wilcox')
                       return(cor_fc)
                     })

save(qlf.res, lrt.res, limma.res, wilcox.res,
     cor_fc_list, cor_p_list, file=paste0(output_dir, '/pbDEG_dat.RData'))

load(paste0(output_dir, '/pbDEG_dat.RData'))
#### Compare DEG profiles for SSc vs. HC and GERD vs. HC for epithelial layers ####
#  (there is probably a more efficient way to do this)  

comp_cols <- lapply(setNames(layers, layers), function(l) {
  cols <- factor(ifelse(qlf.res[[paste('SSc','HC',l,sep='_')]]$p_adj<0.05,
                        ifelse(qlf.res[[paste('GERD','HC',l,sep='_')]]$p_adj<0.05,'cornflowerblue',cond_cols[1]),
                        ifelse(qlf.res[[paste('GERD','HC',l,sep='_')]]$p_adj<0.05,cond_cols[2], 'grey90')),
                 levels=c('grey90',cond_cols[2], cond_cols[1],'cornflowerblue'))
  return(cols)
})

pdf(file=paste0(output_dir, '/PseudoDEG_epi_scatterFC_SScGERD_HCsig.pdf'), width=10, height=4)
xrange<-c(-6,6)
yrange<-c(-6,6)
par(mfrow=c(1,3))
for(l in layers) {
  plot(qlf.res[[paste('GERD','HC',l,sep='_')]][order(comp_cols[[l]]),'logFC'], 
       qlf.res[[paste('SSc','HC',l,sep='_')]][order(comp_cols[[l]]),'logFC'],
       pch=16,col=as.character(sort(comp_cols[[l]])), cex=0.8,
       ylab='SSc vs. HC (logFC)', xlab='GERD vs. HC (logFC)', main=str_to_title(l), xlim=xrange,ylim=yrange)
  text(-5.3,6,paste('r =',round(cor(qlf.res[[paste('GERD','HC',l,sep='_')]]$logFC, qlf.res[[paste('SSc','HC',l,sep='_')]]$logFC),2)))
  abline(0,1,lty=2, col='grey80')
}
dev.off()

comp_cols <- lapply(setNames(layers, layers), function(l) {
  cols <- factor(ifelse(qlf.res[[paste('SSc','HC',l,sep='_')]]$p_adj<0.05 & 
                          qlf.res[[paste('GERD','HC',l,sep='_')]]$p_adj>0.05 & 
                          qlf.res[[paste('SSc','GERD',l,sep='_')]]$p_adj<0.05, cond_cols[1],
                        ifelse(qlf.res[[paste('SSc','HC',l,sep='_')]]$p_adj>0.05 & 
                                 qlf.res[[paste('GERD','HC',l,sep='_')]]$p_adj<0.05 & 
                                 qlf.res[[paste('SSc','GERD',l,sep='_')]]$p_adj<0.05, cond_cols[2],
                               'grey90')),
                 levels=c('grey90',cond_cols[2], cond_cols[1]))
  return(cols)
})

pdf(file=paste0(output_dir, '/PseudoDEG_epi_scatterFC_SSc_GERDsig.pdf'), width=10, height=4)
xrange<-c(-6,6)
yrange<-c(-6,6)
par(mfrow=c(1,3))
for(l in layers) {
  plot(qlf.res[[paste('GERD','HC',l,sep='_')]][order(comp_cols[[l]]),'logFC'], 
       qlf.res[[paste('SSc','HC',l,sep='_')]][order(comp_cols[[l]]),'logFC'],
       pch=16,col=as.character(sort(comp_cols[[l]])), cex=0.8,
       ylab='SSc vs. HC (logFC)', xlab='GERD vs. HC (logFC)', main=str_to_title(l), xlim=xrange,ylim=yrange)
  text(-5.3,6,paste('r =',round(cor(qlf.res[[paste('GERD','HC',l,sep='_')]]$logFC, qlf.res[[paste('SSc','HC',l,sep='_')]]$logFC, method='spearman'),2)))
  abline(0,1,lty=2, col='grey80')
}
dev.off()


# Plot individual genes
count_norm <- cpm(pb_obj, log=T)
count_norm <- as.data.frame(count_norm)

genes=c('NEFL','NEFM','IFI44L','KRT1','KRTDAP','MUC22')
genes=c('CCL19','IL33','CDKN2A','CCL2','PEG3','LCE3D')

for (gene in genes) {
  pdf(file=paste0(output_dir, '/PseudoDEG_epiLayers_',gene,'.pdf'), width=8, height=3)
  par(mfrow=c(1,length(layers)))
  for (l in layers) {
    data <- cbind.data.frame(gene=as.numeric(t(count_norm[gene,])), condition, layer)
    data$condition <- factor(data$condition, c('HC','GERD','SSc'))
    boxplot(gene~condition, data=data[data$layer==l,], ylab='log(CPM)', col=alpha(rev(cond_cols),0.75),
            xlab='', pch=16, cex.main=1.5, ylim=c(min(count_norm[gene,]),max(count_norm[gene,])), outline=F)
    points(gene~jitter(as.numeric(condition), amount=0.1), data=data[data$layer==l,],pch=16, cex=0.8)
    title(main=bquote(italic(.(gene))), line=3, cex.main=1.5, font.main=1.8)
    mtext(str_to_title(l), side=3, padj=-1, las=1 ,cex=0.8, col='grey30')
  }
  dev.off()
}

#### Pathway enrichment Analysis ####
annot_go <- read.table(paste0(output_dir,'../pathway_gene_maps/GO_Biological_Process_2023.map'), 
                       header=F, sep='\t', quote="")
annot_kegg <-  read.table(paste0(output_dir,'../pathway_gene_maps/KEGG_2021_Human.map'), 
                          header=F, sep='\t', quote="")
annot_tf <-  read.table(paste0(output_dir,'../pathway_gene_maps/ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.map'), 
                        header=F, sep='\t', quote="")
annot_gedipnet <-  read.table(paste0(output_dir,'../pathway_gene_maps/GeDiPNet_2023.map'), 
                              header=F, sep='\t', quote="")


go_enrich.res <- lapply(setNames(layers, layers), function(ct) {
  background_genes <- ct.genes[[ct]]
  background_genes <- background_genes[!grepl('^MRP',background_genes)]
  background_genes <- background_genes[!grepl('^RP',background_genes)]
  
  test_genes <- rownames(qlf.res[[paste0('SSc_HC_',ct)]][(qlf.res[[paste0('SSc_HC_',ct)]]$p_adj<0.05) &
                                                           (qlf.res[[paste0('SSc_GERD_',ct)]]$p_adj<0.05) &
                                                           (qlf.res[[paste0('GERD_HC_',ct)]]$p_adj>0.05) &
                                                           (qlf.res[[paste0('SSc_HC_',ct)]]$logFC>0),])
  test_genes <- test_genes[!grepl('^MRP',test_genes)]
  test_genes <- test_genes[!grepl('^RP',test_genes)]
  
  ego <- enricher(gene          = test_genes,
                  universe      = background_genes,
                  pAdjustMethod = "fdr",
                  pvalueCutoff  = 1,
                  TERM2GENE     = annot_go)
  
  #View(as.data.frame(ego@result))
  ego_df <- as.data.frame(ego@result)
  ego_df <- ego_df[ego_df$Count>2,]
  ego_df$p.adjust <- p.adjust(ego_df$pvalue)
  ego_df <- ego_df[,c(1,3:ncol(ego_df))]
  rownames(ego_df) <- NULL
  View(ego_df)
  return(ego_df[ego_df$p.adjust<0.05,])
})

kegg_enrich.res <- lapply(setNames(layers, layers), function(ct) {
  background_genes <- ct.genes[[ct]]
  background_genes <- background_genes[!grepl('^MRP',background_genes)]
  background_genes <- background_genes[!grepl('^RP',background_genes)]
  
  test_genes <- rownames(qlf.res[[paste0('SSc_HC_',ct)]][(qlf.res[[paste0('SSc_HC_',ct)]]$p_adj<0.05) &
                                                           (qlf.res[[paste0('SSc_GERD_',ct)]]$p_adj<0.05) &
                                                           (qlf.res[[paste0('GERD_HC_',ct)]]$p_adj>0.05) &
                                                           (qlf.res[[paste0('SSc_HC_',ct)]]$logFC>0),])
  test_genes <- test_genes[!grepl('^MRP',test_genes)]
  test_genes <- test_genes[!grepl('^RP',test_genes)]
  
  ego <- enricher(gene          = test_genes,
                  universe      = background_genes,
                  pAdjustMethod = "fdr",
                  pvalueCutoff  = 1,
                  TERM2GENE     = annot_kegg)
  
  #View(as.data.frame(ego@result))
  ego_df <- as.data.frame(ego@result)
  ego_df <- ego_df[ego_df$Count>2,]
  ego_df$p.adjust <- p.adjust(ego_df$pvalue)
  ego_df <- ego_df[,c(1,3:ncol(ego_df))]
  rownames(ego_df) <- NULL
  #View(ego_df)
  return(ego_df[ego_df$p.adjust<0.05,])
})


tf_enrich.res <- lapply(setNames(layers, layers), function(ct) {
  background_genes <- ct.genes[[ct]]
  background_genes <- background_genes[!grepl('^MRP',background_genes)]
  background_genes <- background_genes[!grepl('^RP',background_genes)]
  
  test_genes <- rownames(qlf.res[[paste0('SSc_HC_',ct,'_Proximal')]][(qlf.res[[paste0('SSc_HC_',ct,'_Proximal')]]$p_adj<0.05) &
                                                                       (qlf.res[[paste0('GERD_HC_',ct,'_Proximal')]]$p_adj>0.05) &
                                                                       (qlf.res[[paste0('SSc_HC_',ct,'_Proximal')]]$logFC>-200),])
  test_genes <- test_genes[!grepl('^MRP',test_genes)]
  test_genes <- test_genes[!grepl('^RP',test_genes)]
  
  ego <- enricher(gene          = test_genes,
                  universe      = background_genes,
                  pAdjustMethod = "fdr",
                  pvalueCutoff  = 1,
                  TERM2GENE     = annot_tf)
  
  #View(as.data.frame(ego@result))
  ego_df <- as.data.frame(ego@result)
  ego_df <- ego_df[ego_df$Count>2,]
  ego_df$p.adjust <- p.adjust(ego_df$pvalue)
  ego_df <- ego_df[,c(1,3:ncol(ego_df))]
  rownames(ego_df) <- NULL
  #View(ego_df)
  return(ego_df[ego_df$p.adjust<0.05,])
})

for (ct in layers) {
  n_sig_pathways <- nrow(as.data.frame(go_enrich.res[[ct]]))
  if (n_sig_pathways>0){
    pdf(file=paste0(output_dir, '/PseudoDEG_enrich_GO_',ct,'_dotplot.pdf'), width=6, height=8)
    print(dotplot(go_enrich.res[[ct]], showCategory=min(n_sig_pathways,20)) + ggtitle(ct) + 
            theme(axis.text.y=element_text(size=8), legend.title=element_text(size=8, face='bold'),
                  plot.title=element_text(hjust=0.5, size=10,face='bold')) + 
            scale_color_gradient(name='q-value',low='rosybrown2',high='red3'))
    dev.off()
  }
  write.table(as.data.frame(go_enrich.res[[ct]]), 
              paste0(output_dir, '/PseudoDEG_enrich_GO_',ct,'.tsv'), quote=F, row.names=F,col.names=T, sep="\t")
}





# Loop through contrasts for all genes included in each comparison (for plotting purposes)
qlf.res.all <- lapply(setNames(1:ncol(contrasts), names(contrasts[1,])),
                      function(i){
                        pb_obj.ct <- pb_obj[unique(unlist(ct.genes)), , keep=FALSE] 
                        fit <- glmQLFit(pb_obj.ct, design, robust=T)
                        qlf <- glmQLFTest(fit, contrast=contrasts[,i])
                        qlf$comparison <- colnames(contrasts)[i]
                        qlf$table$p_adj <- p.adjust(qlf$table$PValue, method='fdr')
                        return(qlf$table)
                      })

# Make heatmap of top genes
sig_genes <- lapply(1:6, function(i){
  return(rownames(qlf.res[[i]][qlf.res[[i]]$p_adj<0.05,]))
})
sig_genes <- unique(unlist(sig_genes))


lfc <- lapply(setNames(1:ncol(contrasts), names(contrasts[1,])), function(i){
  return(qlf.res.all[[i]][sig_genes,'logFC'])})

lfc <- data.frame(lfc)[,c(1,4,2,5,3,6)]  # save fold change values for heatmap
rownames(lfc) <- sig_genes

annot <- data.frame(Condition=rep(c('SSc','GERD'),3),
                    Layer=c(rep('basal',2), rep('suprabasal',2),rep('superficial',2)))

rownames(annot) <- paste(annot$Condition,'HC', annot$Layer, sep='_')
ann_colors <- list(layer=2:4)
names(ann_colors$layer) <- paste0("layer ", levels(layer))

jpeg(paste0(output_dir, '/pseudoDEG_heatmap_colClust_SScGERD_HC.jpg'), width=7, height=9,units='in',res=150)
pheatmap::pheatmap(lfc[sig_genes,], breaks=seq(-5,5,length.out=101),
                   color=colorRampPalette(c("blue","white","red"))(100), scale="none",
                   cluster_cols=T, border_color="NA", fontsize_row=5,
                   treeheight_row=70, treeheight_col=70, cutree_cols=3,
                   clustering_method="ward.D2", show_colnames=FALSE, show_rownames = F,
                   annotation_col=annot, legend_labels=c('LogFC vs. HC'),
                   annotation_colors=list(Layer=c(basal=layer_cols_3[1], suprabasal=layer_cols_3[2],
                                                  superficial=layer_cols_3[3]),
                                          Condition=c(SSc=cond_cols[1], GERD=cond_cols[2]) ))
dev.off()




# Make heatmap of top genes
sig_genes <- lapply(7:9, function(i){
  return(rownames(qlf.res[[i]][qlf.res[[i]]$p_adj<0.05,]))
})
sig_genes <- unique(unlist(sig_genes))


lfc <- lapply(setNames(1:ncol(contrasts), names(contrasts[1,])), function(i){
  return(qlf.res.all[[i]][sig_genes,'logFC'])})

lfc <- data.frame(lfc)[,c(7,8,9)]  # save fold change values for heatmap
rownames(lfc) <- sig_genes

annot <- data.frame(Layer=c('basal','suprabasal','superficial'))
rownames(annot) <- paste('SSc_GERD', annot$Layer, sep='_')
ann_colors <- list(layer=2:4)
names(ann_colors$layer) <- paste0("layer ", levels(layer))

jpeg(paste0(output_dir, '/pseudoDEG_heatmap_SSc_GERD.jpg'), width=7, height=9,units='in',res=150)
pheatmap::pheatmap(lfc[sig_genes,], breaks=seq(-5,5,length.out=101),
                   color=colorRampPalette(c("blue","white","red"))(100), scale="none",
                   cluster_cols=F, border_color="NA", fontsize_row=5,
                   treeheight_row=70, treeheight_col=70, cutree_cols=3,
                   clustering_method="ward.D2", show_colnames=FALSE, show_rownames = F,
                   annotation_col=annot, legend_labels=c('SSc vs.'),
                   annotation_colors=list(Layer=c(basal=layer_cols_3[1], suprabasal=layer_cols_3[2],
                                                  superficial=layer_cols_3[3])) )
dev.off()



# Check direction of fold change convention
count_norm <- cpm(pb_obj, log=T)
count_norm <- as.data.frame(count_norm)
head(qlf.res[[3]][order(qlf.res[[3]]$logFC),],10)
gene='IL36A'
genes <- rownames(head(qlf.res[[3]][order(qlf.res[[3]]$logFC),]))
par(mfrow=c(1,length(genes)))
for (gene in genes) {
  data <- cbind.data.frame(gene=as.numeric(t(count_norm[gene,])), condition, layer)
  boxplot(gene~condition, data=data[(data$condition %in% c('HC','SSc')) & 
                                      (data$layer=='superficial'),],
          ylab=gene)
  beeswarm(gene~condition, data=data[(data$condition %in% c('HC','SSc')) & 
                                       (data$layer=='superficial'),], ylab=gene, pch=16, add=T)
}




# Aggregate expression by diff bins
epi_aggro_diffbins <- AggregateExpression(epi_set, group.by = c('orig.ident', 'diff.bin'), 
                                          return.seurat=TRUE, assays='RNA', slot='counts')
save(epi_aggro_diffbins, file=paste0(output_dir,'epi_aggro_diffbins.RData'))
load(file=paste0(output_dir,'epi_aggro_diffbins.RData'))

pb_obj <- DGEList(epi_aggro_diffbins@assays$RNA@counts)  # create DGEList

# add metadata to psuedobulk object
pb_obj$samples$sample <- rownames(pb_obj$samples)
pb_obj$samples$condition <- ifelse(grepl('^I', pb_obj$samples$sample), 'GERD',
                                   ifelse(grepl('^P', pb_obj$samples$sample), 'SSc', 'HC'))
pb_obj$samples$diff.bin <- word(rownames(pb_obj$samples),2,sep='\\_')
pb_obj$samples$location <- ifelse(grepl('D_', pb_obj$samples$sample), 'Distal','Proximal')
pb_obj$samples$group <- paste(pb_obj$samples$condition, pb_obj$samples$diff.bin, sep='_')

samples <- factor(sub("\\-.*$","", pb_obj$samples$sample))
sample_cols <- sample(alphabet(26),20)
diffbin_cols <- rev(brewer.pal(10, "Spectral"))

diff.bin <- factor(pb_obj$samples$diff.bin, levels=as.character(1:10))
condition <- factor(pb_obj$samples$condition, levels=c('SSc','GERD','HC'))
location <- factor(pb_obj$samples$location, levels=c('Distal','Proximal'))
group <- factor(pb_obj$samples$group, levels=unique(pb_obj$samples$group))


# remove non-expressed genes
count_norm <- cpm(pb_obj, log=T)
min_cpm <- log(0.1)
b.genes <- lapply(setNames(levels(diff.bin), levels(diff.bin)),
                  function(b){
                    b.df <- as.data.frame(count_norm[,diff.bin==b])
                    b.ssc.df <- as.data.frame(count_norm[,diff.bin==b & condition=='SSc'])
                    b.gerd.df <- as.data.frame(count_norm[,diff.bin==b & condition=='GERD'])
                    b.hc.df <- as.data.frame(count_norm[,diff.bin==b & condition=='HC'])
                    b.ssc.genes <- rownames(b.ssc.df[rowSums(b.ssc.df>min_cpm)>ncol(b.ssc.df)*0.75,])
                    b.gerd.genes <- rownames(b.gerd.df[rowSums(b.gerd.df>min_cpm)>ncol(b.gerd.df)*0.75,])
                    b.hc.genes <- rownames(b.hc.df[rowSums(b.hc.df>min_cpm)>ncol(b.hc.df)*0.75,])
                    b.outlier.genes <- rownames(b.df[rowSums(b.df>(rowMeans(b.df)+3*rowSds(as.matrix(b.df))))==1,])
                    b.outlier.genes2 <- rownames(b.df[rowSums(b.df>(rowMeans(b.df)+3*rowSds(as.matrix(b.df))))==2,])
                    for (gene in b.outlier.genes2) {
                      top_two <- sub("\\-.*$","", names(b.df[gene,order(as.numeric(b.df[gene,]),decreasing=T)][1:2]))
                      if (top_two[1]==top_two[2]) {
                        b.outlier.genes <- c(b.outlier.genes, gene)}}
                    print(c(length(b.ssc.genes),length(b.gerd.genes),length(b.hc.genes),length(b.outlier.genes),length(b.outlier.genes2) ))
                    b.genes <- unique(c(b.ssc.genes,b.gerd.genes,b.hc.genes))
                    b.genes <- b.genes[!b.genes %in% b.outlier.genes]
                    return(b.genes)
                  })

pb_obj <- pb_obj[unique(unlist(b.genes)), , keep=FALSE] 

##### Panel of MDS plots by different variables ####
pdf(file=paste0(output_dir, '/PseudoDEG_epiDiffBins_MDS.pdf'), width=11, height=11)
par(mfrow=c(2,2))
plotMDS(pb_obj, pch=21, bg=sample_cols[samples], main="MDS by Sample", cex=1.3)
legend("topright", legend=levels(samples),
       pch=21, pt.bg=sample_cols, cex=0.65)

plotMDS(pb_obj, pch=21, bg=diffbin_cols[diff.bin], main="MDS by Epithelial Layer", cex=1.3)
legend("topright", legend=levels(diff.bin),
       pch=21, pt.bg=diffbin_cols, cex=0.8)

plotMDS(pb_obj, pch=21, bg=cond_cols[condition], main="MDS by Condition", cex=1.3)
legend("topright", legend=levels(condition),
       pch=21, pt.bg=cond_cols, cex=0.8)

plotMDS(pb_obj, pch=21, bg=c(6,13)[location], main="MDS by Biopsy Location", cex=1.3)
legend("topright", legend=levels(location),
       pch=21, pt.bg=c(6,13), cex=0.8)
dev.off()
#####

# create design matrix
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))

# fit EdgeR model
pb_obj <- estimateDisp(pb_obj, design, robust=TRUE)
pb_obj$common.dispersion
plotBCV(pb_obj)

fit <- glmQLFit(pb_obj, design, robust=T)
plotQLDisp(fit)

fit_lrt <- glmFit(pb_obj, design, robust=TRUE)

res <- glmQLFTest(fit, coef=ncol(design)) #?
summary(decideTests(res)) #?

# Make contrasts matrix (for DEG comparisons)
contrast_comps <- c(); contrast_names <- c()
for (c in comps) {
  for (l in 1:10) {
    contrast_names <- c(contrast_names, paste(c[1],c[2],l,sep='_'))
    contrast_comps <- c(contrast_comps, paste0(paste(c[1],l,sep='_'), '-', paste(c[2],l, sep='_')))
  }
}

contrasts <- makeContrasts(contrasts=contrast_comps, levels=design)
colnames(contrasts) <- contrast_names

# Loop through contrasts and compute DEGs, save to list
qlf.res <- lapply(setNames(1:ncol(contrasts), names(contrasts[1,])),
                  function(i){
                    print(c(i,names(contrasts[1,])[i]))
                    pb_obj.b <- pb_obj[b.genes[[strsplit(colnames(contrasts)[i],'_')[[1]][3]]], , keep=FALSE] 
                    pb_obj.b <- calcNormFactors(pb_obj.b, method='RLE')
                    fit <- glmQLFit(pb_obj.b, design, robust=T)
                    qlf <- glmQLFTest(fit, contrast=contrasts[,i])
                    qlf$comparison <- colnames(contrasts)[i]
                    qlf$table$p_adj <- p.adjust(qlf$table$PValue, method='fdr')
                    return(qlf$table)
                  })


comp_cols <- cond_cols; cond_cols[3] <- '#B26183'

### For combined
for (c in comps) {
  strip_df <- lapply(1:10, function(i) {
    tmp_df <- qlf.res[[paste(c[1],c[2],i,sep='_')]]
    tmp_df$bin <- diff.bin_labs[i]
    rownames(tmp_df) <- paste(rownames(tmp_df),i,sep='_')
    return(tmp_df)
  })
  strip_df <- as.data.frame(do.call(rbind, strip_df))
  
  prop_degs <- sapply(1:10, function(i) {
    x <- paste(c[1],c[2],i,sep='_')
    return(round(nrow(qlf.res[[x]][qlf.res[[x]]$p_adj<0.05,])/nrow(qlf.res[[x]]), 5))
  })
  
  pdf(file=paste0(output_dir, '/PseudoDEG_epiDiffBins_jitterFC_RLE_',c[1],'_',c[2],'.pdf'), width=8, height=6)
  par(mar=c(4.5,4.5,2,4.5))
  stripchart(logFC ~ bin, strip_df[strip_df$p_adj>=0.05,], method='jitter', ylim=c(-7,7), cex=0.5, cex.axis=0.8,
             cex.lab=1.1,pch=16, jitter=0.3, vertical=T, col='grey90', main=paste(c[1],'vs.',c[2],'DEGs'),
             cex.main=1.2, xlab='Epithelial Cell Differentiation Decile', yaxt='n')
  axis(2, at=seq(-6,6,2), labels=seq(-6,6,2), las=2, cex.axis=0.8)
  lines(x=1:10, y=prop_degs*(14)-7, col='cornflowerblue',lwd=1.5)
  points(x=1:10, y=prop_degs*(14)-7, col='cornflowerblue', pch=18, cex=0.7)
  par(new=T)
  stripchart(logFC ~ bin, strip_df[strip_df$p_adj<0.05,], method='jitter', ylim=c(-7,7), cex=0.5, 
             xaxt='n', yaxt='n', ylab='',pch=16, jitter=0.3, vertical=T,col=cond_cols[which(comps %in% list(c))])
  axis(4, at=seq(-7,7,by=2.8), labels=seq(0,100,length=6), col='cornflowerblue',lwd=0,lwd.ticks=1,
       las=2, cex.axis=0.8, col.axis='cornflowerblue')
  mtext('Percent of Genes DE', side=4, padj=5, las=3 ,cex=1.1, col='cornflowerblue')
  dev.off()
}



lrt.res <- lapply(setNames(1:ncol(contrasts), names(contrasts[1,])),
                  function(i){
                    pb_obj.b <- pb_obj[b.genes[[strsplit(colnames(contrasts)[i],'_')[[1]][3]]], , keep=FALSE] 
                    fit_lrt <- glmFit(pb_obj.b, design, robust=TRUE)
                    lrt <- glmLRT(fit_lrt, contrast=contrasts[,i])
                    lrt$comparison <- colnames(contrasts)[i]
                    lrt$table$p_adj <- p.adjust(lrt$table$PValue, method='fdr')
                    return(lrt$table)
                  })


# Alternative approach, using limma voom
limma.res <- lapply(setNames(1:ncol(contrasts), names(contrasts[1,])),
                    function(i) {
                      pb_obj.b <- pb_obj[b.genes[[strsplit(colnames(contrasts)[i],'_')[[1]][3]]], , keep=FALSE] 
                      voom_weights <- voom(pb_obj.b, design, plot=T)
                      limma_fit <- lmFit(voom_weights, design)
                      limma_results <- contrasts.fit(limma_fit, contrasts)
                      limma_results <- eBayes(limma_results)
                      r <- cbind(as.data.frame(limma_results$coefficients)[colnames(contrasts)[i]],
                                 as.data.frame(limma_results$p.value)[colnames(contrasts)[i]],
                                 p.adjust(as.data.frame(limma_results$p.value)[colnames(contrasts)[i]][[1]], method='fdr'))
                      names(r) <- c('logFC','p','p_adj')
                      return(r)})

# compute manually using Wilcoxon test and evaluate correlation

## Perform TMM normalization and convert to CPM (Counts Per Million)
#d <- calcNormFactors(d, method = "TMM")  # already did this
count_norm <- cpm(pb_obj)
count_norm <- as.data.frame(log(count_norm)+1)

# Run the Wilcoxon rank-sum test for each gene
wilcox.res <- list()
comps <- list(c('SSc','HC'), c('GERD','HC'), c('SSc', 'GERD'))
for (c in comps) {
  print(c)
  comp_results <- list()
  for (l in levels(diff.bin)) {
    print(l)
    gene.results <- lapply(1:nrow(pb_obj$counts), 
                           function(i){
                             gene.dat <- cbind.data.frame(gene=as.numeric(t(count_norm[i,])), condition, diff.bin)
                             gene.test <- wilcox.test(jitter(gene, amount=0.0000001)~condition, 
                                                      gene.dat[(gene.dat$condition %in% c) & (gene.dat$diff.bin==l),])
                             gene.stats <- data.frame(gene=rownames(count_norm)[i],
                                                      fc=mean(gene.dat[gene.dat$condition==c[1] & gene.dat$diff.bin==l,'gene']) -
                                                        mean(gene.dat[gene.dat$condition==c[2] & gene.dat$diff.bin==l,'gene']),
                                                      p=gene.test$p.value, stringsAsFactors = F)
                             return(gene.stats)
                           })
    layer_results <- as.data.frame(do.call(rbind, gene.results))
    comp_results[[paste(c[1],c[2],l,sep='_')]] <- layer_results
  }
  wilcox.res <- c(wilcox.res, comp_results)
}

cor_fc_list <- lapply(setNames(names(contrasts[1,]), names(contrasts[1,])),
                      function(c){
                        x <- cbind(qlf.res[[c]]$logFC, lrt.res[[c]]$logFC, limma.res[[c]]$logFC, wilcox.res[[c]]$fc)
                        cor_fc <- cor(x, method='spearman')
                        colnames(cor_fc) <- c('edgeR_QLF', 'edgeR_LRT','limma_voom','wilcox')
                        rownames(cor_fc) <- c('edgeR_QLF', 'edgeR_LRT','limma_voom','wilcox')
                        return(cor_fc)
                      })
cor_p_list <- lapply(setNames(names(contrasts[1,]), names(contrasts[1,])),
                     function(c){
                       x <- cbind(qlf.res[[c]]$PValue, lrt.res[[c]]$PValue, limma.res[[c]]$p, wilcox.res[[c]]$p)
                       cor_fc <- cor(x, method='spearman')
                       colnames(cor_fc) <- c('edgeR_QLF', 'edgeR_LRT','limma_voom','wilcox')
                       rownames(cor_fc) <- c('edgeR_QLF', 'edgeR_LRT','limma_voom','wilcox')
                       return(cor_fc)
                     })

save(qlf.res, lrt.res, limma.res, wilcox.res,
     cor_fc_list, cor_p_list, file=paste0(output_dir, '/pbDEG_dat_bins.RData'))

load(paste0(output_dir, '/pbDEG_dat_bins.RData'))


######### SPLIT BY LOCATION #########
pb_obj$samples$group <- paste(pb_obj$samples$condition, pb_obj$samples$location, pb_obj$samples$diff.bin, sep='_')

samples <- factor(sub("\\-.*$","", pb_obj$samples$sample))
sample_cols <- sample(alphabet(26),20)
diffbin_cols <- rev(brewer.pal(10, "Spectral"))

diff.bin <- factor(pb_obj$samples$diff.bin, levels=as.character(1:10))
condition <- factor(pb_obj$samples$condition, levels=c('SSc','GERD','HC'))
location <- factor(pb_obj$samples$location, levels=c('Distal','Proximal'))
group <- factor(pb_obj$samples$group, levels=unique(pb_obj$samples$group))


pb_obj <- pb_obj[unique(unlist(b.genes)), , keep=FALSE] 


# create design matrix
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))

# fit EdgeR model
pb_obj <- estimateDisp(pb_obj, design, robust=TRUE)
pb_obj$common.dispersion
plotBCV(pb_obj)

fit <- glmQLFit(pb_obj, design, robust=T)
plotQLDisp(fit)

fit_lrt <- glmFit(pb_obj, design, robust=TRUE)

res <- glmQLFTest(fit, coef=ncol(design)) #?
summary(decideTests(res)) #?

# Make contrasts matrix (for DEG comparisons)
contrast_comps <- c(); contrast_names <- c()
for (loc in rev(levels(location))) {
  for (c in comps) {
    for (l in 1:10) {
      contrast_names <- c(contrast_names, paste(c[1],c[2],loc,l,sep='_'))
      contrast_comps <- c(contrast_comps, paste0(paste(c[1],loc,l,sep='_'), '-', paste(c[2],loc,l, sep='_')))
    }
  }
}

contrasts <- makeContrasts(contrasts=contrast_comps, levels=design)
colnames(contrasts) <- contrast_names

# Loop through contrasts and compute DEGs, save to list
qlf.res <- lapply(setNames(1:ncol(contrasts), names(contrasts[1,])),
                  function(i){
                    print(c(i,names(contrasts[1,])[i]))
                    pb_obj.b <- pb_obj[b.genes[[strsplit(colnames(contrasts)[i],'_')[[1]][4]]], , keep=FALSE] 
                    pb_obj.b <- calcNormFactors(pb_obj.b, method='RLE')
                    #print(head(pb_obj.b$samples$norm.factors))
                    #plotBCV(pb_obj.b)
                    fit <- glmQLFit(pb_obj.b, design, robust=T)
                    plotQLDisp(fit)
                    qlf <- glmQLFTest(fit, contrast=contrasts[,i])
                    qlf$comparison <- colnames(contrasts)[i]
                    qlf$table$p_adj <- p.adjust(qlf$table$PValue, method='fdr')
                    return(qlf$table)
                  })


#### Compare DEG profiles for SSc vs. HC and GERD vs. HC for epithelial layers ####
#  (there is probably a more efficient way to do this)  

load(paste0(output_dir, '/endo_pb.RData'))


#### Pathway enrichment Analysis ####

tf_enrich.bins.res <- lapply(setNames(levels(diff.bin), diff.bin_labs), function(b) {
  print(b)
  background_genes <- b.genes[[b]]
  test_genes <- rownames(qlf.res[[paste0('SSc_HC_',b)]][(qlf.res[[paste0('SSc_HC_',b)]]$p_adj<0.05) &
                                                          (qlf.res[[paste0('SSc_GERD_',b)]]$p_adj<0.05) &
                                                          (qlf.res[[paste0('SSc_HC_',b)]]$logFC>0),])
  ego <- enricher(gene          = test_genes,
                  universe      = background_genes,
                  pAdjustMethod = "fdr",
                  pvalueCutoff  = 1,
                  TERM2GENE     = annot_tf)
  
  #View(as.data.frame(ego@result))
  if(!is.null(ego)) {
    ego_df <- as.data.frame(ego@result)
    ego_df <- ego_df[ego_df$Count>2,]
    ego_df$p.adjust <- p.adjust(ego_df$pvalue)
    ego_df <- ego_df[,c(1,3:ncol(ego_df))]
    rownames(ego_df) <- NULL
    if(nrow(ego_df)>0) {
      return(ego_df)
    }
  }
})


kegg_enrich.bins.res <- lapply(setNames(levels(diff.bin), diff.bin_labs), function(b) {
  print(b)
  background_genes <- b.genes[[b]]
  test_genes <- rownames(qlf.res[[paste0('SSc_HC_',b)]][(qlf.res[[paste0('SSc_HC_',b)]]$p_adj<0.05) &
                                                          (qlf.res[[paste0('SSc_GERD_',b)]]$p_adj<0.05),])
  ego <- enricher(gene          = test_genes,
                  universe      = background_genes,
                  pAdjustMethod = "fdr",
                  pvalueCutoff  = 1,
                  TERM2GENE     = annot_kegg)
  
  #View(as.data.frame(ego@result))
  if(!is.null(ego)) {
    ego_df <- as.data.frame(ego@result)
    ego_df <- ego_df[ego_df$Count>2,]
    ego_df$p.adjust <- p.adjust(ego_df$pvalue)
    ego_df <- ego_df[,c(1,3:ncol(ego_df))]
    rownames(ego_df) <- NULL
    if(nrow(ego_df)>0) {
      return(ego_df)
    }
  }
})

fc <- gene.df[gene.df$SYMBOL %in% test_genes, 'logFC']
names(fc) <- gene.df[gene.df$SYMBOL %in% test_genes, 'ENTREZID']
if (nrow(as.data.frame(ego)) > 200) {
  cnetplot(ego, showCategory=sum(ego@result$p.adjust<0.05), layout='kk',
           color.params=list(foldChange=fc,category=alpha("#E5C494",0.8)), 
           cex.params=list(category_label=0.4, gene_label=0.55)) + scale_size(name='N genes') +
    theme(legend.title=element_text(size=8, face='bold'), plot.title=element_text(hjust=0.5, size=10,face='bold')) + 
    ggtitle(b) + scale_colour_gradient2(name='LogFC',low='blue3',high='red3')
}
return(ego)
})

for (b in layers) {
  pdf(file=paste0(output_dir, '/PseudoDEG_enrich_GO_',b,'_cnet.pdf'), width=6, height=6)
  print(attributes(go_enrich.res[[b]])$plot)
  dev.off()
  n_sig_pathways <- nrow(as.data.frame(go_enrich.res[[b]]))
  if (n_sig_pathways>0){
    pdf(file=paste0(output_dir, '/PseudoDEG_enrich_GO_',b,'_dotplot.pdf'), width=6, height=8)
    print(dotplot(go_enrich.res[[b]], showCategory=min(n_sig_pathways,20)) + ggtitle(b) + 
            theme(axis.text.y=element_text(size=8), legend.title=element_text(size=8, face='bold'),
                  plot.title=element_text(hjust=0.5, size=10,face='bold')) + 
            scale_color_gradient(name='q-value',low='rosybrown2',high='red3'))
    dev.off()
  }
  write.table(as.data.frame(go_enrich.res[[b]]), 
              paste0(output_dir, '/PseudoDEG_enrich_GO_',b,'.tsv'), quote=F, row.names=F,col.names=T, sep="\t")
}







diff.bin_labs <- c('[0-10%]','(10-20%]','(20-30%]','(30-40%]','(40-50%]','(50-60%]','(60-70%]','(70-80%]','(80-90%]','(90-100%]')
diff.bin_labs <- factor(diff.bin_labs, levels=diff.bin_labs)
comp_cols <- lapply(setNames(levels(diff.bin), levels(diff.bin)), function(l) {
  cols <- factor(ifelse(qlf.res[[paste('SSc','HC',l,sep='_')]]$p_adj<0.05,
                        ifelse(qlf.res[[paste('GERD','HC',l,sep='_')]]$p_adj<0.05,'cornflowerblue',cond_cols[1]),
                        ifelse(qlf.res[[paste('GERD','HC',l,sep='_')]]$p_adj<0.05,cond_cols[2], 'grey90')),
                 levels=c('grey90',cond_cols[2], cond_cols[1],'cornflowerblue'))
  return(cols)
})

pdf(file=paste0(output_dir, '/PseudoDEG_epiDiffBins_scatterFC_SScGERD_HCsig.pdf'), width=12, height=6)
xrange<-c(-6,6)
yrange<-c(-6,6)
par(mfrow=c(2,5))
for(l in levels(diff.bin)) {
  plot(qlf.res[[paste('GERD','HC',l,sep='_')]][order(comp_cols[[l]]),'logFC'], 
       qlf.res[[paste('SSc','HC',l,sep='_')]][order(comp_cols[[l]]),'logFC'],
       pch=16,col=as.character(sort(comp_cols[[l]])), cex=0.8,
       ylab='SSc vs. HC (logFC)', xlab='GERD vs. HC (logFC)', main=diff.bin_labs[as.numeric(l)], xlim=xrange,ylim=yrange)
  text(-4.5,6,paste('r =',round(cor(qlf.res[[paste('GERD','HC',l,sep='_')]]$logFC, qlf.res[[paste('SSc','HC',l,sep='_')]]$logFC),2)))
  abline(0,1,lty=2, col='grey80')
}
dev.off()

comp_cols <- lapply(setNames(levels(diff.bin), levels(diff.bin)), function(l) {
  cols <- factor(ifelse(qlf.res[[paste('SSc','HC',l,sep='_')]]$p_adj<0.05 & 
                          qlf.res[[paste('GERD','HC',l,sep='_')]]$p_adj>0.05 & 
                          qlf.res[[paste('SSc','GERD',l,sep='_')]]$p_adj<0.05, cond_cols[1],
                        ifelse(qlf.res[[paste('SSc','HC',l,sep='_')]]$p_adj>0.05 & 
                                 qlf.res[[paste('GERD','HC',l,sep='_')]]$p_adj<0.05 & 
                                 qlf.res[[paste('SSc','GERD',l,sep='_')]]$p_adj<0.05, cond_cols[2],
                               'grey90')),
                 levels=c('grey90',cond_cols[2], cond_cols[1]))
  return(cols)
})

pdf(file=paste0(output_dir, '/PseudoDEG_epiDiffBins_scatterFC_SSc_GERDsig.pdf'), width=12, height=6)
xrange<-c(-6,6)
yrange<-c(-6,6)
par(mfrow=c(2,5))
for(l in levels(diff.bin)) {
  plot(qlf.res[[paste('GERD','HC',l,sep='_')]][order(comp_cols[[l]]),'logFC'], 
       qlf.res[[paste('SSc','HC',l,sep='_')]][order(comp_cols[[l]]),'logFC'],
       pch=16,col=as.character(sort(comp_cols[[l]])), cex=0.8,
       ylab='SSc vs. HC (logFC)', xlab='GERD vs. HC (logFC)', main=diff.bin_labs[as.numeric(l)], xlim=xrange,ylim=yrange)
  abline(0,1,lty=2, col='grey80')
  text(-4.5,6,paste('r =',round(cor(qlf.res[[paste('GERD','HC',l,sep='_')]]$logFC, qlf.res[[paste('SSc','HC',l,sep='_')]]$logFC),2)))
}
dev.off()
#####

comp_cols <- cond_cols; cond_cols[3] <- '#B26183'

for (c in comps) {
  for (loc in rev(levels(location))) {
    strip_df <- lapply(1:10, function(i) {
      tmp_df <- qlf.res[[paste(c[1],c[2],loc,i,sep='_')]]
      tmp_df$bin <- diff.bin_labs[i]
      rownames(tmp_df) <- paste(rownames(tmp_df),i,sep='_')
      return(tmp_df)
    })
    strip_df <- as.data.frame(do.call(rbind, strip_df))
    
    prop_degs <- sapply(1:10, function(i) {
      x <- paste(c[1],c[2],loc,i,sep='_')
      return(round(nrow(qlf.res[[x]][qlf.res[[x]]$p_adj<0.05,])/nrow(qlf.res[[x]]), 5))
    })
    
    pdf(file=paste0(output_dir, '/PseudoDEG_epiDiffBins_jitterFC_',c[1],'_',c[2],'_',loc,'.pdf'), width=8, height=6)
    par(mar=c(4.5,4.5,2,4.5))
    stripchart(logFC ~ bin, strip_df[strip_df$p_adj>=0.05,], method='jitter', ylim=c(-7,7), cex=0.4, cex.axis=0.7,
               cex.lab=0.9,pch=16, jitter=0.3, vertical=T, col='grey90', main=paste(c[1],'vs.',c[2],'DEGs -',loc),
               cex.main=0.9, xlab='Epithelial Cell Differentiation Decile', yaxt='n')
    axis(2, at=seq(-6,6,2), labels=seq(-6,6,2), las=2, cex.axis=0.7)
    lines(x=1:10, y=prop_degs*(14)-7, col='cornflowerblue',lwd=1.5)
    points(x=1:10, y=prop_degs*(14)-7, col='cornflowerblue', pch=18, cex=0.7)
    par(new=T)
    stripchart(logFC ~ bin, strip_df[strip_df$p_adj<0.05,], method='jitter', ylim=c(-7,7), cex=0.4, 
               xaxt='n', yaxt='n', ylab='',pch=16, jitter=0.3, vertical=T,col=cond_cols[which(comps %in% list(c))])
    axis(4, at=seq(-7,7,by=2.8), labels=seq(0,100,length=6), col='cornflowerblue',lwd=0,lwd.ticks=1,
         las=2, cex.axis=0.7, col.axis='cornflowerblue')
    mtext('Percent of Genes DE', side=4, padj=5, las=3 ,cex=0.9, col='cornflowerblue')
    dev.off()
    
  }
}



# Check direction of fold change convention
count_norm <- cpm(pb_obj)
count_norm <- as.data.frame(count_norm)
head(qlf.res[[3]][order(qlf.res[[3]]$logFC),],10)

genes=c('NEFL','KRT1')
#genes <- rownames(head(qlf.res[[3]][order(qlf.res[[3]]$logFC, decreasing=T),]))
par(mfrow=c(length(genes),3))
for (gene in genes) {
  for (l in layers) {
    data <- cbind.data.frame(gene=as.numeric(t(count_norm[gene,])), condition, layer)
    boxplot(gene~condition, data=data[(data$layer==l),], main=str_to_title(l),
            ylab=paste(gene, '(cpm)'),xlab='', col=scales::alpha(cond_cols, .5))
    beeswarm(gene~condition, data=data[(data$layer==l),], ylab=gene, pch=16, add=T)
  }
}


dt <- lapply(lapply(qlf.res, decideTestsDGE), summary)
dt.all <- do.call("cbind", dt)
dt.all


# Make heatmap of top genes
top <- 1000
topMarkers <- list()
lfc <- data.frame(row.names=rownames(pb_obj$counts))

for(i in 1:length(qlf.res)) {
  print(i)
  sig_genes <- rownames(qlf.res[[i]][qlf.res[[i]]$PValue<(0.05/nrow(qlf.res[[i]])),])
  #ord <- order(qlf[[i]]$table$PValue, decreasing=FALSE)
  #topMarkers[[i]] <- rownames(y)[ord[1:top]]
  topMarkers[[i]] <- sig_genes
  if (i==1){
    lfc[,names(qlf.res)[i]] <- qlf.res[[i]]$logFC
  } else {
    lfc <- cbind(lfc, qlf.res[[i]]$logFC)
    names(lfc) <- c(names(lfc)[1:(i-1)], names(qlf.res)[i])
  }
}
topMarkers <- unique(unlist(topMarkers))

lfc <- lfc[,c(1,4,2,5,3,6)]  # save fold change values for heatmap

y2 <- Seurat2PB(epi_set, sample='condition',cluster='layers')
colnames(y2$counts) <- gsub("cluster", "", colnames(y2$counts))
rownames(y2$samples) <- gsub("cluster", "", rownames(y2$samples))
lcpm <- cpm(y2, log=T)[,4:6]  # save log cpm values for heatmap (alternative to lfc)

annot <- data.frame(Condition=rep(c('SSc','GERD'),3),
                    Layer=c(rep('basal',2), rep('suprabasal',2),rep('superficial',2)))

rownames(annot) <- paste(annot$Condition,'HC', annot$Layer, sep='_')
ann_colors <- list(layer=2:4)
names(ann_colors$layer) <- paste0("layer ", levels(layer))

jpeg(paste0(output_dir, '/pseudoDEG_heatmap_SScGERD_HC.jpg'), width=7, height=9,units='in',res=150)
pheatmap::pheatmap(lfc[topMarkers,], breaks=seq(-7,7,length.out=101),
                   color=colorRampPalette(c("blue","white","red"))(100), scale="none",
                   cluster_cols=F, border_color="NA", fontsize_row=5,
                   treeheight_row=70, treeheight_col=70, cutree_cols=3,
                   clustering_method="ward.D2", show_colnames=FALSE, show_rownames = F,
                   annotation_col=annot, legend_labels=c('LogFC vs. HC'),
                   annotation_colors=list(Layer=c(basal=layer_cols_3[1], suprabasal=layer_cols_3[2],
                                                  superficial=layer_cols_3[3]),
                                          Condition=c(SSc=cond_cols[1], GERD=cond_cols[2]) ))
dev.off()



# Check direction of fold change convention
count_norm <- cpm(pb_obj)
count_norm <- as.data.frame(count_norm)
head(qlf.res[[3]][order(qlf.res[[3]]$logFC),],10)
gene='IL36A'
genes <- rownames(head(qlf.res[[3]][order(qlf.res[[3]]$logFC),]))
genes <- c('AADAC','NEFL')
par(mfrow=c(1,length(genes)))
for (gene in genes) {
  data <- cbind.data.frame(gene=as.numeric(t(count_norm[gene,])), condition, layer)
  boxplot(gene~condition, data=data[(data$condition %in% c('HC','SSc')) & 
                                      (data$layer=='superficial'),],
          ylab=gene)
  beeswarm(gene~condition, data=data[(data$condition %in% c('HC','SSc')) & 
                                       (data$layer=='superficial'),], ylab=gene, pch=16, add=T)
}

genes <- c('NEFL', 'MUC22','LCE3D','STK31','AADAC','CCL19','SLC8A1-AS1')
#layer='suprabasal'

genes <- c('HOXA9', 'PITX2','IL33','KRT1','NTS','TLE6','PRR4','HOXC8')


for (gene in genes) {
  pdf(file=paste0(output_dir, '/PseudoDEG_epiLayers_',gene,'.pdf'), width=8, height=3)
  par(mfrow=c(1,length(layers)))
  for (l in layers) {
    data <- cbind.data.frame(gene=as.numeric(t(count_norm[gene,])), condition, layer)
    boxplot(gene~condition, data=data[data$layer==l,], ylab='CPM', col=alpha(cond_cols,0.75),
            xlab='', pch=16, cex.main=1.5, ylim=c(0,max(count_norm[gene,])), outline=F)
    points(gene~jitter(as.numeric(condition), amount=0.1), data=data[data$layer==layer,],pch=16, cex=0.8)
    title(main=bquote(italic(.(gene))), line=3, cex.main=1.5, font.main=1.8)
    mtext(str_to_title(l), side=3, padj=-1, las=1 ,cex=0.8, col='grey30')
  }
  dev.off()
}


beeswarm(gene~condition, data=data[data$layer==layer,],pch=16, add=T)



pb_ssc_hc <- FindMarkers(subset(epi_aggro, layer=='superficial'), ident.1='SSc',ident.2='HC',
                         min.pct=0, logfc.threshold=0)
pb_gerd_hc<- FindMarkers(subset(epi_aggro, layer=='superficial'), ident.1='GERD',ident.2='HC',
                         min.pct=0, logfc.threshold=0)


#epi_aggro <- AggregateExpression(epi_set, group.by = c("orig.ident",'diff.bin'), 
#                                 return.seurat=TRUE, assays='RNA', slot='counts')
load(paste0(output_dir,'epi_aggro_diffbins.RData'))


Idents(epi_aggro) <- 'condition'
epi_aggro <- SetIdent(epi_aggro, value = "condition",)

library(stringr)
epi_aggro@meta.data$diff.bin <- word(rownames(epi_aggro@meta.data),2,sep='\\_')
epi_aggro@meta.data$location <- word(epi_aggro@meta.data$orig.ident,2,sep='\\-')
epi_aggro@meta.data$condition <- ifelse(str_starts(epi_aggro@meta.data$orig.ident, 'HC'),'HC',
                                        ifelse(str_starts(epi_aggro@meta.data$orig.ident,'I'),
                                               'GERD','SSc'))
epi_aggro <- NormalizeData(epi_aggro)

x <- FindMarkers(subset(epi_aggro, diff.bin=='9'), ident.1='SSc',ident.2='HC')

pdf(file=paste0(output_dir, '/Epi_IL6_IL6R_HCvSSc_pseudobulk.pdf'), width=6, height=4)
plot_grid(ggdraw() + draw_label('Epi cells (pseudobulk)', fontface='bold', size=12),
          VlnPlot(epi_aggro, features=c('IL6','IL6R'), group.by='condition', 
                  pt.size=0.4, cols=c('skyblue','#F8766D')) & 
            theme(axis.title.x = element_blank(), axis.title = element_text(size=9), 
                  axis.text =element_text(size=8), title = element_text(size=10)) & 
            geom_boxplot(width=0.1, fill="white", outlier.shape=NA),
          ncol=1, rel_heights=c(0.1, 1), scale=0.98)
dev.off()

#epi_set$diff.score <- (epi_set$diff.score-min(epi_set$diff.score)) / (max(epi_set$diff.score)-min(epi_set$diff.score))

epi_set$diff.bin <- NA
for (loc in c('Proximal','Distal')) {
  print(loc)
  epi_subset <- subset(epi_set, subset=location==loc)
  epi_subset$diff.bin <- cut(epi_subset$diff.score, 
                             breaks=unique(quantile(epi_subset$diff.score,
                                                    probs=seq.int(0,1,by=1/10))), include.lowest=T)
  epi_set@meta.data[epi_set@meta.data$location==loc, 'diff.bin'] <- epi_subset$diff.bin
}


epi_aggro <- AggregateExpression(epi_set, group.by = c("orig.ident"), 
                                 return.seurat=F, assays='RNA', slot='counts')

epi_aggro <- SetIdent(epi_aggro, value = "condition",)

epi_aggro@meta.data$orig.ident <- word(rownames(epi_aggro@meta.data),1,sep='\\_')
epi_aggro@meta.data$condition <- factor(word(rownames(epi_aggro@meta.data),2,sep='\\_'), 
                                        levels = c('SSc','GERD','HC'))
epi_aggro@meta.data$layer <- factor(word(rownames(epi_aggro@meta.data),3,sep='\\_'),
                                    levels=c('basal','suprabasal','superficial'))

pdf(file=paste0(output_dir, '/Myeloid_IL6_IL6R_HCvSSc_pseudobulk.pdf'), width=6, height=4)
plot_grid(ggdraw() + draw_label('Epithelial cells (pseudobulk)', fontface='bold', size=12),
          VlnPlot(epi_aggro, features=c('IL33','H19'), group.by='layer', split.by='condition', 
                  pt.size=0.4, cols=rep(cond_cols,10)) & 
            theme(axis.title.x = element_blank(), axis.title = element_text(size=9), 
                  axis.text =element_text(size=8), title = element_text(size=10)) & 
            geom_boxplot(width=0.1, fill="white", outlier.shape=NA),
          ncol=1, rel_heights=c(0.1, 1), scale=0.98)
dev.off()

VlnPlot(epi_aggro, features=c('IL33','H19'),  split.by='condition', group.by='layer',
        pt.size=0.4, cols=rep(cond_cols,10)) & 
  theme(axis.title.x = element_blank(), axis.title = element_text(size=9), 
        axis.text =element_text(size=8), title = element_text(size=10)) #& 
#geom_boxplot(epi_aggro@meta.data, width=0.1, fill="white", outlier.shape=NA)

epi_aggro <- SetIdent(epi_aggro, value = "condition",)
deg_epi_pseudobulk <- list()
for (l in c('basal','suprabasal','superficial')) {
  print(l)
  deg_epi_pseudobulk[[l]] <- FindMarkers(subset(epi_aggro, layer==l), ident.1='SSc', ident.2='HC', 
                                         logfc.threshold=0,min.pct=0.001)
}


deg_epi_pseudo <- deg_epi_pseudo[,c(3,4,2,1)]
names(deg_epi_pseudo) <- c('pct.Hc','pct.SSc', 'avg_log2FC', 'p_val')
deg_epi_pseudo$avg_log2FC <- deg_epi_pseudo$avg_log2FC*-1
write.table(deg_epi_pseudo, paste0(root_dir, '/Epi_pseudobulkDEG_IL6_IL6R.tsv'), quote=F, row.names=T,col.names=T, sep="\t")



#### ANALYSIS BY BIOPSY LOCATION ####
pb_obj$samples$group <- paste(pb_obj$samples$condition, pb_obj$samples$layer, pb_obj$samples$location, sep='_')

group <- factor(pb_obj$samples$group, levels=unique(pb_obj$samples$group))


# get correlations - across conditions
df <- data.frame(matrix(ncol=6, nrow=0)); i <- 1
for (comp in comps) {
  for (layer in layers_rep) {
    for (loc in rev(levels(location))) {
      corP.P <- round(cor(rowMeans(log2(pb_obj[, pb_obj$samples$condition==comp[1] & pb_obj$samples$layer==layer & 
                                                 pb_obj$samples$location==loc, keep=FALSE]$counts[ct.genes[[layer]],]+1)), 
                          rowMeans(log2(pb_obj[, pb_obj$samples$condition==comp[2] & pb_obj$samples$layer==layer & 
                                                 pb_obj$samples$location==loc, keep=FALSE]$counts[ct.genes[[layer]],]+1))),4)
      corS.P <- round(cor(rowMeans(log2(pb_obj[, pb_obj$samples$condition==comp[1] & pb_obj$samples$layer==layer & 
                                                 pb_obj$samples$location==loc, keep=FALSE]$counts[ct.genes[[layer]],]+1)), 
                          rowMeans(log2(pb_obj[, pb_obj$samples$condition==comp[2] & pb_obj$samples$layer==layer & 
                                                 pb_obj$samples$location==loc, keep=FALSE]$counts[ct.genes[[layer]],]+1)),
                          method='spearman'),4)
      df[i,] <- c(comp, layer, loc, corP.P, corS.P)
      i <- i+1
    }
  }
}
write.table(df, paste0(output_dir, '/cor_Xconds_byLayerRep.tsv'), quote=F, row.names=F,col.names=F, sep="\t")



pdf(file=paste0(output_dir, '/PseudoDEG_epi_scatter_cors_SuperDist.pdf'), width=10, height=8)
loc='Distal'
layer='superficial'
par(mfrow=c(2,2))
# ssc vs gerd
plot(rowMeans(count_norm[ct.genes[[layer]], pb_obj$samples$condition==comp[2] & pb_obj$samples$layer==layer & 
                           pb_obj$samples$location==loc]), 
     rowMeans(count_norm[ct.genes[[layer]], pb_obj$samples$condition==comp[1] & pb_obj$samples$layer==layer & 
                           pb_obj$samples$location==loc]),
     xlab='log2(counts+1)', ylab='log2(counts+1)', main='SSc vs. GERD',pch=16,col=alpha('black',0.3), cex=0.5)
abline(a=0,b=1,col='green',lty=3,lwd=1.5)
# ssc-hc vs gerd-hc
plot(rowMeans(count_norm[ct.genes[[layer]], pb_obj$samples$condition==comp[2] & pb_obj$samples$layer==layer & 
                           pb_obj$samples$location==loc])-
       rowMeans(count_norm[ct.genes[[layer]], pb_obj$samples$condition=='HC' & pb_obj$samples$layer==layer & 
                             pb_obj$samples$location==loc]), 
     rowMeans(count_norm[ct.genes[[layer]], pb_obj$samples$condition==comp[1] & pb_obj$samples$layer==layer & 
                           pb_obj$samples$location==loc])-
       rowMeans(count_norm[ct.genes[[layer]], pb_obj$samples$condition=='HC' & pb_obj$samples$layer==layer & 
                             pb_obj$samples$location==loc]),
     xlab='log2(FC)', ylab='log2(FC)', main='SSc-HC vs. GERD-HC',pch=16,col=alpha('black',0.3), cex=0.5)
abline(a=0,b=1,col='green',lty=3,lwd=1.5)
# ssc vs hc
plot(rowMeans(count_norm[ct.genes[[layer]], pb_obj$samples$condition=='HC' & pb_obj$samples$layer==layer & 
                           pb_obj$samples$location==loc]), 
     rowMeans(count_norm[ct.genes[[layer]], pb_obj$samples$condition==comp[1] & pb_obj$samples$layer==layer & 
                           pb_obj$samples$location==loc]),
     xlab='log2(counts+1)', ylab='log2(counts+1)', main='SSc vs. HC',pch=16,col=alpha('black',0.3), cex=0.5)
abline(a=0,b=1,col='green',lty=3,lwd=1.5)
# gerd vs hc
plot(rowMeans(count_norm[ct.genes[[layer]], pb_obj$samples$condition=='HC' & pb_obj$samples$layer==layer & 
                           pb_obj$samples$location==loc]), 
     rowMeans(count_norm[ct.genes[[layer]], pb_obj$samples$condition==comp[2] & pb_obj$samples$layer==layer & 
                           pb_obj$samples$location==loc]),
     xlab='log2(counts+1)', ylab='log2(counts+1)', main='GERD vs. HC',pch=16,col=alpha('black',0.3), cex=0.5)
abline(a=0,b=1,col='green',lty=3,lwd=1.5)
dev.off()

pdf(file=paste0(output_dir, '/PseudoDEG_epi_geneExpCorrs_residuals.pdf'), width=11, height=6.6)
par(mfrow=c(2,3))
for (loc in rev(levels(location))){
  for (layer in layers) {
    plot(density((rowMeans(count_norm[ct.genes[[layer]], pb_obj$samples$condition==comp[1] & pb_obj$samples$layer==layer & 
                                        pb_obj$samples$location==loc])-
                    rowMeans(count_norm[ct.genes[[layer]], pb_obj$samples$condition==comp[2] & pb_obj$samples$layer==layer & 
                                          pb_obj$samples$location==loc]))/sqrt(2)),
         main='Relative Gene Expression', xlim=c(-1.5,1.5), ylim=c(0,3.2), col='#B26183',lwd=3,
         xlab='Distance from X=Y [log2(counts+1)]')
    abline(v=0,col='grey',lty=2,lwd=2)
    par(new=T)
    plot(density((rowMeans(count_norm[ct.genes[[layer]], pb_obj$samples$condition==comp[1] & pb_obj$samples$layer==layer & 
                                        pb_obj$samples$location==loc])-
                    rowMeans(count_norm[ct.genes[[layer]], pb_obj$samples$condition=='HC' & pb_obj$samples$layer==layer & 
                                          pb_obj$samples$location==loc]))/sqrt(2)),
         main='', xlim=c(-1.5,1.5), ylim=c(0,3.2), col='#8156B3',lwd=3, xlab='')
    par(new=T)
    plot(density((rowMeans(count_norm[ct.genes[[layer]], pb_obj$samples$condition==comp[2] & pb_obj$samples$layer==layer & 
                                        pb_obj$samples$location==loc])-
                    rowMeans(count_norm[ct.genes[[layer]], pb_obj$samples$condition=='HC' & pb_obj$samples$layer==layer & 
                                          pb_obj$samples$location==loc]))/sqrt(2)),
         main='', xlim=c(-1.5,1.5), ylim=c(0,3.2), col='#E26C53',lwd=3, xlab='')
    legend(x=0.5,y=3.2,legend=c('SSc vs. HC','GERD vs. HC','SSc vs. GERD'), col=c('#8156B3','#E26C53','#B26183'),
           bty='n',lty=1, lwd=3, cex=1, seg.len=0.4, x.intersp=0.5, y.intersp=0.8)
    mtext(paste(loc, layer,sep=', '),3,padj=-0.5, cex=0.8)
  }
}
dev.off()

# get correlations - across layers
df <- data.frame(matrix(ncol=6, nrow=0)); i <- 1
for (cond in levels(condition)) {
  for (layer in layers) {
    for (otherLayer in layers[!layers %in% layer]) {
      for (loc in rev(levels(location))) {
        corP <- round(cor(rowMeans(pb_obj[, pb_obj$samples$condition==cond & pb_obj$samples$layer==layer & 
                                            pb_obj$samples$location==loc, keep=FALSE]$counts), 
                          rowMeans(pb_obj[, pb_obj$samples$condition==cond & pb_obj$samples$layer==otherLayer & 
                                            pb_obj$samples$location==loc, keep=FALSE]$counts)),4)
        corS <- round(cor(rowMeans(pb_obj[, pb_obj$samples$condition==cond & pb_obj$samples$layer==layer & 
                                            pb_obj$samples$location==loc, keep=FALSE]$counts), 
                          rowMeans(pb_obj[, pb_obj$samples$condition==cond & pb_obj$samples$layer==otherLayer & 
                                            pb_obj$samples$location==loc, keep=FALSE]$counts),
                          method='spearman'),4)
        df[i,] <- c(cond, layer, otherLayer, loc, corP, corS)
        i <- i+1
      }
    }
  }
}
write.table(df, paste0(output_dir, '/../epi_cells/cor_Xlayers_byCond.tsv'), quote=F, row.names=F,col.names=F, sep="\t")


# get correlations - across locations
df <- data.frame(matrix(ncol=4, nrow=0)); i <- 1
for (cond in levels(condition)) {
  for (layer in layers) {
    corP <- round(cor(rowMeans(pb_obj[, pb_obj$samples$condition==cond & pb_obj$samples$layer==layer & 
                                        pb_obj$samples$location=='Proximal', keep=FALSE]$counts), 
                      rowMeans(pb_obj[, pb_obj$samples$condition==cond & pb_obj$samples$layer==layer & 
                                        pb_obj$samples$location=='Distal', keep=FALSE]$counts)),4)
    corS <- round(cor(rowMeans(pb_obj[, pb_obj$samples$condition==cond & pb_obj$samples$layer==layer & 
                                        pb_obj$samples$location=='Proximal', keep=FALSE]$counts), 
                      rowMeans(pb_obj[, pb_obj$samples$condition==cond & pb_obj$samples$layer==layer & 
                                        pb_obj$samples$location=='Distal', keep=FALSE]$counts),
                      method='spearman'),4)
    df[i,] <- c(cond, layer, corP, corS)
    i <- i+1
  }
}

write.table(df, paste0(output_dir, '/../epi_cells/cor_Xloc_byCond.tsv'), quote=F, row.names=F,col.names=F, sep="\t")


# get correlations - across locations
df <- data.frame(matrix(ncol=3, nrow=0)); i <- 1
for (cond in levels(condition)) {
  corP <- round(cor(rowMeans(pb_obj[, pb_obj$samples$condition==cond & pb_obj$samples$location=='Proximal', keep=FALSE]$counts), 
                    rowMeans(pb_obj[, pb_obj$samples$condition==cond & pb_obj$samples$location=='Distal', keep=FALSE]$counts)), 4)
  corS <- round(cor(rowMeans(pb_obj[, pb_obj$samples$condition==cond & pb_obj$samples$location=='Proximal', keep=FALSE]$counts), 
                    rowMeans(pb_obj[, pb_obj$samples$condition==cond & pb_obj$samples$location=='Distal', keep=FALSE]$counts),
                    method='spearman'), 4)
  df[i,] <- c(cond, corP, corS)
  i <- i+1
}
write.table(df, paste0(output_dir, '/../epi_cells/cor_conds.tsv'), quote=F, row.names=F,col.names=F, sep="\t")


# create design matrix
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))


# fit EdgeR model
pb_obj <- estimateDisp(pb_obj, design, robust=TRUE)
#pb_obj$common.dispersion
plotBCV(pb_obj)

#fit <- glmQLFit(pb_obj, design, robust=T)
plotQLDisp(fit)

#res <- glmQLFTest(fit, coef=ncol(design)) #?
summary(decideTests(res)) #?

# Make contrasts matrix (for DEG comparisons)
contrast_comps <- c(); contrast_names <- c()
for (loc in rev(levels(location))) {
  for (c in comps) {
    for (l in layers_rep) {
      contrast_names <- c(contrast_names, paste(c[1],c[2],l,loc,sep='_'))
      contrast_comps <- c(contrast_comps, paste0(paste(c[1],l,loc,sep='_'), '-', paste(c[2],l,loc, sep='_')))
    }
  }
}

contrasts <- makeContrasts(contrasts=contrast_comps, levels=design)
colnames(contrasts) <- contrast_names

# Loop through contrasts and compute DEGs, save to list
qlf.res <- lapply(setNames(1:ncol(contrasts), names(contrasts[1,])),
                  function(i){
                    pb_obj.ct <- pb_obj[ct.genes[[strsplit(colnames(contrasts)[i],'_')[[1]][3]]], , keep.lib.sizes=FALSE] 
                    pb_obj.ct <- calcNormFactors(pb_obj.ct)
                    fit <- glmQLFit(pb_obj.ct, design, robust=T)
                    qlf <- glmQLFTest(fit, contrast=contrasts[,i])
                    qlf$comparison <- colnames(contrasts)[i]
                    qlf$table$p_adj <- p.adjust(qlf$table$PValue, method='fdr')
                    return(qlf$table)
                  })


# Alternative approach, using limma voom
limma.res <- lapply(setNames(1:ncol(contrasts), names(contrasts[1,])),
                    function(i) {
                      pb_obj.ct <- pb_obj[ct.genes[[strsplit(colnames(contrasts)[i],'_')[[1]][3]]], , keep=FALSE] 
                      voom_weights <- voom(pb_obj.ct, design, plot=T)
                      limma_fit <- lmFit(voom_weights, design)
                      limma_results <- contrasts.fit(limma_fit, contrasts)
                      limma_results <- eBayes(limma_results)
                      r <- cbind(as.data.frame(limma_results$coefficients)[colnames(contrasts)[i]],
                                 as.data.frame(limma_results$p.value)[colnames(contrasts)[i]],
                                 p.adjust(as.data.frame(limma_results$p.value)[colnames(contrasts)[i]][[1]], method='fdr'))
                      names(r) <- c('logFC','p','p_adj')
                      return(r)})


# compute manually using Wilcoxon test and evaluate correlation

# Run the Wilcoxon rank-sum test for each gene
wilcox.res <- list()
for (c in comps) {
  print(c)
  comp_results <- list()
  for (l in layers_rep) {
    print(l)
    pb_obj.ct <- pb_obj[ct.genes[[l]], , keep=FALSE] 
    count_norm <- cpm(pb_obj.ct, log=T)
    count_norm <- as.data.frame(count_norm)
    for (loc in rev(levels(location))) {
      print(loc)
      gene.results <- lapply(1:nrow(pb_obj.ct$counts), 
                             function(i){
                               gene.dat <- cbind.data.frame(gene=as.numeric(t(count_norm[i,])), condition, layer, location)
                               gene.test <- wilcox.test(jitter(gene, amount=0.0000001)~condition, 
                                                        gene.dat[(gene.dat$condition %in% c) & (gene.dat$layer==l)
                                                                 & (gene.dat$location==loc),])
                               gene.stats <- data.frame(gene=rownames(count_norm)[i],
                                                        fc=mean(gene.dat[gene.dat$condition==c[1] & gene.dat$layer==l
                                                                         & gene.dat$location==loc,'gene']) -
                                                          mean(gene.dat[gene.dat$condition==c[2] & gene.dat$layer==l
                                                                        & gene.dat$location==loc,'gene']),
                                                        p=gene.test$p.value, 
                                                        p_adj=p.adjust(gene.test$p.value, method='fdr'), stringsAsFactors = F)
                               return(gene.stats)
                             })
      layer_loc_results <- as.data.frame(do.call(rbind, gene.results))
      comp_results[[paste(c[1],c[2],l,loc,sep='_')]] <- layer_loc_results
    }
  }
  wilcox.res <- c(wilcox.res, comp_results)
}

comp_cols <- lapply(setNames(as.vector(outer(layers, levels(location),paste,sep='_')), 
                             as.vector(outer(layers, levels(location),paste,sep='_'))), function(l) {
                               cols <- factor(ifelse(qlf.res[[paste('SSc','HC',l,sep='_')]]$p_adj<0.05 & 
                                                       qlf.res[[paste('GERD','HC',l,sep='_')]]$p_adj>0.05 &
                                                       qlf.res[[paste('SSc','GERD',l,sep='_')]]$p_adj<0.05, cond_cols[1],
                                                     ifelse(qlf.res[[paste('SSc','HC',l,sep='_')]]$p_adj>0.05 & 
                                                              qlf.res[[paste('GERD','HC',l,sep='_')]]$p_adj<0.05 &
                                                              qlf.res[[paste('SSc','GERD',l,sep='_')]]$p_adj<0.05, cond_cols[2],
                                                            'grey90')),
                                              levels=c('grey90',cond_cols[2], cond_cols[1]))
                               return(cols)
                             })


comp_cols <- lapply(setNames(as.vector(outer(layers, levels(location),paste,sep='_')), 
                             as.vector(outer(layers, levels(location),paste,sep='_'))), function(l) {
                               cols <- factor(ifelse(qlf.res[[paste('SSc','GERD',l,sep='_')]]$p_adj<0.05 , cond_cols[1],
                                                     ifelse(qlf.res[[paste('GERD','HC',l,sep='_')]]$p_adj<0.05, cond_cols[2],
                                                            'grey90')),
                                              levels=c('grey90',cond_cols[2], cond_cols[1]))
                               return(cols)
                             })

comp_cols <- lapply(setNames(as.vector(outer(layers, levels(location),paste,sep='_')), 
                             as.vector(outer(layers, levels(location),paste,sep='_'))), function(l) {
                               cols <- factor(ifelse(qlf.res[[paste('SSc','HC',l,sep='_')]]$p_adj<0.05,
                                                     ifelse(qlf.res[[paste('GERD','HC',l,sep='_')]]$p_adj<0.05,'cornflowerblue',cond_cols[1]),
                                                     ifelse(qlf.res[[paste('GERD','HC',l,sep='_')]]$p_adj<0.05,cond_cols[2], 'grey90')),
                                              levels=c('grey90',cond_cols[1], cond_cols[2],'cornflowerblue'))
                               return(cols)
                             })

pdf(file=paste0(output_dir, '/PseudoDEG_epi_scatterFC_SScGERD_HCloc.pdf'), width=12, height=8)
xrange<-c(-6.7,6.7)
yrange<-c(-6.7,6.7)
par(mfrow=c(2,3))
for(loc in rev(levels(location))) {
  for (l in layers) {
    temp_data <- merge(qlf.res[[paste('GERD','HC',l,loc,sep='_')]], qlf.res[[paste('SSc','HC',l,loc,sep='_')]],
                       by='row.names')
    model <- lm(I(logFC.y)~0+logFC.x, data=temp_data)
    mod_seq <- seq(-8,8, length.out=100)
    preds <- predict(model, newdata = data.frame(logFC.x=mod_seq), interval = 'confidence')
    plot(qlf.res[[paste('GERD','HC',l,loc,sep='_')]][,'logFC'], qlf.res[[paste('SSc','HC',l,loc,sep='_')]][,'logFC'],
         type='n',ylab='', xlab='', xlim=xrange, ylim=yrange, xaxt='n',yaxt='n')
    polygon(c(rev(mod_seq), mod_seq), c(rev(preds[ ,3]), preds[ ,2]), col = alpha('lightpink',0.5), border = NA)
    abline(0,model$coefficients,lty=3, col='pink')
    par(new=T)
    plot(qlf.res[[paste('GERD','HC',l,loc,sep='_')]][order(comp_cols[[paste(l, loc,sep='_')]]),'logFC'], 
         qlf.res[[paste('SSc','HC',l,loc,sep='_')]][order(comp_cols[[paste(l, loc,sep='_')]]),'logFC'],
         pch=16,col=as.character(sort(comp_cols[[paste(l, loc,sep='_')]])), cex=0.85,
         ylab='SSc vs. HC (logFC)', xlab='GERD vs. HC (logFC)', 
         main=str_to_title(paste(loc,l)), xlim=xrange,ylim=yrange)
    abline(0,1,lty=2, col='grey80')
    text(-5.3,6,paste('r =',round(cor(qlf.res[[paste('GERD','HC',l,loc,sep='_')]]$logFC, 
                                      qlf.res[[paste('SSc','HC',l,loc,sep='_')]]$logFC, method='spearman'),2)))
    text(6.7,-5.5, paste0('m = ',round(model$coefficients,2),'\n',
                          "r\u00b2 = ", round(summary(model)$r.squared,2)),col='palevioletred3', adj=c(1,1))
  }}
dev.off()


comp_cols <- lapply(setNames(as.vector(outer(layers, levels(location),paste,sep='_')), 
                             as.vector(outer(layers, levels(location),paste,sep='_'))), function(l) {
                               cols <- factor(ifelse(qlf.res[[paste('SSc','HC',l,sep='_')]]$p_adj<0.05,
                                                     ifelse(qlf.res[[paste('GERD','HC',l,sep='_')]]$p_adj<0.05,'cornflowerblue',cond_cols[1]),
                                                     ifelse(qlf.res[[paste('GERD','HC',l,sep='_')]]$p_adj<0.05,cond_cols[2], 'grey90')),
                                              levels=c('grey90',cond_cols[1], cond_cols[2],'cornflowerblue'))
                               return(cols)
                             })
#"#8156B3" "#e26b53" "#A8A39d"
pdf(file=paste0(output_dir, '/PseudoDEG_epi_scatterFC_SSc_GERDsigLocSLOPE.pdf'), width=12, height=8)
xrange<-c(-6.7,6.7)
yrange<-c(-6.7,6.7)
par(mfrow=c(2,3))
for(loc in rev(levels(location))) {
  for (l in layers) {
    temp_data <- merge(qlf.res[[paste('GERD','HC',l,loc,sep='_')]], qlf.res[[paste('SSc','HC',l,loc,sep='_')]],
                       by='row.names')
    model <- lm(I(logFC.y)~0+logFC.x, data=temp_data)
    mod_seq <- seq(-8,8, length.out=100)
    preds <- predict(model, newdata = data.frame(logFC.x=mod_seq), interval = 'confidence')
    plot(qlf.res[[paste('GERD','HC',l,loc,sep='_')]][,'logFC'], qlf.res[[paste('SSc','HC',l,loc,sep='_')]][,'logFC'],
         type='n',ylab='', xlab='', xlim=xrange, ylim=yrange, xaxt='n',yaxt='n')
    polygon(c(rev(mod_seq), mod_seq), c(rev(preds[ ,3]), preds[ ,2]), col = alpha('lightpink',0.5), border = NA)
    abline(0,model$coefficients,lty=3, col='pink')
    par(new=T)
    plot(qlf.res[[paste('GERD','HC',l,loc,sep='_')]][order(comp_cols[[paste(l, loc,sep='_')]]),'logFC'], 
         qlf.res[[paste('SSc','HC',l,loc,sep='_')]][order(comp_cols[[paste(l, loc,sep='_')]]),'logFC'],
         pch=16,col=as.character(sort(comp_cols[[paste(l, loc,sep='_')]])), cex=0.85,
         ylab='SSc vs. HC (logFC)', xlab='GERD vs. HC (logFC)', 
         main=str_to_title(paste(loc,l)), xlim=xrange,ylim=yrange)
    abline(0,1,lty=2, col='grey80')
    text(-5.3,6,paste('r =',round(cor(qlf.res[[paste('GERD','HC',l,loc,sep='_')]]$logFC, 
                                      qlf.res[[paste('SSc','HC',l,loc,sep='_')]]$logFC, method='spearman'),2)))
    text(6.7,-5.5, paste0('m = ',round(model$coefficients,2),'\n',
                          "r\u00b2 = ", round(summary(model)$r.squared,2)),col='palevioletred3', adj=c(1,1))
  }}
dev.off()

pdf(file=paste0(output_dir, '/PseudoDEG_epi_scatterFC_SScGERD_HClocBlack.pdf'), width=12, height=8)
xrange<-c(-4,4)
yrange<-c(-4,4)
par(mfrow=c(2,3))
for(loc in rev(levels(location))) {
  for (l in layers) {
    temp_data <- merge(qlf.res[[paste('GERD','HC',l,loc,sep='_')]], qlf.res[[paste('SSc','HC',l,loc,sep='_')]],
                       by='row.names')
    model <- lm(I(logFC.y)~0+logFC.x, data=temp_data)
    mod_seq <- seq(-8,8, length.out=100)
    preds <- predict(model, newdata = data.frame(logFC.x=mod_seq), interval = 'confidence')
    plot(qlf.res[[paste('GERD','HC',l,loc,sep='_')]][,'logFC'], qlf.res[[paste('SSc','HC',l,loc,sep='_')]][,'logFC'],
         type='n',ylab='', xlab='', xlim=xrange, ylim=yrange, xaxt='n',yaxt='n')
    polygon(c(rev(mod_seq), mod_seq), c(rev(preds[ ,3]), preds[ ,2]), col = alpha('lightpink',0.5), border = NA)
    abline(0,model$coefficients,lty=3, col='pink')
    par(new=T)
    plot(qlf.res[[paste('GERD','HC',l,loc,sep='_')]][order(comp_cols[[paste(l, loc,sep='_')]]),'logFC'], 
         qlf.res[[paste('SSc','HC',l,loc,sep='_')]][order(comp_cols[[paste(l, loc,sep='_')]]),'logFC'],
         pch=16,col=alpha('midnightblue',0.5), cex=0.15,
         ylab='SSc vs. HC (logFC)', xlab='GERD vs. HC (logFC)', 
         main=str_to_title(paste(loc,l)), xlim=xrange,ylim=yrange)
    abline(0,1,lty=2, col='grey80')
    text(-5.3,6,paste('r =',round(cor(qlf.res[[paste('GERD','HC',l,loc,sep='_')]]$logFC, 
                                      qlf.res[[paste('SSc','HC',l,loc,sep='_')]]$logFC, method='spearman'),2)))
    text(6.7,-5.5, paste0('m = ',round(model$coefficients,2),'\n',
                          "r\u00b2 = ", round(summary(model)$r.squared,2)),col='palevioletred3', adj=c(1,1))
  }}
dev.off()

## Now will explore gene expression differences between Proximal and Distal regions for each condition, layer

# Make contrasts matrix (for DEG comparisons)
contrast_comps <- c(); contrast_names <- c()
for (c in levels(condition)) {
  for (l in layers) {
    contrast_names <- c(contrast_names, paste(c,l,sep='_'))
    contrast_comps <- c(contrast_comps, paste0(paste(c,l,'Proximal',sep='_'), '-', paste(c,l,'Distal', sep='_')))
  }
}

contrasts <- makeContrasts(contrasts=contrast_comps, levels=design)
colnames(contrasts) <- contrast_names

# Loop through contrasts and compute DEGs, save to list
qlf.res <- lapply(setNames(1:ncol(contrasts), names(contrasts[1,])),
                  function(i){
                    pb_obj.ct <- pb_obj[ct.genes[[strsplit(colnames(contrasts)[i],'_')[[1]][2]]], , keep=FALSE] 
                    fit <- glmQLFit(pb_obj.ct, design, robust=T)
                    qlf <- glmQLFTest(fit, contrast=contrasts[,i])
                    qlf$comparison <- colnames(contrasts)[i]
                    qlf$table$p_adj <- p.adjust(qlf$table$PValue, method='fdr')
                    return(qlf$table)
                  })

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

ssc_hc_genes <- get_genes(deg_ssc_hc_p)
enrich_ssc_hc_p <- fgseaMultilevel(fgsea_c2_sets, stats=ssc_hc_genes,
                                   minSize=10, scoreType='pos', nPermSimple=10000)


#########################################################################################################
# okay, let's test whether the surplus number of DEGs in scleroderma is due to the larger number of samples...
# I will run some bootstraps of the samples by condition with equal N's and compare counts...

deg_df_gerd_hc <- data.frame(matrix(ncol=14, nrow=0))
col_names <- c()
for (loc in rev(levels(location))) {
  for (c in comps[2]) {
    for (l in layers) {col_names <- c(col_names, paste(c[1],c[2],l,loc,'nSig',sep='_'), paste(c[1],c[2],l,loc,'absFC',sep='_'))}}}
names(deg_df_gerd_hc) <- c('n','gerd_samples', col_names)
i=1
for (n in c(1:4)) {
  pb_obj.temp <- pb_obj
  pb_obj.temp$samples$aggregate <- pb_obj.temp$samples$sample
  pb_obj.temp$samples$sample <- word(pb_obj.temp$samples$aggregate,1,sep='\\_')
  pb_obj.temp$samples$indv <- word(pb_obj.temp$samples$aggregate,1,sep='\\-')
  pb_obj.temp$samples <- pb_obj.temp$samples
  # for every combination of gerd
  gerd_samples <- unique(pb_obj.temp$samples$indv[pb_obj.temp$samples$condition=='GERD'])
  hc_samples <- unique(pb_obj.temp$samples$indv[pb_obj.temp$samples$condition=='HC'])
  gerd_combos <- combn(gerd_samples, n, simplify=F)
  for (gerd_subset in gerd_combos) {
    print(c(n, i, gerd_subset))
    pb_obj.temp.subset <- pb_obj.temp[, pb_obj.temp$samples$indv %in% c(hc_samples, gerd_subset), keep=FALSE] 
    
    pb_obj.temp.subset$samples$group <- paste(pb_obj.temp.subset$samples$condition, pb_obj.temp.subset$samples$layer, 
                                              pb_obj.temp.subset$samples$location, sep='_')
    group <- factor(pb_obj.temp.subset$samples$group, levels=unique(pb_obj.temp.subset$samples$group))
    
    design <- model.matrix(~0+group)
    colnames(design) <- gsub("group", "", colnames(design))
    
    pb_obj.temp.subset <- estimateDisp(pb_obj.temp.subset, design, robust=TRUE)
    #fit <- glmQLFit(pb_obj.temp.subset, design, robust=T)
    #res <- glmQLFTest(fit, coef=ncol(design)) #?
    
    contrast_comps <- c(); contrast_names <- c()
    if((sum(gerd_subset=='P0080')/length(gerd_subset))==1) {
      locs <- 'Distal'
    } else {
      locs <- rev(levels(location))
    }
    for (loc in locs) {
      for (c in comps[2]) {
        for (l in layers) {
          contrast_names <- c(contrast_names, paste(c[1],c[2],l,loc,sep='_'))
          contrast_comps <- c(contrast_comps, paste0(paste(c[1],l,loc,sep='_'), '-', paste(c[2],l,loc, sep='_')))
        }
      }
    }
    
    contrasts <- makeContrasts(contrasts=contrast_comps, levels=design)
    colnames(contrasts) <- contrast_names
    
    # Loop through contrasts and compute DEGs, save to list
    qlf.res <- lapply(setNames(1:ncol(contrasts), names(contrasts[1,])),
                      function(i){
                        pb_obj.ct <- pb_obj.temp.subset[ct.genes[[strsplit(colnames(contrasts)[i],'_')[[1]][3]]], , keep=FALSE] 
                        pb_obj.ct <- calcNormFactors(pb_obj.ct)
                        fit <- glmQLFit(pb_obj.ct, design, robust=T)
                        qlf <- glmQLFTest(fit, contrast=contrasts[,i])
                        qlf$comparison <- colnames(contrasts)[i]
                        qlf$table$p_adj <- p.adjust(qlf$table$PValue, method='fdr')
                        return(qlf$table)
                      })
    
    degs <- lapply(setNames(1:ncol(contrasts), names(contrasts[1,])), function(l) {
      n <- sum(qlf.res[[l]]$p_adj<0.05)
      return(n)
    })
    
    absFC <- lapply(setNames(1:ncol(contrasts), names(contrasts[1,])), function(l) {
      a <- mean(abs(qlf.res[[l]]$logFC), na.rm=T)
      return(a)
    })
    
    # store degs
    deg_df_gerd_hc[i,c('n', 'gerd_samples', col_names)] <- c(n, paste(gerd_subset, collapse='_'), unlist(rbind(degs, absFC)))
    i <- i + 1
  }
}

write.table(deg_df_gerd_hc, paste0(output_dir, '/PseudoDEG_permuteResults_GERD.tsv'), quote=F, row.names=F,col.names=T, sep="\t")


deg_df_ssc_hc <- data.frame(matrix(ncol=26, nrow=0))
col_names <- c()
for (loc in rev(levels(location))) {
  for (c in comps[1]) {
    for (l in layers) {col_names <- c(col_names, paste(l,loc,'nSig',sep='_'), 
                                      paste(l,loc,'absFC',sep='_'),
                                      paste(l,loc,'slope',sep='_'),
                                      paste(l,loc,'r2',sep='_'))}}}
names(deg_df_ssc_hc) <- c('n','ssc_samples', col_names)
i=1
for (n in c(1:10)) {
  pb_obj.temp <- pb_obj
  pb_obj.temp$samples$aggregate <- pb_obj.temp$samples$sample
  pb_obj.temp$samples$sample <- word(pb_obj.temp$samples$aggregate,1,sep='\\_')
  pb_obj.temp$samples$indv <- word(pb_obj.temp$samples$aggregate,1,sep='\\-')
  pb_obj.temp$samples <- pb_obj.temp$samples
  # for every combination of gerd
  ssc_samples <- unique(pb_obj.temp$samples$indv[pb_obj.temp$samples$condition=='SSc'])
  gerd_samples <- unique(pb_obj.temp$samples$indv[pb_obj.temp$samples$condition=='GERD'])
  hc_samples <- unique(pb_obj.temp$samples$indv[pb_obj.temp$samples$condition=='HC'])
  ssc_combos <- combn(ssc_samples, n, simplify=F)
  for (ssc_subset in ssc_combos) {
    print(c(n, i, ssc_subset))
    pb_obj.temp.subset <- pb_obj.temp[, pb_obj.temp$samples$indv %in% c(hc_samples, gerd_samples, ssc_subset), keep=FALSE] 
    
    pb_obj.temp.subset$samples$group <- paste(pb_obj.temp.subset$samples$condition, pb_obj.temp.subset$samples$layer, 
                                              pb_obj.temp.subset$samples$location, sep='_')
    
    group <- factor(pb_obj.temp.subset$samples$group, levels=unique(pb_obj.temp.subset$samples$group))
    
    design <- model.matrix(~0+group)
    colnames(design) <- gsub("group", "", colnames(design))
    
    pb_obj.temp.subset <- estimateDisp(pb_obj.temp.subset, design, robust=TRUE)
    #fit <- glmQLFit(pb_obj.temp.subset, design, robust=T)
    #res <- glmQLFTest(fit, coef=ncol(design)) #?
    
    contrast_comps <- c(); contrast_names <- c()
    if((sum(ssc_subset=='P0080')/length(ssc_subset))==1) {
      locs <- 'Distal'
    } else {
      locs <- rev(levels(location))
    }
    for (loc in locs) {
      for (c in comps[1:2]) {
        for (l in layers) {
          contrast_names <- c(contrast_names, paste(c[1],c[2],l,loc,sep='_'))
          contrast_comps <- c(contrast_comps, paste0(paste(c[1],l,loc,sep='_'), '-', paste(c[2],l,loc, sep='_')))
        }
      }
    }
    
    contrasts <- makeContrasts(contrasts=contrast_comps, levels=design)
    colnames(contrasts) <- contrast_names
    
    # Loop through contrasts and compute DEGs, save to list
    qlf.res <- lapply(setNames(1:ncol(contrasts), names(contrasts[1,])),
                      function(i){
                        pb_obj.ct <- pb_obj.temp.subset[ct.genes[[strsplit(colnames(contrasts)[i],'_')[[1]][3]]], , keep=FALSE] 
                        pb_obj.ct <- calcNormFactors(pb_obj.ct)
                        fit <- glmQLFit(pb_obj.ct, design, robust=T)
                        qlf <- glmQLFTest(fit, contrast=contrasts[,i])
                        qlf$comparison <- colnames(contrasts)[i]
                        qlf$table$p_adj <- p.adjust(qlf$table$PValue, method='fdr')
                        return(qlf$table)
                      })
    
    ssc_comps <- grep("^SSc", names(contrasts[1,]), value=T)
    
    stats <- lapply(setNames(ssc_comps, ssc_comps), function(l) {
      n <- sum(qlf.res[[l]]$p_adj<0.05)
      a <- round(mean(abs(qlf.res[[l]]$logFC), na.rm=T),2)
      lay <- word(l,3,sep='\\_')
      loc <- word(l,4,sep='\\_')
      temp_data <- merge(qlf.res[[paste('GERD','HC',lay,loc,sep='_')]], qlf.res[[l]],by='row.names')
      model <- lm(I(logFC.y)~0+logFC.x, data=temp_data)
      m <- round(model$coefficients[[1]],2)
      r2 <- round(summary(model)$r.squared,2)
      return(c(n,a,m,r2))
    })
    
    if((sum(ssc_subset=='P0080')/length(ssc_subset))==1) {
      for (d in gsub('Distal','Proximal',names(stats))) {
        stats[[d]] <- rep(NA, 4)
      }
    }
    
    # store degs
    deg_df_ssc_hc[i,c('n', 'ssc_samples', col_names)] <- c(n, paste(ssc_subset, collapse='_'), unlist(stats))
    i <- i + 1
  }
}

write.table(deg_df_ssc_hc, paste0(output_dir, '/PseudoDEG_permuteResults_SSc.tsv'), quote=F, row.names=F,col.names=T, sep="\t")

#deg_df_ssc_hc_og <- read.delim(paste0(output_dir, '/PseudoDEG_permuteResults_SSc.tsv'))

deg_df_ssc_hc$n <- as.numeric(deg_df_ssc_hc$n)
deg_df_ssc_hc[,c(3:26)] <- lapply(deg_df_ssc_hc[,c(3:26)], as.numeric)


pdf(file=paste0(output_dir, '/SSc_HC_samplePermute_Superficial.pdf'), width=6, height=8)
par(mfrow=c(2,1), mgp=c(1.8,0.5,0))
vioplot(SSc_HC_superficial_Proximal_nSig ~ n, deg_df_ssc_hc, outline=F, xlab='N SSc samples', ylab='N DEGs', col=alpha(cond_cols[1],0.25),
        main='SSc vs. HC\nSuperficial Epithelial Cells\n(Proximal)', cex.main=0.95, cex.axis=0.8, cex.lab=0.8, colMed=cond_cols[1])
vioplot(SSc_HC_superficial_Distal_nSig ~ n, deg_df_ssc_hc, outline=F, xlab='N SSc samples', ylab='N DEGs', col=alpha(cond_cols[1],0.25),
        main='SSc vs. HC\nSuperficial Epithelial Cells\n(Distal)', cex.main=0.95, cex.axis=0.8, cex.lab=0.8, colMed=cond_cols[1])
dev.off()

pdf(file=paste0(output_dir, '/SSc_HC_samplePermute_Suprabasal.pdf'), width=6, height=8)
par(mfrow=c(2,1), mgp=c(1.8,0.5,0))
vioplot(SSc_HC_suprabasal_Proximal_nSig ~ n, deg_df_ssc_hc, outline=F, xlab='N SSc samples', ylab='N DEGs', ylim=c(0,40),col=alpha(cond_cols[1],0.25),
        main='SSc vs. HC\nSuprabasal Epithelial Cells\n(Proximal)', cex.main=0.95, cex.axis=0.8, cex.lab=0.8, colMed=cond_cols[1])
vioplot(SSc_HC_suprabasal_Distal_nSig ~ n, deg_df_ssc_hc, outline=F, xlab='N SSc samples', ylab='N DEGs',ylim=c(0,40),col=alpha(cond_cols[1],0.25),
        main='SSc vs. HC\nSuprabasal Epithelial Cells\n(Distal)', cex.main=0.95, cex.axis=0.8, cex.lab=0.8, colMed=cond_cols[1])
dev.off()

pdf(file=paste0(output_dir, '/SSc_HC_samplePermute_Basal.pdf'), width=6, height=8)
par(mfrow=c(2,1), mgp=c(1.8,0.5,0))
vioplot(SSc_HC_basal_Proximal_nSig ~ n, deg_df_ssc_hc, outline=F, xlab='N SSc samples', ylab='N DEGs', ylim=c(0,30),col=alpha(cond_cols[1],0.25),
        main='SSc vs. HC\nBasal Epithelial Cells\n(Proximal)', cex.main=0.95, cex.axis=0.8, cex.lab=0.8, colMed=cond_cols[1])
vioplot(SSc_HC_basal_Distal_nSig ~ n, deg_df_ssc_hc, outline=F, xlab='N SSc samples', ylab='N DEGs', ylim=c(0,30),col=alpha(cond_cols[1],0.25),
        main='SSc vs. HC\nBasal Epithelial Cells\n(Distal)', cex.main=0.95, cex.axis=0.8, cex.lab=0.8, colMed=cond_cols[1])
dev.off()


pdf(file=paste0(output_dir, '/GERD_HC_samplePermute_Superficial.pdf'), width=4, height=8)
par(mfrow=c(2,1), mgp=c(1.8,0.5,0))
vioplot(GERD_HC_superficial_Proximal_nSig ~ n, deg_df_gerd_hc, outline=F, xlab='N GERD samples', ylab='N DEGs', ylim=c(0,600),col=alpha(cond_cols[2],0.25),
        main='GERD vs. HC\nSuperficial Epithelial Cells\n(Proximal)', cex.main=0.95, cex.axis=0.8, cex.lab=0.8, colMed=cond_cols[2])
vioplot(GERD_HC_superficial_Distal_nSig ~ n, deg_df_gerd_hc, outline=F, xlab='N GERD samples', ylab='N DEGs', ylim=c(0,600),col=alpha(cond_cols[2],0.25),
        main='GERD vs. HC\nSuperficial Epithelial Cells\n(Distal)', cex.main=0.95, cex.axis=0.8, cex.lab=0.8, colMed=cond_cols[2])
dev.off()

pdf(file=paste0(output_dir, '/GERD_HC_samplePermute_Suprabasal.pdf'), width=4, height=8)
par(mfrow=c(2,1), mgp=c(1.8,0.5,0))
vioplot(GERD_HC_suprabasal_Proximal_nSig ~ n, deg_df_gerd_hc, outline=F, xlab='N GERD samples', ylab='N DEGs', ylim=c(0,20),col=alpha(cond_cols[2],0.25),
        main='GERD vs. HC\nSuprabasal Epithelial Cells\n(Proximal)', cex.main=0.95, cex.axis=0.8, cex.lab=0.8, colMed=cond_cols[2])
vioplot(GERD_HC_suprabasal_Distal_nSig ~ n, deg_df_gerd_hc, outline=F, xlab='N GERD samples', ylab='N DEGs',ylim=c(0,20),col=alpha(cond_cols[2],0.25),
        main='GERD vs. HC\nSuprabasal Epithelial Cells\n(Distal)', cex.main=0.95, cex.axis=0.8, cex.lab=0.8, colMed=cond_cols[2])
dev.off()

pdf(file=paste0(output_dir, '/GERD_HC_samplePermute_Basal.pdf'), width=4, height=8)
par(mfrow=c(2,1), mgp=c(1.8,0.5,0))
vioplot(GERD_HC_basal_Proximal_nSig ~ n, deg_df_gerd_hc, outline=F, xlab='N GERD samples', ylab='N DEGs', ylim=c(0,12),col=alpha(cond_cols[2],0.25),
        main='GERD vs. HC\nBasal Epithelial Cells\n(Proximal)', cex.main=0.95, cex.axis=0.8, cex.lab=0.8, colMed=cond_cols[2])
vioplot(GERD_HC_basal_Distal_nSig ~ n, deg_df_gerd_hc, outline=F, xlab='N GERD samples', ylab='N DEGs', ylim=c(0,12),col=alpha(cond_cols[2],0.25),
        main='GERD vs. HC\nBasal Epithelial Cells\n(Distal)', cex.main=0.95, cex.axis=0.8, cex.lab=0.8, colMed=cond_cols[2])
dev.off()

deg_df_ssc_hc <- deg_df_ssc_hc[,c(1:14)]

deg_df_ssc_hc$n <- as.numeric(deg_df_ssc_hc$n)
deg_df_ssc_hc[,c(3:14)] <- lapply(deg_df_ssc_hc[,c(3:14)], as.numeric)

pdf(file=paste0(output_dir, '/SSc_HC_samplePermute_Superficial_FC.pdf'), width=6, height=8)
par(mfrow=c(2,1), mgp=c(1.8,0.5,0))
vioplot(SSc_HC_superficial_Proximal_absFC ~ n, deg_df_ssc_hc, outline=F, xlab='N SSc samples', ylab='Mean abs(logFC)', col=alpha(cond_cols[1],0.25),
        main='SSc vs. HC\nSuperficial Epithelial Cells\n(Proximal)', cex.main=0.95, cex.axis=0.8, cex.lab=0.8, colMed=cond_cols[1], ylim=c(0,3))
vioplot(SSc_HC_superficial_Distal_absFC ~ n, deg_df_ssc_hc, outline=F, xlab='N SSc samples', ylab='Mean abs(logFC)', col=alpha(cond_cols[1],0.25),
        main='SSc vs. HC\nSuperficial Epithelial Cells\n(Distal)', cex.main=0.95, cex.axis=0.8, cex.lab=0.8, colMed=cond_cols[1], ylim=c(0,3))
dev.off()

pdf(file=paste0(output_dir, '/SSc_HC_samplePermute_Suprabasal_FC.pdf'), width=6, height=8)
par(mfrow=c(2,1), mgp=c(1.8,0.5,0))
vioplot(SSc_HC_suprabasal_Proximal_absFC ~ n, deg_df_ssc_hc, outline=F, xlab='N SSc samples', ylab='Mean abs(logFC)', ylim=c(0,1),col=alpha(cond_cols[1],0.25),
        main='SSc vs. HC\nSuprabasal Epithelial Cells\n(Proximal)', cex.main=0.95, cex.axis=0.8, cex.lab=0.8, colMed=cond_cols[1])
vioplot(SSc_HC_suprabasal_Distal_absFC ~ n, deg_df_ssc_hc, outline=F, xlab='N SSc samples', ylab='Mean abs(logFC)',ylim=c(0,1),col=alpha(cond_cols[1],0.25),
        main='SSc vs. HC\nSuprabasal Epithelial Cells\n(Distal)', cex.main=0.95, cex.axis=0.8, cex.lab=0.8, colMed=cond_cols[1])
dev.off()

pdf(file=paste0(output_dir, '/SSc_HC_samplePermute_Basal_FC.pdf'), width=6, height=8)
par(mfrow=c(2,1), mgp=c(1.8,0.5,0))
vioplot(SSc_HC_basal_Proximal_absFC ~ n, deg_df_ssc_hc, outline=F, xlab='N SSc samples', ylab='Mean abs(logFC)', ylim=c(0,1),col=alpha(cond_cols[1],0.25),
        main='SSc vs. HC\nBasal Epithelial Cells\n(Proximal)', cex.main=0.95, cex.axis=0.8, cex.lab=0.8, colMed=cond_cols[1])
vioplot(SSc_HC_basal_Distal_absFC ~ n, deg_df_ssc_hc, outline=F, xlab='N SSc samples', ylab='Mean abs(logFC)', ylim=c(0,1),col=alpha(cond_cols[1],0.25),
        main='SSc vs. HC\nBasal Epithelial Cells\n(Distal)', cex.main=0.95, cex.axis=0.8, cex.lab=0.8, colMed=cond_cols[1])
dev.off()


pdf(file=paste0(output_dir, '/SSc_HC_samplePermute_Superficial_slope.pdf'), width=6, height=8)
par(mfrow=c(2,1), mgp=c(1.8,0.5,0))
vioplot(superficial_Proximal_slope ~ n, deg_df_ssc_hc, outline=F, xlab='N SSc samples', ylab='Slope (SSc-HC)/(GERD-HC)', col=alpha(cond_cols[1],0.25),
        main='SSc vs. HC\nSuperficial Epithelial Cells\n(Proximal)', cex.main=0.95, cex.axis=0.8, cex.lab=0.8, colMed=cond_cols[1], ylim=c(-2,2))
vioplot(superficial_Distal_slope ~ n, deg_df_ssc_hc, outline=F, xlab='N SSc samples', ylab='Slope (SSc-HC)/(GERD-HC)', col=alpha(cond_cols[1],0.25),
        main='SSc vs. HC\nSuperficial Epithelial Cells\n(Distal)', cex.main=0.95, cex.axis=0.8, cex.lab=0.8, colMed=cond_cols[1], ylim=c(-2,2))
dev.off()


pdf(file=paste0(output_dir, '/SSc_HC_samplePermute_Superficial_r2.pdf'), width=6, height=8)
par(mfrow=c(2,1), mgp=c(1.8,0.5,0))
vioplot(superficial_Proximal_r2 ~ n, deg_df_ssc_hc, outline=F, xlab='N SSc samples', ylab='r^2 (SSc-HC)/(GERD-HC)', col=alpha(cond_cols[1],0.25),
        main='SSc vs. HC\nSuperficial Epithelial Cells\n(Proximal)', cex.main=0.95, cex.axis=0.8, cex.lab=0.8, colMed=cond_cols[1], ylim=c(0,1))
vioplot(superficial_Distal_r2 ~ n, deg_df_ssc_hc, outline=F, xlab='N SSc samples', ylab='r^2 (SSc-HC)/(GERD-HC)', col=alpha(cond_cols[1],0.25),
        main='SSc vs. HC\nSuperficial Epithelial Cells\n(Distal)', cex.main=0.95, cex.axis=0.8, cex.lab=0.8, colMed=cond_cols[1], ylim=c(0,1))
dev.off()


pdf(file=paste0(output_dir, '/GERD_HC_samplePermute_Superficial_FC.pdf'), width=4, height=8)
par(mfrow=c(2,1), mgp=c(1.8,0.5,0))
vioplot(GERD_HC_superficial_Proximal_absFC ~ n, deg_df_gerd_hc, outline=F, xlab='N GERD samples', ylab='Mean abs(logFC)', col=alpha(cond_cols[2],0.25),
        main='GERD vs. HC\nSuperficial Epithelial Cells\n(Proximal)', cex.main=0.95, cex.axis=0.8, cex.lab=0.8, colMed=cond_cols[2], ylim=c(0,3))
vioplot(GERD_HC_superficial_Distal_absFC ~ n, deg_df_gerd_hc, outline=F, xlab='N GERD samples', ylab='Mean abs(logFC)', col=alpha(cond_cols[2],0.25),
        main='GERD vs. HC\nSuperficial Epithelial Cells\n(Distal)', cex.main=0.95, cex.axis=0.8, cex.lab=0.8, colMed=cond_cols[2], ylim=c(0,3))
dev.off()

pdf(file=paste0(output_dir, '/GERD_HC_samplePermute_Suprabasal_FC.pdf'), width=4, height=8)
par(mfrow=c(2,1), mgp=c(1.8,0.5,0))
vioplot(GERD_HC_suprabasal_Proximal_absFC ~ n, deg_df_gerd_hc, outline=F, xlab='N GERD samples', ylab='Mean abs(logFC)', ylim=c(0,1),col=alpha(cond_cols[2],0.25),
        main='GERD vs. HC\nSuprabasal Epithelial Cells\n(Proximal)', cex.main=0.95, cex.axis=0.8, cex.lab=0.8, colMed=cond_cols[2])
vioplot(GERD_HC_suprabasal_Distal_absFC ~ n, deg_df_gerd_hc, outline=F, xlab='N GERD samples', ylab='Mean abs(logFC)',ylim=c(0,1),col=alpha(cond_cols[2],0.25),
        main='GERD vs. HC\nSuprabasal Epithelial Cells\n(Distal)', cex.main=0.95, cex.axis=0.8, cex.lab=0.8, colMed=cond_cols[2])
dev.off()

pdf(file=paste0(output_dir, '/GERD_HC_samplePermute_Basal_FC.pdf'), width=4, height=8)
par(mfrow=c(2,1), mgp=c(1.8,0.5,0))
vioplot(GERD_HC_basal_Proximal_absFC ~ n, deg_df_gerd_hc, outline=F, xlab='N GERD samples', ylab='Mean abs(logFC)', ylim=c(0,1),col=alpha(cond_cols[2],0.25),
        main='GERD vs. HC\nBasal Epithelial Cells\n(Proximal)', cex.main=0.95, cex.axis=0.8, cex.lab=0.8, colMed=cond_cols[2])
vioplot(GERD_HC_basal_Distal_absFC ~ n, deg_df_gerd_hc, outline=F, xlab='N GERD samples', ylab='Mean abs(logFC)', ylim=c(0,1),col=alpha(cond_cols[2],0.25),
        main='GERD vs. HC\nBasal Epithelial Cells\n(Distal)', cex.main=0.95, cex.axis=0.8, cex.lab=0.8, colMed=cond_cols[2])
dev.off()



#########################


go_enrich.up.res <- lapply(setNames(c('superficial_Proximal','superficial_Distal'),
                                    c('superficial_Proximal','superficial_Distal')), function(ct) {
                                      background_genes <- ct.genes[['superficial']]
                                      background_genes <- background_genes[!grepl('^MRP',background_genes)]
                                      background_genes <- background_genes[!grepl('^RP',background_genes)]
                                      
                                      test_genes <- rownames(qlf.res[[paste0('SSc_HC_',ct)]][(qlf.res[[paste0('SSc_HC_',ct)]]$p_adj<0.05) &
                                                                                               (qlf.res[[paste0('SSc_HC_',ct)]]$logFC)<200,])
                                      test_genes <- test_genes[!grepl('^MRP',test_genes)]
                                      test_genes <- test_genes[!grepl('^RP',test_genes)]
                                      
                                      ego <- enricher(gene          = test_genes,
                                                      universe      = background_genes,
                                                      pAdjustMethod = "fdr",
                                                      pvalueCutoff  = 1,
                                                      TERM2GENE     = annot_go)
                                      
                                      ego_df <- as.data.frame(ego@result)
                                      ego_df <- ego_df[ego_df$Count>2,]
                                      ego_df$p.adjust <- p.adjust(ego_df$pvalue)
                                      ego_df <- ego_df[,c(1,3:ncol(ego_df))]
                                      rownames(ego_df) <- NULL
                                      View(ego_df)
                                      return(ego_df[ego_df$p.adjust<0.05,])
                                    })

kegg_enrich.up.res <- lapply(setNames(c('superficial_Proximal','superficial_Distal'),
                                      c('superficial_Proximal','superficial_Distal')), function(ct) {
                                        background_genes <- ct.genes[['superficial']]
                                        background_genes <- background_genes[!grepl('^MRP',background_genes)]
                                        background_genes <- background_genes[!grepl('^RP',background_genes)]
                                        
                                        test_genes <- rownames(qlf.res[[paste0('SSc_HC_',ct)]][(qlf.res[[paste0('SSc_HC_',ct)]]$p_adj<0.05) &
                                                                                                 (qlf.res[[paste0('SSc_HC_',ct)]]$logFC>0),])
                                        test_genes <- test_genes[!grepl('^MRP',test_genes)]
                                        test_genes <- test_genes[!grepl('^RP',test_genes)]
                                        
                                        ego <- enricher(gene          = test_genes,
                                                        universe      = background_genes,
                                                        pAdjustMethod = "fdr",
                                                        pvalueCutoff  = 1,
                                                        TERM2GENE     = annot_kegg)
                                        
                                        #View(as.data.frame(ego@result))
                                        ego_df <- as.data.frame(ego@result)
                                        ego_df <- ego_df[ego_df$Count>2,]
                                        ego_df$p.adjust <- p.adjust(ego_df$pvalue)
                                        ego_df <- ego_df[,c(1,3:ncol(ego_df))]
                                        rownames(ego_df) <- NULL
                                        #View(ego_df)
                                        return(ego_df[ego_df$p.adjust<0.05,])
                                      })


annot_tf$V1 <- sapply(annot_tf$V1, function(x) {strsplit(x,' ')[[1]][1]})



tf_enrich.up.res <- lapply(setNames(c('superficial_Proximal','superficial_Distal'),
                                    c('superficial_Proximal','superficial_Distal')), 
                           function(ct) {
                             background_genes <- ct.genes[['superficial']]
                             background_genes <- background_genes[!grepl('^MRP',background_genes)]
                             background_genes <- background_genes[!grepl('^RP',background_genes)]
                             
                             test_genes <- rownames(qlf.res[[paste0('SSc_HC_',ct)]][(qlf.res[[paste0('SSc_HC_',ct)]]$p_adj<0.05) &
                                                                                      (qlf.res[[paste0('SSc_HC_',ct)]]$logFC)>0,])
                             test_genes <- test_genes[!grepl('^MRP',test_genes)]
                             test_genes <- test_genes[!grepl('^RP',test_genes)]
                             ego <- enricher(gene          = test_genes,
                                             universe      = background_genes,
                                             pAdjustMethod = "fdr",
                                             pvalueCutoff  = 1,
                                             TERM2GENE     = annot_tf[annot_tf$V1 %in% background_genes,])
                             
                             #View(as.data.frame(ego@result))
                             ego_df <- as.data.frame(ego@result)
                             ego_df <- ego_df[ego_df$Count>2,]
                             ego_df$p.adjust <- p.adjust(ego_df$pvalue)
                             ego_df <- ego_df[,c(1,3:ncol(ego_df))]
                             rownames(ego_df) <- NULL
                             #View(ego_df)
                             return(ego_df)
                           })


tf_enrich.up.res <- lapply(setNames(c('superficial_Proximal','superficial_Distal'),
                                    c('superficial_Proximal','superficial_Distal')), 
                           function(ct) {
                             background_genes <- ct.genes[['superficial']]
                             background_genes <- background_genes[!grepl('^MRP',background_genes)]
                             background_genes <- background_genes[!grepl('^RP',background_genes)]
                             
                             test_genes <- rownames(qlf.res[[paste0('SSc_HC_',ct)]][(qlf.res[[paste0('SSc_HC_',ct)]]$p_adj<0.05) &
                                                                                      (qlf.res[[paste0('SSc_HC_',ct)]]$logFC)<100,])
                             test_genes <- test_genes[!grepl('^MRP',test_genes)]
                             test_genes <- test_genes[!grepl('^RP',test_genes)]
                             ego <- enricher(gene          = test_genes,
                                             universe      = background_genes,
                                             pAdjustMethod = "fdr",
                                             pvalueCutoff  = 1,
                                             TERM2GENE     = annot_tf[annot_tf$V1 %in% background_genes,])
                             
                             #View(as.data.frame(ego@result))
                             ego_df <- as.data.frame(ego@result)
                             ego_df <- ego_df[ego_df$Count>2,]
                             ego_df$p.adjust <- p.adjust(ego_df$pvalue)
                             ego_df <- ego_df[,c(1,3:ncol(ego_df))]
                             rownames(ego_df) <- NULL
                             #View(ego_df)
                             return(ego_df)
                           })

tf_enrich.down.res <- lapply(setNames(c('superficial_Proximal','superficial_Distal'),
                                      c('superficial_Proximal','superficial_Distal')), 
                             function(ct) {
                               background_genes <- ct.genes[['superficial']]
                               background_genes <- background_genes[!grepl('^MRP',background_genes)]
                               background_genes <- background_genes[!grepl('^RP',background_genes)]
                               
                               test_genes <- rownames(qlf.res[[paste0('SSc_HC_',ct)]][(qlf.res[[paste0('SSc_HC_',ct)]]$p_adj<0.05) &
                                                                                        (qlf.res[[paste0('SSc_HC_',ct)]]$logFC)<0,])
                               test_genes <- test_genes[!grepl('^MRP',test_genes)]
                               test_genes <- test_genes[!grepl('^RP',test_genes)]
                               ego <- enricher(gene          = test_genes,
                                               universe      = background_genes,
                                               pAdjustMethod = "fdr",
                                               pvalueCutoff  = 1,
                                               TERM2GENE     = annot_tf[annot_tf$V1 %in% background_genes,])
                               
                               #View(as.data.frame(ego@result))
                               ego_df <- as.data.frame(ego@result)
                               ego_df <- ego_df[ego_df$Count>2,]
                               ego_df$p.adjust <- p.adjust(ego_df$pvalue)
                               ego_df <- ego_df[,c(1,3:ncol(ego_df))]
                               rownames(ego_df) <- NULL
                               #View(ego_df)
                               return(ego_df)
                             })




for (ct in layers) {
  n_sig_pathways <- nrow(as.data.frame(tf_enrich.up.res[[ct]]))
  if (n_sig_pathways>0){
    pdf(file=paste0(output_dir, '/PseudoDEG_enrich_GO_',ct,'_dotplot.pdf'), width=6, height=8)
    print(dotplot(tf_enrich.up.res[[ct]], showCategory=min(n_sig_pathways,20)) + ggtitle(ct) + 
            theme(axis.text.y=element_text(size=8), legend.title=element_text(size=8, face='bold'),
                  plot.title=element_text(hjust=0.5, size=10,face='bold')) + 
            scale_color_gradient(name='q-value',low='rosybrown2',high='blue3'))
    dev.off()
  }
  write.table(as.data.frame(go_enrich.res[[ct]]), 
              paste0(output_dir, '/PseudoDEG_enrich_GO_',ct,'.tsv'), quote=F, row.names=F,col.names=T, sep="\t")
}





#### Diff.bins by location
pb_obj$samples$group <- paste(pb_obj$samples$condition, pb_obj$samples$diff.bin, pb_obj$samples$location, sep='_')
group <- factor(pb_obj$samples$group, levels=unique(pb_obj$samples$group))

diffbin_cols <- rev(brewer.pal(10, "Spectral"))
diff.bin <- factor(pb_obj$samples$diff.bin, levels=as.character(1:10))
diff.bin_labs <- c('[0-10%]','(10-20%]','(20-30%]','(30-40%]','(40-50%]','(50-60%]','(60-70%]','(70-80%]','(80-90%]','(90-100%]')
diff.bin_labs <- factor(diff.bin_labs, levels=diff.bin_labs)

# create design matrix
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))

# fit EdgeR model
pb_obj <- estimateDisp(pb_obj, design, robust=TRUE)
pb_obj$common.dispersion
plotBCV(pb_obj)

fit <- glmQLFit(pb_obj, design, robust=T)
plotQLDisp(fit)

# Make contrasts matrix (for DEG comparisons)
contrast_comps <- c(); contrast_names <- c()
for (loc in rev(levels(location))) {
  for (c in comps) {
    for (b in levels(diff.bin)) {
      contrast_names <- c(contrast_names, paste(c[1],c[2],b,loc,sep='_'))
      contrast_comps <- c(contrast_comps, paste0(paste(c[1],b,loc,sep='_'), '-', paste(c[2],b,loc, sep='_')))
    }
  }
}

contrasts <- makeContrasts(contrasts=contrast_comps, levels=design)
colnames(contrasts) <- contrast_names

# Loop through contrasts and compute DEGs, save to list
qlf.res <- lapply(setNames(1:ncol(contrasts), names(contrasts[1,])),
                  function(i){
                    print(i)
                    pb_obj.b <- pb_obj[b.genes[[strsplit(colnames(contrasts)[i],'_')[[1]][3]]], , keep=FALSE] 
                    fit <- glmQLFit(pb_obj.b, design, robust=T)
                    qlf <- glmQLFTest(fit, contrast=contrasts[,i])
                    qlf$comparison <- colnames(contrasts)[i]
                    qlf$table$p_adj <- p.adjust(qlf$table$PValue, method='fdr')
                    return(qlf$table)
                  })

comp_cols <- cond_cols; cond_cols[3] <- '#B26183'
for (c in comps) {
  for (loc in rev(levels(location))) {
    strip_df <- lapply(1:10, function(i) {
      tmp_df <- qlf.res[[paste(c[1],c[2],i,loc,sep='_')]]
      tmp_df$bin <- diff.bin_labs[i]
      rownames(tmp_df) <- paste(rownames(tmp_df),i,sep='_')
      return(tmp_df)
    })
    
    strip_df <- as.data.frame(do.call(rbind, strip_df))
    
    prop_degs <- sapply(1:10, function(i) {
      x <- paste(c[1],c[2],i,loc,sep='_')
      return(round(nrow(qlf.res[[x]][qlf.res[[x]]$p_adj<0.05,])/nrow(qlf.res[[x]]), 5))
    })
    
    pdf(file=paste0(output_dir, '/PseudoDEG_epiDiffBins_jitterFC_',c[1],'_',c[2],'.',loc,'.pdf'), width=8, height=6)
    par(mar=c(4.5,4.5,2,4.5))
    stripchart(logFC ~ bin, strip_df[strip_df$p_adj>=0.05,], method='jitter', ylim=c(-7,7), 
               cex.lab=0.9,pch=16, jitter=0.3, vertical=T, col='grey90', cex=0.4, cex.axis=0.7,
               main=paste(c[1],'vs.',c[2],'-',loc,'- DEGs'),
               cex.main=0.9, xlab='Epithelial Cell Differentiation Decile', yaxt='n')
    axis(2, at=seq(-6,6,2), labels=seq(-6,6,2), las=2, cex.axis=0.7)
    lines(x=1:10, y=prop_degs*(14)-7, col='cornflowerblue',lwd=1.5)
    points(x=1:10, y=prop_degs*(14)-7, col='cornflowerblue', pch=18, cex=0.7)
    par(new=T)
    stripchart(logFC ~ bin, strip_df[strip_df$p_adj<0.05,], method='jitter', ylim=c(-7,7), cex=0.4, 
               xaxt='n', yaxt='n', ylab='',pch=16, jitter=0.3, vertical=T,col=cond_cols[which(comps %in% list(c))])
    axis(4, at=seq(-7,7,by=2.8), labels=seq(0,100,length=6), col='cornflowerblue',lwd=0,lwd.ticks=1,
         las=2, cex.axis=0.7, col.axis='cornflowerblue')
    mtext('Percent of Genes DE', side=4, padj=5, las=3 ,cex=0.9, col='cornflowerblue')
    dev.off()
  }
}









genes <- c('HOXA9', 'PITX2','IL33','KRT1','NTS','TLE6','PRR4','HOXC8')

genes <- 'TNNI3'

count_norm <- cpm(pb_obj)
count_norm <- as.data.frame(log(count_norm+1))

genes=c('CCL19','IL33','CDKN2A','CCL2','PEG3','LCE3D')

genes=c("MT1A","MT1E","MT1F","MT1G","MT1H","MT1M","MT1X")

genes=c("BHLHE40","TAF7","SMAD4","FOXM1")
genes <- c('FOXM1')

genes <- c('CDC20','CDKN2C','GPSM2','NUF2','UBE2T','NUCKS1')
genes <- c('CENPF','LBR','CKAP2L','CCNA2','CCNB1','LMNB1')

foxm1_target_genes <- c('CDC20', 'CDKN2C', 'GPSM2', 'NUF2', 'UBE2T', 'NUCKS1', 'CENPF', 'LBR', 'CKAP2L', 
                        'CCNA2', 'CCNB1', 'LMNB1', 'KIF20A', 'PTTG1', 'MDC1', 'RHEB', 'KCTD9', 'CDCA2', 
                        'TUBB4B', 'CDK1', 'KIF20B', 'KIF11', 'CEP55', 'MKI67', 'CKAP5', 'TROAP', 'RACGAP1', 
                        'NUP37', 'PARPBP', 'HSP90B1', 'CKAP2', 'MZT1', 'BORA', 'CDKN3', 'KNSTRN', 'CCNB2', 
                        'PRC1', 'PLK1', 'FZR1', 'CDC25B', 'TPX2', 'ASXL1', 'UBE2C', 'GTSE1', 'NEK2')

for (i in seq(1,length(foxm1_target_genes),6)) {
  genes <- foxm1_target_genes[c(i:(i+5))]
  genes <- genes[!is.na(genes)]
  print(genes)
  pdf(file=paste0(output_dir, '/PseudoDEG_FOXM1_targets_',i,'.pdf'), width=14, height=8.5)
  par(mfrow=c(2,length(genes)), mar=c(2.5,2,3.5,1), oma=c(0,2.5,0,0))
  for (loc in rev(levels(location))) {
    for (gene in genes) {
      data <- cbind.data.frame(gene=as.numeric(t(count_norm[gene,])), condition, layer, location)
      data$condition <- factor(data$condition, c('HC','GERD','SSc'))
      boxplot(gene~condition, data=data[data$layer=='superficial' & data$location==loc,], ylab='', col=alpha(rev(cond_cols),0.75),
              xlab='', pch=16, cex.main=1.5,  outline=F)
      points(gene~jitter(as.numeric(condition), amount=0.1), data=data[data$layer==l & data$location==loc,],pch=16, cex=0.8)
      if (loc=='Proximal') {title(main=bquote(italic(.(gene))), line=3, cex.main=1.5, font.main=1.8)}
      if (gene==genes[1]) {mtext('log(CPM+1)', 2, line=2.8,cex=0.9)}
      mtext(str_to_title(loc), side=3, padj=-1, las=1 ,cex=0.8, col='grey30')
      text(.15 * par('usr')[2], .95 * par('usr')[4], adj=0,
           labels = paste("SSc vs. HC P=", round(qlf.res[[paste0('SSc_HC_','superficial_',loc)]][gene,'PValue'],4)))
    }
  }
  dev.off()
}




### Correlation with motility
load(paste0(output_dir,'epi_aggro.RData'))
load(paste0(output_dir, '/pbDEG_dat.RData'))

diff.genes <- unique(c(rownames(qlf.res[[1]][qlf.res[[1]]$p_adj<0.05 & abs(qlf.res[[1]]$logFC)>0,]),
                       rownames(qlf.res[[2]][qlf.res[[2]]$p_adj<0.05 & abs(qlf.res[[2]]$logFC)>0,]),
                       rownames(qlf.res[[3]][qlf.res[[3]]$p_adj<0.05 & abs(qlf.res[[3]]$logFC)>0,]),
                       rownames(qlf.res[[10]][qlf.res[[10]]$p_adj<0.05 & abs(qlf.res[[1]]$logFC)>0,]),
                       rownames(qlf.res[[11]][qlf.res[[11]]$p_adj<0.05 & abs(qlf.res[[2]]$logFC)>0,]),
                       rownames(qlf.res[[12]][qlf.res[[12]]$p_adj<0.05 & abs(qlf.res[[3]]$logFC)>0,])))

pb_obj <- DGEList(epi_aggro@assays$RNA@counts)  # create DGEList
pb_obj <- calcNormFactors(pb_obj)

# add metadata to psuedobulk object
pb_obj$samples$sample <- rownames(pb_obj$samples)
pb_obj$samples$condition <- ifelse(grepl('^I', pb_obj$samples$sample), 'GERD',
                                   ifelse(grepl('^P', pb_obj$samples$sample), 'SSc', 'HC'))
pb_obj$samples$layer <- word(rownames(pb_obj$samples),2,sep='\\_')
pb_obj$samples$location <- ifelse(grepl('D_', pb_obj$samples$sample), 'Distal','Proximal')

# filter non-expressed genes by subtype (logCPM >0 in >50%)
layers <- c('basal','suprabasal','superficial')
layer <- factor(pb_obj$samples$layer, levels=c('basal','suprabasal','superficial'))
condition <- factor(pb_obj$samples$condition, levels=c('SSc','GERD','HC'))
location <- factor(pb_obj$samples$location, levels=c('Distal','Proximal'))

#pb_obj <- pb_obj[unique(unlist(diff.genes)), , keep=FALSE] 

clin_dat <- read.delim(paste0(output_dir,'../clinical_data.txt'), na.strings=c('NA','NULL','null',''))

pb_obj$samples$indv <- sapply(pb_obj$samples$sample, function(x) {strsplit(x,'-')[[1]][1]})

pb_obj$samples <- merge(pb_obj$samples, clin_dat[,c('Participant_ID','median_irp','mean_DCI','pressure_60ml',
                                                    'egd_di_60ml','s_basal_ee_egjp','egj_cont_index',
                                                    'cc_v4_class_hrm','cc_v4_class_peristalsis',
                                                    'flip_contract_pattern_v2','NM_model_motility')],
                        by.x='indv', by.y='Participant_ID')

pb_obj$samples$cc_v4_class_hrm <- toupper(pb_obj$samples$cc_v4_class_hrm)
library(plyr)
hrm_vals <- c('Normal','Ineffective','Absent')
flip_vals <- c('Normal/Diminished','Normal/Diminished', 'Impaired','Absent')
nm_vals <- c('Normal','Ineffective', 'Absent')

pb_obj$samples$cc_v4_class_peristalsis[pb_obj$samples$cc_v4_class_peristalsis=='premature'] <- NA
pb_obj$samples$cc_v4_class_peristalsis <- mapvalues(pb_obj$samples$cc_v4_class_peristalsis, 
                                                    c('normal', 'ineffective', 'absent'), hrm_vals)
pb_obj$samples$flip_contract_pattern_v2 <- sub(" CONTRACTILE RESPONSE","", 
                                               pb_obj$samples$flip_contract_pattern_v2)
pb_obj$samples$flip_contract_pattern_v2 <- mapvalues(pb_obj$samples$flip_contract_pattern_v2, 
                                                     c('NORMAL','BORDERLINE/DIMINISHED', 'IMPAIRED/DISORDERED',
                                                       'ABSENT'), flip_vals)
pb_obj$samples$NM_model_motility <- mapvalues(pb_obj$samples$NM_model_motility, 
                                              c('Normal','Stage I  - Ineffective', 'Stage II - Ineffective',
                                                'Stage III - Absent'), nm_vals[c(1,2,2,3)])

diffScore_dat <- read.delim(paste0(output_dir,'../epi_cells/diffScoreStats.tsv'))

pb_obj$samples$sample <- sapply(pb_obj$samples$sample, function(x) {strsplit(x,'_')[[1]][1]})

pb_obj$samples <- merge(pb_obj$samples, diffScore_dat, by.x='sample', by.y='orig.ident')

#pb_obj <- pb_obj[, pb_obj$samples$condition=='SSc', keep=FALSE] 

pheno_flip <- factor(pb_obj$samples$flip_contract_pattern_v2, ordered=T, levels=flip_vals[2:4])
pheno_hrm <- factor(pb_obj$samples$cc_v4_class_peristalsis, ordered=T, levels=hrm_vals)
pheno_nm <- factor(pb_obj$samples$NM_model_motility, ordered=T, levels=nm_vals)

layer <- factor(pb_obj$samples$layer, levels=c('basal','suprabasal','superficial'))
condition <- factor(pb_obj$samples$condition, levels=c('SSc','GERD','HC'))
location <- factor(pb_obj$samples$location, levels=c('Distal','Proximal'))

library(MASS)
library(rms)
# Run an ordered logistic regression on each gene
library(ACAT)

ord_reg.res <- list()
phenos <- list(flip=pheno_flip, hrm=pheno_hrm, nm=pheno_nm)
for (l in layers[3]) {
  print(l)
  x <- paste0('SSc_HC_',l)
  #diff.genes <- rownames(qlf.res[[x]][qlf.res[[x]]$p_adj<0.05,])
  #pb_obj.deg <- pb_obj[diff.genes, , keep=FALSE] 
  pb_obj.deg <- pb_obj
  count_norm <- cpm(pb_obj.deg)
  count_norm <- as.data.frame(log(count_norm+1))
  count_norm['BAK1/BCL2',] <- count_norm['BAK1',]-count_norm['BCL2',]
  gene.results <- lapply(1:(nrow(pb_obj.deg$counts)+1), 
                         function(i){
                           print(c(i,rownames(count_norm[i,])))
                           gene.dat <- cbind.data.frame(gene=as.numeric(t(count_norm[i,])), layer, location,
                                                        pheno_flip,pheno_hrm, pheno_nm)
                           res <- c(rownames(count_norm[i,]))
                           names(res) <- c('gene')
                           for (p in 1:3) {
                             pheno <- phenos[[p]]
                             pheno_name <- paste0('pheno_',names(phenos)[p])
                             for (loc in c('Proximal','Distal')) {
                               suffix <- paste(names(phenos)[p],substr(loc,1,1),sep='.')
                               gene.stats <- tryCatch({
                                 gene.test <- orm(as.formula(paste0(pheno_name,'~gene')), gene.dat[gene.dat$layer==l & gene.dat$location==loc,])
                                 b_test <- coef(gene.test)['gene']
                                 p_test <- pchisq(b_test^2/vcov(gene.test)[2,2], 1,lower.tail=F)
                                 z_test <- sign(b_test)*qnorm(p_test/2)
                                 n_test <- sum(!is.na(pheno) & location==loc & layer==l)
                                 gene.stats <- c(b_test, p_test, z_test, n_test)
                               },
                               error=function(e){
                                 
                                 gene.stats <- c(rep(NA,4))
                               })
                               names(gene.stats) <- c(paste0('b_',suffix), paste0('p_',suffix),
                                                      paste0('z_',suffix), paste0('n_',suffix))
                               res <- c(res, gene.stats)
                             }
                           }
                           gene.stats <- data.frame(t(res), stringsAsFactors=F)
                           #print(res)
                           return(gene.stats)
                         })
  
  layer_results <- as.data.frame(do.call(rbind, gene.results))
  ord_reg.res[[l]] <- layer_results
}
View(ord_reg.res[[3]])

df <- ord_reg.res[[1]]
df[-1] <- lapply(df[-1], as.numeric)
df <- na.omit(df)

df$p_flip.acat <- ACAT(rbind(as.numeric(df[,'p_flip.P']),as.numeric(df[,'p_flip.D'])))
df$p_hrm.acat <- ACAT(rbind(as.numeric(df[,'p_hrm.P']),as.numeric(df[,'p_hrm.D'])))
df$p_nm.acat <- ACAT(rbind(as.numeric(df[,'p_nm.P']),as.numeric(df[,'p_nm.D'])))

pdf(file=paste0(output_dir, '/PseudoDEG_super_BAK1.pdf'), width=7, height=4)
par(mfrow=c(1,2), mar=c(3,2,4.5,1), oma=c(0.5,2.5,0,0))
for (test in c('FLIP','HRM')) {
  if (test=='FLIP') {pheno <- pheno_flip} else if (test=='HRM') {pheno <- pheno_hrm} else {pheno <- pheno_nm}
  data <- cbind.data.frame(gene=as.numeric(t(count_norm['BAK1',])), layer, pheno, condition, location)
  p <- boxplot(gene~pheno, data=data[data$layer==l,], ylab='', yaxt='n',
               xlab='', pch=16, cex.main=1.5,ylim=c(0,3))
  beeswarm(gene~pheno, data=data[(data$layer==l),], ylab=gene, add=T,
           pwbg=cond_cols[condition], pwpch=c(25,24)[location], cex=1.3)
  tick <- seq_along(p$names)
  axis(1, at = tick, labels = FALSE)
  axis(2, at=c(0,1,2,3), labels=c('0','1','2','3'))
  #text(tick, par("usr")[3] - 1, p$names, srt = 0, xpd = TRUE, cex=0.65)
  #points(gene~jitter(as.numeric(condition), amount=0.1), data=data[data$layer==l,],pch=16, cex=0.8)
  title(main=test, line=3, cex.main=1.2, font.main=2)
  if (test=='FLIP') {mtext('log(CPM+1)', 2, line=2.8,cex=1)}
  mtext(bquote(bolditalic(.('BAK1'))), side=3, padj=-1.5, las=1 ,cex=1, col='grey30')
}
dev.off()

#z_meta <- (n_flip.p*z_flip.p + n_hrm.p*z_hrm.p + n_flip.d*z_flip.d + n_hrm.d*z_hrm.d)/
#  sqrt(n_flip.p^2 + n_hrm.p^2 + n_flip.d^2 + n_hrm.d^2)
#p_stouffer <- 1-pnorm(abs(z_meta))
#gene.stats$p_stouffer <- p_stouffer

# Compute meta p-values for correlated variables
for (l in layers[3]) {
  df <- ord_reg.res[[l]]
  df[-1] <- lapply(df[-1], as.numeric)
  df <- na.omit(df)
  #df <- df
  p.flip.acat <- ACAT(rbind(as.numeric(df$p_flip.P), as.numeric(df$p_flip.D)))
  p.hrm.acat <- ACAT(rbind(as.numeric(df$p_hrm.P), as.numeric(df$p_hrm.D)))
  print(length(p_acat))
  n <- nrow(df)
  n_traits <- 4
  
  z_flip.p <- sign(df$b_flip.P)*qnorm(df$p_flip.P/2)
  z_hrm.p <- sign(df$b_hrm.P)*qnorm(df$p_hrm.P/2)
  z_flip.d <- sign(df$b_flip.D)*qnorm(df$p_flip.D/2)
  z_hrm.d <- sign(df$b_hrm.D)*qnorm(df$p_hrm.D/2)
  
  mu <- 2*n_traits
  sum_cov <- 0
  c_1 <- 3.9081; c_2 <- 0.0313; c_3 <- 0.1022; c_4 <- -0.1378; c_5 <- 0.0941
  
  fisher_s <- rep(0, n)
  
  for (i in 1:4) {
    fisher_s <- fisher_s -2*log(df[,c('p_flip.p','p_flip.d','p_hrm.p','p_hrm.d')][,i])
    
    for (i_o in c(1:4)[!1:4 %in% i]) {
      cor_fh <- cor(cbind(z_flip.p,z_flip.d,z_hrm.p,z_hrm.d)[,c(i,i_o)])[1,2]
      r <- cor_fh * (1 + ((1 - (cor_fh^2)) / (2 * (n - 3))))
      cov_fh <- (c_1*(r^2)) + (c_2*(r^4)) + (c_3*(r^6)) + (c_4*(r^8)) + (c_5*(r^10)) - 
        ((c_1 / n)*((1 - (r^2))^2))
      sum_cov <- sum_cov + cov_fh
    }
  }
  
  sigma_2_fh <- 4*n_traits*sum_cov  # 4 * len(traits) * sum_cov
  
  # get meta p-value
  fisher_s <- -2*log(df$p_flip.p) - 2*log(df$p_hrm.p) -2*log(df$p_flip.d) - 2*log(df$p_hrm.d)
  df$p_fish <- pchisq(fisher_s, df=2*n_traits, lower.tail=F)
  df$p_acat <- p_acat
  df$p_fh <- pgamma(q=fisher_s, shape=(mu^2)/sigma_2_fh, scale=sigma_2_fh/mu, lower.tail=F)
  
  df$p_adj <- p.adjust(df$p_fh, method='bonferroni')
  
  ord_reg.res[[l]] <- df[order(df$p_fh),]
}


View(ord_reg.res[['superficial']])

cor(ord_reg.res[['superficial']][,c('p_stouffer','p_fish','p_acat','p_fh')])
cor(ord_reg.res[['superficial']][,c('p_stouffer','p_fish','p_acat','p_fh')], method='spearman')

View(ord_reg.res[['suprabasal']])
View(ord_reg.res[['basal']])


genes <- c('BAK1','NRG4', 'CASP3')

genes <- c('CCL19','IL33', 'CARD18')

genes <- c('ADNP2','BAK1', 'CARD18')

count_norm <- cpm(pb_obj, log=F)
count_norm <- as.data.frame(log(count_norm+1))

l <- 'superficial'
#l <- 'basal'

#"#8156B3" "#e26b53" "#A8A39d"
genes <- 'BAK1'
#pdf(file=paste0(output_dir, '/PseudoDEG_apop_super.pdf'), width=7, height=8)
par(mfrow=c(length(genes),3), mar=c(3,2,4.5,1), oma=c(0.5,2.5,0,0))
for (gene in genes) {
  for (test in c('FLIP','HRM', 'NM')) {
    if (test=='FLIP') {pheno <- pheno_flip} else if (test=='HRM') {pheno <- pheno_hrm} else {pheno <- pheno_nm}
    data <- cbind.data.frame(gene=as.numeric(t(count_norm[gene,])), layer, pheno, condition, location)
    p <- boxplot(gene~pheno+location, data=data[data$layer==l,], ylab='', 
                 xlab='', pch=16, cex.main=1.5)
    beeswarm(gene~pheno+location, data=data[(data$layer==l),], ylab=gene, add=T,
             pwbg=cond_cols[condition], pwpch=c(25,24)[location], cex=1.2, method='compact')
    tick <- seq_along(p$names)
    axis(1, at = tick, labels = FALSE)
    #text(tick, par("usr")[3] - 1, p$names, srt = 0, xpd = TRUE, cex=0.65)
    #points(gene~jitter(as.numeric(condition), amount=0.1), data=data[data$layer==l,],pch=16, cex=0.8)
    title(main=bquote(bolditalic(.(gene))), line=3, cex.main=1.2, font.main=2)
    if (test=='FLIP') {mtext('log(CPM+1)', 2, line=2.8,cex=0.8)}
    mtext(test, side=3, padj=-1, las=1 ,cex=0.8, col='grey30')
  }
}
dev.off()

library(fgsea)
library(msigdbr)

c2_df <- msigdbr(species='Homo sapiens', category='C2')
c2_df <- c2_df[grepl('^CP', c2_df$gs_subcat),]
c7_df <- msigdbr(species='Homo sapiens', category='C7')

fgsea_c2_sets<- c2_df %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_c7_sets<- c7_df %>% split(x = .$gene_symbol, f = .$gs_name)


genes <- -log10(as.numeric(ord_reg.res[[l]][order(ord_reg.res[[l]]$p_nm.D),'p_nm.D']))

names(genes) <- ord_reg.res[[l]][order(ord_reg.res[[l]]$p_nm.D),'gene']

genes <- genes[!is.na(genes)]
genes <- genes[!grepl('^MT-', names(genes))]
genes <- genes[!grepl('^RP', names(genes))]
genes <- genes[!grepl('^MRP', names(genes))]

enrich_NM <- fgseaMultilevel(fgsea_c2_sets, stats=genes,
                             minSize=10, scoreType='pos', nPermSimple=10000)
View(enrich_NM)

pdf(file=paste0(output_dir, '/PseudoDEG_apopSuper_BAKratio.pdf'), width=7, height=4)
par(mfrow=c(1,2), mar=c(3,2,4.5,1), oma=c(0.5,2.5,0,0))
for (test in c('FLIP','HRM')) {
  if (test=='FLIP') {pheno <- pheno_flip} else {pheno <- pheno_hrm}
  data <- cbind.data.frame(gene=as.numeric(t(count_norm['BAK1',]-count_norm['BCL2',])), layer, pheno, condition, location)
  p <- boxplot(gene~pheno, data=data[data$layer==l,], ylab='', 
               xlab='', pch=16, cex.main=1.5,)
  beeswarm(gene~pheno, data=data[(data$layer==l),], ylab=gene, add=T,
           pwbg=cond_cols[condition], pwpch=c(25,24)[location], cex=1.2)
  tick <- seq_along(p$names)
  axis(1, at = tick, labels = FALSE)
  #text(tick, par("usr")[3] - 1, p$names, srt = 0, xpd = TRUE, cex=0.65)
  #points(gene~jitter(as.numeric(condition), amount=0.1), data=data[data$layer==l,],pch=16, cex=0.8)
  title(main=bquote(bolditalic(.('BAK1'))), line=3, cex.main=1.2, font.main=2)
  if (test=='FLIP') {mtext('log(CPM)', 2, line=2.8,cex=0.8)}
  mtext(test, side=3, padj=-1, las=1 ,cex=0.8, col='grey30')
}
dev.off()

levels(condition) <- c('HC','GERD','SSc')

pdf(file=paste0(output_dir, '/PseudoDEG_diffScores_cond.pdf'), width=6, height=6)
par(mfrow=c(2,2), mar=c(3,2,2.5,1), oma=c(0.5,2.5,0.25,0.25))
for (m in c('Mean', 'Median')) {
  for (loc in c('Proximal','Distal')) {
    data <- cbind.data.frame(med=pb_obj$samples$med, mean=pb_obj$samples$mean, layer, q=pb_obj$samples[,test], condition, location)
    #points(gene~jitter(as.numeric(condition), amount=0.1), data=data[data$layer==l,],pch=16, cex=0.8)
    #title(main=test, line=3, cex.main=1.2, font.main=2)
    if (m=='Median') {
      p <- boxplot(med~condition, data=data[data$layer==l & location==loc,], ylab='', 
                   xlab='', pch=16, cex.main=1.5, col=alpha(rev(cond_cols)),0.75, ylim=c(0.1, 0.7))
      beeswarm(med~condition, data=data[(data$layer==l) & location==loc,], ylab=gene, add=T, cex=1.2, pch=16)
    } else {
      p <- boxplot(mean~condition, data=data[data$layer==l & location==loc,], ylab='', 
                   xlab='', pch=16, cex.main=1.5, col=alpha(rev(cond_cols)),0.75, ylim=c(0.2, 0.6))
      beeswarm(mean~condition, data=data[(data$layer==l) & location==loc,], ylab=gene, add=T, cex=1.2, pch=16)
    }
    if (loc=='Proximal') {mtext(paste(m,'Diff. Score'), 2, line=2.8,cex=1)}
    mtext(loc, side=3, padj=-1, las=1 ,cex=0.8, col='grey30', font=2)
  }
}
dev.off()


### Test quant traits
quant_clin_traits <- c('median_irp','mean_DCI','pressure_60ml','egd_di_60ml',
                       's_basal_ee_egjp','egj_cont_index')



pdf(file=paste0(output_dir, '../epi_cells/quant_trait_hist.pdf'), width=7, height=7)
par(mfrow=c(3,2), mar=c(3,3,4.5,3), oma=c(0.5,2.5,0,0))
for (trait in quant_clin_traits) {
  n <- length(unique(pb_obj$samples[!is.na(pb_obj$samples[,trait]) & pb_obj$samples$layer==l, 'indv']))
  hist(pb_obj$samples[,trait], breaks=20, main=trait, ylab='count')
  mtext(paste('n=',n), side=3, adj=1,padj=-1, las=1 ,cex=0.8, col='grey30')
}
dev.off()

lm.res <- list()
for (loc in c('Proximal','Distal')) {
  lc <- ifelse(loc=='Proximal','p','d')
  for (l in layers[3]) {
    print(l)
    x <- paste('SSc','HC',l,loc, sep='_')
    diff.genes <- rownames(qlf.res[[paste('SSc','HC',l,loc,sep='_')]][qlf.res[[paste('SSc','HC',l,loc,sep='_')]]$p_adj<0.05,])
    #pb_obj.deg <- pb_obj[diff.genes, , keep=FALSE] 
    pb_obj.deg <- pb_obj
    count_norm <- cpm(pb_obj.deg)
    count_norm <- as.data.frame(log(count_norm+1))
    gene.results <- lapply(1:(nrow(pb_obj.deg$counts)), 
                           function(i){
                             gene.stats <- data.frame(row.names=rownames(count_norm)[i])
                             for (trait in quant_clin_traits) {
                               g <- rownames(count_norm)[i]
                               gene.dat <- cbind.data.frame(gene=as.numeric(t(count_norm[i,])), layer, location,
                                                            condition, qt=pb_obj.deg$samples[,trait])
                               gene.test <- lm(qt~gene, gene.dat[gene.dat$layer==l & gene.dat$location==loc
                                                                 & gene.dat$condition=='SSc',])
                               test_sum <- summary(gene.test)$coefficients
                               if(g=='LINC02494'){print(c(l,loc,i))
                                 print(test_sum)
                                 print(gene.dat[gene.dat$layer==l & gene.dat$location==loc
                                                & gene.dat$condition=='SSc',])}
                               corS <- cor(gene.dat[gene.dat$layer==l & gene.dat$location==loc,c('qt','gene')], 
                                           use='pairwise.complete.obs', method='spearman')[1,2]
                               b <- test_sum['gene',1]
                               e <- test_sum['gene',2]
                               p <- test_sum['gene',4]
                               gene.stats.trait <- data.frame(corS, b, e, p, stringsAsFactors = F)
                               res <- c(paste('corS',lc,trait,sep='.'),
                                        paste('b',lc,trait, sep='.'),
                                        paste('e',lc,trait, sep='.'),
                                        paste('p',lc,trait, sep='.'))
                               names(gene.stats.trait) <- res
                               rownames(gene.stats.trait) <- rownames(count_norm)[i]
                               gene.stats[rownames(count_norm)[i],res] <- gene.stats.trait
                             }
                             return(gene.stats)
                           })
    layer_results <- as.data.frame(do.call(rbind, gene.results))
    lm.res[[paste(l,loc,sep='_')]] <- layer_results
  }
}


# Compute meta p-values for correlated variables
for (l in layers) {
  df <- lm.res[[1]]
  for (trait in quant_clin_traits) {
    df[,paste0('p.acat.', trait)] <- ACAT(rbind(df[,paste('p.p',trait, sep='.')], df[,paste('p.d',trait, sep='.')]))
  }
  lm.res[[l]] <- df
}

library(fgsea)
library(msigdbr)

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

genes <- -log10(lm.res[[1]][order(lm.res[[1]]$p.p.mean_DCI),'p.p.mean_DCI'])

names(genes) <- rownames(lm.res[[1]][order(lm.res[[1]]$p.p.mean_DCI),])

genes <- genes[!is.na(genes)]
genes <- genes[!grepl('^MT-', names(genes))]
genes <- genes[!grepl('^RP', names(genes))]
genes <- genes[!grepl('^MRP', names(genes))]


enrich_DCI <- fgseaMultilevel(fgsea_c2_sets, stats=genes,
                              minSize=10, scoreType='pos', nPermSimple=10000)


fgenes <- c('SLC38A6','MT1X', 'BAK1')
genes <- c('MT1A','MT1X', 'MT1E')

count_norm <- cpm(pb_obj)
count_norm <- log(as.data.frame(count_norm+1))
l <- 'superficial'

genes <- c('CTBP1-DT','USPL1', 'IRX2')
test <- 'pressure_60ml'
#"#8156B3" "#e26b53" "#A8A39d"

#pdf(file=paste0(output_dir, '/PseudoDEG_apop_super.pdf'), width=7, height=8)
par(mfrow=c(length(genes),2), mar=c(3,2,4.5,1), oma=c(0.5,2.5,0,0))
for (gene in genes) {
  for (loc in c('Proximal','Distal')) {
    data <- cbind.data.frame(gene=as.numeric(t(count_norm[gene,])), layer, q=pb_obj$samples[,test], condition, location)
    plot(gene~q, data=data[data$layer==l & location==loc & condition=='SSc',], pch=16, cex.main=1.5, xlab=test)
    #points(gene~jitter(as.numeric(condition), amount=0.1), data=data[data$layer==l,],pch=16, cex=0.8)
    title(main=bquote(bolditalic(.(gene))), line=3, cex.main=1.2, font.main=2)
    if (loc=='Proximal') {mtext('log(CPM+1)', 2, line=2.8,cex=0.8)}
    mtext(loc, side=3, padj=-1, las=1 ,cex=0.8, col='grey30')
  }
}
dev.off()


par(mfrow=c(2,length(quant_clin_traits)), mar=c(3,2,4.5,1), oma=c(0.5,2.5,0,0))
for (loc in c('Proximal','Distal')) {
  for (test in quant_clin_traits) {
    data <- cbind.data.frame(med=pb_obj$samples$med, mean=pb_obj$samples$mean, layer, q=pb_obj$samples[,test], condition, location)
    plot(mean~q, data=data[data$layer==l & location==loc,], pch=16, cex.main=1.5, xlab=test)
    #points(gene~jitter(as.numeric(condition), amount=0.1), data=data[data$layer==l,],pch=16, cex=0.8)
    title(main=test, line=3, cex.main=1.2, font.main=2)
    if (test==quant_clin_traits[1]) {mtext('Med. Diff. Score', 2, line=2.8,cex=0.8)}
    mtext(loc, side=3, padj=-1, las=1 ,cex=0.8, col='grey30')
  }
}
dev.off()

## Plot differentation score vs quant traits
library(hrbrthemes)

plot_list <- list()
for (i in 1:length(quant_clin_traits)) {
  trait <- quant_clin_traits[[i]]
  trait_name <- names(quant_clin_traits[i])
  for (loc in c('Proximal','Distal')) {
    data <- pb_obj$samples[pb_obj$samples$location==loc,]
    plot_list[[paste(trait,loc,sep='_')]] <- ggplot(data, aes(x=.data[[trait]], y=med, group=.data[['condition']])) + 
      geom_point(aes(col=.data[['condition']])) + theme_ipsum() + scale_color_manual(name='Condition', values=cond_cols[c(2,3,1)]) +
      ggtitle(loc) +
      scale_fill_manual(name='Condition', values=alpha(cond_cols[c(2,3,1)],0.3)) + 
      theme(plot.title=element_text(size=12, face='italic', hjust=0.5),
            axis.title.x=element_text(size=11, face='bold'),
            axis.title.y=element_text(size=10.5, face='bold'))
  }
}
library(ggpubr)
for (i in seq(1, length(plot_list),2)) {
  trait <- quant_clin_traits[(i+1)/2]
  cairo_pdf(file=paste0(output_dir, paste('/QuantClinTrait','medDiff',trait,'regression.pdf', sep='_')), width=11, height=5)
  print(ggarrange(plot_list[[i]], plot_list[[i+1]], ncol=2, common.legend=T, legend='right'))
  dev.off()
}


for (var in quant_clin_traits) {
  pdf(file=paste0(output_dir, '/Quant_medDiffScore_Corr_',var,'.pdf'), width=10, height=6)
  par(mfrow=c(1,2))
  for (loc in rev(levels(location))) {
    plot_dat <- pb_obj$samples[pb_obj$samples$layer=='superficial' & pb_obj$samples$location==loc,]
    if (var=='mean_DCI') {plot_dat[,var] <- log(plot_dat[,var]+1)}
    model <- lm(as.formula(paste(var,'~ med')), data=plot_dat[plot_dat$condition=='HC',])
    mod_seq <- seq(0,1, length.out=100)
    preds <- predict(model, newdata = data.frame(med=mod_seq), interval = 'confidence')
    plot(as.formula(paste(var,'~ med')), data=plot_dat, type='n', ylab=var, xlim=c(0.1,0.6),
         xlab='Median Differentiation Score', pch=16, col=cond_cols[3], main=loc)
    polygon(c(rev(mod_seq), mod_seq), c(rev(preds[ ,3]), preds[ ,2]), col=alpha(cond_cols[3],0.2), border = NA)
    abline(model$coefficients[1],model$coefficients[2],lty=3, col=cond_cols[3])
    p_ssc <- summary(model)$coefficients[2,4]
    
    model <- lm(as.formula(paste(var,'~ med')), data=plot_dat[plot_dat$condition=='SSc',])
    mod_seq <- seq(0,1, length.out=100)
    preds <- predict(model, newdata = data.frame(med=mod_seq), interval = 'confidence')
    polygon(c(rev(mod_seq), mod_seq), c(rev(preds[ ,3]), preds[ ,2]), col=alpha(cond_cols[1],0.2), border = NA)
    abline(model$coefficients[1],model$coefficients[2],lty=3, col=cond_cols[1])
    p_hc <- summary(model)$coefficients[2,4]
    points(as.formula(paste(var,'~ med')), data=plot_dat[plot_dat$condition=='HC',], pch=18,col=cond_cols[3],main=loc,xlab='')
    points(as.formula(paste(var,'~ med')), data=plot_dat[plot_dat$condition=='SSc',], pch=15,col=cond_cols[1],main=loc,xlab='')
    text(0.6,max(plot_dat[,var],na.rm=T)*.99, paste('P =', round(p_ssc,3)), adj=1,col=cond_cols[1],cex=0.5)
    text(0.6,max(plot_dat[,var],na.rm=T)*.95,paste('P =', round(p_hc,3)), adj=1,col=cond_cols[3],cex=0.5)
  }
  dev.off()
}





### Analysis for Eso cells

load(paste0(output_dir, '/eso_pb.RData'))

pb_obj <- DGEList(sc_aggro@assays$RNA@counts)  # create DGEList

pb_obj <- calcNormFactors(pb_obj)

# add metadata to psuedobulk object
pb_obj$samples$sample <- rownames(pb_obj$samples)
pb_obj$samples$condition <- ifelse(grepl('^I', pb_obj$samples$sample), 'GERD',
                                   ifelse(grepl('^P', pb_obj$samples$sample), 'SSc', 'HC'))
pb_obj$samples$cellType <- word(rownames(pb_obj$samples),2,sep='\\_')
pb_obj$samples$location <- ifelse(grepl('D_', pb_obj$samples$sample), 'Distal','Proximal')
pb_obj$samples$group <- paste(pb_obj$samples$condition, pb_obj$samples$cellType, sep='_')

samples <- factor(sub("\\-.*$","", pb_obj$samples$sample))
sample_cols <- sample(alphabet(26),20)

cellType <- factor(pb_obj$samples$cellType)
condition <- factor(pb_obj$samples$condition, levels=c('SSc','GERD','HC'))
location <- factor(pb_obj$samples$location, levels=c('Distal','Proximal'))
group <- factor(pb_obj$samples$group, levels=unique(pb_obj$samples$group))

# filter non-expressed genes by subtype (logCPM >0 in >50%)
count_norm <- cpm(pb_obj, log=T)
min_cpm <- log(0.5)
ct.genes <- lapply(setNames(levels(cellType), levels(cellType)),
                   function(ct){
                     ct.df <- as.data.frame(count_norm[,cellType==ct])
                     ct.ssc.df <- as.data.frame(count_norm[,cellType==ct & condition=='SSc'])
                     ct.gerd.df <- as.data.frame(count_norm[,cellType==ct & condition=='GERD'])
                     ct.hc.df <- as.data.frame(count_norm[,cellType==ct & condition=='HC'])
                     ct.ssc.genes <- rownames(ct.ssc.df[rowSums(ct.ssc.df>min_cpm)>ncol(ct.ssc.df)*0.75,])
                     ct.gerd.genes <- rownames(ct.gerd.df[rowSums(ct.gerd.df>min_cpm)>ncol(ct.gerd.df)*0.75,])
                     ct.hc.genes <- rownames(ct.hc.df[rowSums(ct.hc.df>min_cpm)>ncol(ct.hc.df)*0.75,])
                     ct.outlier.genes <- rownames(ct.df[rowSums(ct.df>(rowMeans(ct.df)+3*rowSds(as.matrix(ct.df))))==1,])
                     ct.outlier.genes2 <- rownames(ct.df[rowSums(ct.df>(rowMeans(ct.df)+3*rowSds(as.matrix(ct.df))))==2,])
                     for (gene in ct.outlier.genes2) {
                       top_two <- sub("\\-.*$","", names(ct.df[gene,order(as.numeric(ct.df[gene,]),decreasing=T)][1:2]))
                       if (top_two[1]==top_two[2]) {
                         ct.outlier.genes <- c(ct.outlier.genes, gene)}}
                     ct.genes <- unique(c(ct.ssc.genes,ct.gerd.genes,ct.hc.genes))
                     ct.genes <- ct.genes[!ct.genes %in% ct.outlier.genes]
                     return(ct.genes)
                   })

pb_obj <- pb_obj[unique(unlist(ct.genes)), , keep=FALSE] 

clin_dat <- read.delim(paste0(output_dir,'../clinical_data.txt'), na.strings=c('NA','NULL','null'))

pb_obj$samples$indv <- sapply(pb_obj$samples$sample, function(x) {strsplit(x,'-')[[1]][1]})

pb_obj$samples <- merge(pb_obj$samples, clin_dat[,c('Participant_ID','median_irp','mean_DCI',
                                                    'pressure_60ml','egd_di_60ml','s_basal_ee_egjp','egj_cont_index',
                                                    'cc_v4_class_hrm','flip_contract_pattern_v2')],
                        by.x='indv', by.y='Participant_ID')

pb_obj$samples$cc_v4_class_hrm <- toupper(pb_obj$samples$cc_v4_class_hrm)
library(plyr)
hrm_vals <- c('Normal','IEM','Absent','EGJOO/Alachasia','EGJOO/Alachasia')
flip_vals <- c('Normal','Diminished', 'Impaired','Absent')
pb_obj$samples$cc_v4_class_hrm <- mapvalues(pb_obj$samples$cc_v4_class_hrm, c('NORMAL ESOPHAGEAL MOTILITY',
                                                                              'INEFFECTIVE ESOPHAGEAL MOTILITY (IEM)',
                                                                              'ABSENT CONTRACTILITY','EGJOO', 'TYPE II ACHALASIA'), hrm_vals)
pb_obj$samples$flip_contract_pattern_v2 <- sub(" CONTRACTILE RESPONSE","", 
                                               pb_obj$samples$flip_contract_pattern_v2)
pb_obj$samples$flip_contract_pattern_v2 <- mapvalues(pb_obj$samples$flip_contract_pattern_v2, 
                                                     c('NORMAL','BORDERLINE/DIMINISHED', 'IMPAIRED/DISORDERED',
                                                       'ABSENT'), flip_vals)

#pb_obj <- pb_obj[, pb_obj$samples$condition=='SSc', keep=FALSE] 

prop_dat <- read.delim(paste0(output_dir,'../eso_props.tsv'), na.strings=c('NA','NULL','null'))
prop_dat$group <- paste(prop_dat$Sample, prop_dat$Cells, sep='_')
sample_dat <- pb_obj$samples
sample_dat <- merge(sample_dat, prop_dat[,c('Sample','Cells','Count',
                                            'Total','Proportion','group')],
                    by.x='sample', by.y='group')

sample_dat$flip_contract_pattern_v2 <- factor(sample_dat$flip_contract_pattern_v2, ordered=T, 
                                              exclude=NA,levels=flip_vals)
pheno_hrm <- factor(pb_obj$samples$cc_v4_class_hrm, ordered=T, levels=hrm_vals[1:4])


condition <- factor(sample_dat$condition, levels=c('SSc','GERD','HC'))

ct <- 'En'
par(mfrow=c(2,2))
for (loc in c('Proximal','Distal')) {
  dat <- sample_dat[!is.na(sample_dat$flip_contract_pattern_v2) & sample_dat$cellType==ct & sample_dat$location==loc,]
  beeswarm(Proportion ~ flip_contract_pattern_v2, data=dat, pch=21, cex=1.5, method='swarm',
           pwbg=cond_cols[factor(dat$condition, levels=c('SSc','GERD','HC'))], main=paste(ct, '-',loc))
  dat <- sample_dat[!is.na(sample_dat$cc_v4_class_hrm) & sample_dat$cellType==ct & sample_dat$location==loc,]
  beeswarm(Proportion ~ cc_v4_class_hrm, data=dat, pch=21, cex=1.5, method='swarm',
           pwbg=cond_cols[factor(dat$condition, levels=c('SSc','GERD','HC'))], main=paste(ct, '-',loc))
  
}

cols=c('plum3','darkseagreen','cadetblue3','lightpink2','grey')
library(ggpubr)
pdf(file=paste0(output_dir, '/PseudoDEG_cellType_x_flip_ssc.pdf'), width=9, height=4)
print(ggbarplot(sample_dat[!is.na(sample_dat$flip_contract_pattern_v2) & sample_dat$condition=='SSc' &
                             sample_dat$cellType %in% c('En','Ep','L','My'),], 
                x='flip_contract_pattern_v2', y='Proportion', add=c('mean_sd','jitter'), fill='cellType', 
                position=position_dodge(0.8), add.params=list(color='cellType'), palette=cols) + 
        scale_color_manual(values=rep('black',4))+ 
        theme(axis.title.x=element_text(size=12, face='bold'), axis.text.y=element_text(size=8),
              plot.title=element_text(face='bold',hjust=0.5)))
dev.off()




load(paste0(output_dir,'epi_aggro.RData'))

pb_obj <- DGEList(epi_aggro@assays$RNA@counts)  # create DGEList
pb_obj <- calcNormFactors(pb_obj)

# add metadata to psuedobulk object
pb_obj$samples$sample <- rownames(pb_obj$samples)
pb_obj$samples$condition <- ifelse(grepl('^I', pb_obj$samples$sample), 'GERD',
                                   ifelse(grepl('^P', pb_obj$samples$sample), 'SSc', 'HC'))
pb_obj$samples$layer <- word(rownames(pb_obj$samples),2,sep='\\_')
pb_obj$samples$location <- ifelse(grepl('D_', pb_obj$samples$sample), 'Distal','Proximal')
pb_obj$samples$group <- paste(pb_obj$samples$condition, pb_obj$samples$layer, sep='_')

samples <- factor(sub("\\-.*$","", pb_obj$samples$sample))

layer <- factor(pb_obj$samples$layer, levels=c('basal','suprabasal','superficial'))
condition <- factor(pb_obj$samples$condition, levels=c('SSc','GERD','HC'))
location <- factor(pb_obj$samples$location, levels=c('Distal','Proximal'))
group <- factor(pb_obj$samples$group, levels=unique(pb_obj$samples$group))
pb_obj <- pb_obj[unique(unlist(ct.genes)), , keep=FALSE] 

pb_obj <- pb_obj[, pb_obj$samples$condition=='SSc', keep=FALSE] 

clin_dat <- read.delim(paste0(output_dir,'../clinical_data.txt'), na.strings=c('NA','NULL','null'))

pb_obj$samples$indv <- sapply(pb_obj$samples$sample, function(x) {strsplit(x,'-')[[1]][1]})

pb_obj$samples <- merge(pb_obj$samples, clin_dat[,c('Participant_ID','median_irp','mean_DCI',
                                                    'pressure_60ml','egd_di_60ml','s_basal_ee_egjp','egj_cont_index',
                                                    'cc_v4_class_hrm','flip_contract_pattern_v2')],
                        by.x='indv', by.y='Participant_ID')
pb_obj$samples$cc_v4_class_hrm <- toupper(pb_obj$samples$cc_v4_class_hrm)
#library(plyr)
hrm_vals <- c('Normal','IEM','Absent','EGJOO/Alachasia','EGJOO/Alachasia')
flip_vals <- c('Normal','Diminished', 'Impaired','Absent')
pb_obj$samples$cc_v4_class_hrm <- mapvalues(pb_obj$samples$cc_v4_class_hrm, c('NORMAL ESOPHAGEAL MOTILITY',
                                                                              'INEFFECTIVE ESOPHAGEAL MOTILITY (IEM)',
                                                                              'ABSENT CONTRACTILITY','EGJOO', 'TYPE II ACHALASIA'), hrm_vals)
pb_obj$samples$flip_contract_pattern_v2 <- sub(" CONTRACTILE RESPONSE","", 
                                               pb_obj$samples$flip_contract_pattern_v2)
pb_obj$samples$flip_contract_pattern_v2 <- mapvalues(pb_obj$samples$flip_contract_pattern_v2, 
                                                     c('NORMAL','BORDERLINE/DIMINISHED', 'IMPAIRED/DISORDERED',
                                                       'ABSENT'), flip_vals)

prop_dat <- read.delim(paste0(output_dir,'../epi_cells/epi_layer_props.tsv'), na.strings=c('NA','NULL','null'))
prop_dat$group <- paste(prop_dat$Sample, tolower(prop_dat$Layer), sep='_')
sample_dat <- pb_obj$samples
sample_dat <- merge(sample_dat, prop_dat[,c('Sample','Layer','Proportion','group')],
                    by.x='sample', by.y='group')
sample_dat$Layer <- factor(sample_dat$Layer, levels=c('Basal','Suprabasal','Superficial'))

sample_dat$flip_contract_pattern_v2 <- factor(sample_dat$flip_contract_pattern_v2, ordered=T, 
                                              exclude=NA,levels=flip_vals)
sample_dat$cc_v4_class_hrm <- factor(sample_dat$cc_v4_class_hrm, ordered=T, levels=hrm_vals[1:4])


cols=c('plum3','darkseagreen','cadetblue3','lightpink2','grey')
layer_cols_3 <- c('#4ca5b1','#e9b85d','#c72f4c')

#library(ggpubr)
pdf(file=paste0(output_dir, '../PseudoDEG_layer_x_flip_ssc.pdf'), width=8, height=4)
print(ggbarplot(sample_dat[!is.na(sample_dat$flip_contract_pattern_v2),], 
                x='flip_contract_pattern_v2', y='Proportion', add=c('mean_sd'), fill='Layer', 
                position=position_dodge(0.8), palette=layer_cols_3) + 
        geom_jitter(aes(flip_contract_pattern_v2, Proportion, color=condition,
                        group=sample_dat[!is.na(sample_dat$flip_contract_pattern_v2),'Layer']),
                    fill=cond_cols[c(3,1)][as.factor(sample_dat[!is.na(sample_dat$flip_contract_pattern_v2),'condition'])], 
                    position=position_dodge(0.8), pch=21, size=2) +
        labs(fill='Layer', color='Condition') + scale_color_manual(values=c(1,1)) + 
        guides(color=guide_legend(override.aes=list(fill=cond_cols[c(3,1)]), 
                                  title.theme=element_text(face='bold', size=11)),
               fill=guide_legend(title.theme=element_text(face='bold',size=11),
                                 label.theme=element_text(size=9))) +
        theme(axis.title.x=element_text(size=12, face='bold'), axis.text.y=element_text(size=8),
              plot.title=element_text(face='bold',hjust=0.5)))
dev.off()

pdf(file=paste0(output_dir, '../PseudoDEG_layer_x_hrm_ssc.pdf'), width=8, height=4)
print(ggbarplot(sample_dat[!is.na(sample_dat$cc_v4_class_hrm),], 
                x='cc_v4_class_hrm', y='Proportion', add=c('mean_sd'), fill='Layer', 
                position=position_dodge(0.8), palette=layer_cols_3) + 
        geom_jitter(aes(cc_v4_class_hrm, Proportion, color=factor(condition, levels=c('HC','GERD','SSc')),
                        group=sample_dat[!is.na(sample_dat$cc_v4_class_hrm),'Layer']),
                    fill=cond_cols[c(3,2,1)][factor(sample_dat[!is.na(sample_dat$cc_v4_class_hrm),'condition'], levels=c('HC','GERD','SSc'))], 
                    position=position_dodge(0.8), pch=21,size=1.5) +
        labs(fill='Layer', color='Condition') + scale_color_manual(values=c(1,1,1)) + 
        guides(color=guide_legend(override.aes=list(fill=cond_cols[c(3,2,1)]), 
                                  title.theme=element_text(face='bold', size=11)),
               fill=guide_legend(title.theme=element_text(face='bold',size=11),
                                 label.theme=element_text(size=9))) +
        theme(axis.title.x=element_text(size=12, face='bold'), axis.text.y=element_text(size=8),
              plot.title=element_text(face='bold',hjust=0.5)))
dev.off()


prop_dat <- read.delim(paste0(output_dir,'../epi_cells/epi_layerRep_props.tsv'), na.strings=c('NA','NULL','null'))
prop_dat$group <- paste(prop_dat$Sample, tolower(prop_dat$LayerRep), sep='_')
sample_dat <- pb_obj$samples
sample_dat$sample <- word(sample_dat$sample,1,sep='\\_')

sample_dat <- merge(prop_dat[,c('Sample','LayerRep','Proportion','group')],
                    sample_dat[sample_dat$layer=='basal',c('sample','condition','cc_v4_class_hrm','flip_contract_pattern_v2')],
                    by.y='sample', by.x='Sample')
sample_dat$LayerRep <- factor(sample_dat$LayerRep, levels=c('Basal','Replicating_basal',
                                                            'Replicating_suprabasal','Suprabasal','Superficial'))

sample_dat$flip_contract_pattern_v2 <- factor(sample_dat$flip_contract_pattern_v2, ordered=T, 
                                              exclude=NA,levels=flip_vals)
sample_dat$cc_v4_class_hrm <- factor(sample_dat$cc_v4_class_hrm, ordered=T, levels=hrm_vals[1:4])


layer_cols_5 <- brewer.pal(n=11,'Spectral')[c(10,9,7,4,2)]

#library(ggpubr)
pdf(file=paste0(output_dir, '../PseudoDEG_layerRep_x_flip_ssc.pdf'), width=9.5, height=4)
print(ggbarplot(sample_dat[!is.na(sample_dat$flip_contract_pattern_v2),], 
                x='flip_contract_pattern_v2', y='Proportion', add=c('mean_sd'), fill='LayerRep', 
                position=position_dodge(0.8), palette=layer_cols_5) + 
        geom_jitter(aes(flip_contract_pattern_v2, Proportion, color=factor(condition, levels=c('HC','GERD','SSc')),
                        group=sample_dat[!is.na(sample_dat$flip_contract_pattern_v2),'LayerRep']),
                    fill=cond_cols[c(3,2,1)][factor(sample_dat[!is.na(sample_dat$flip_contract_pattern_v2),'condition'], levels=c('HC','GERD','SSc'))], 
                    position=position_dodge(0.8), pch=21,size=1.5) +
        labs(fill='Layer', color='Condition') + scale_color_manual(values=c(1,1)) + 
        guides(color=guide_legend(override.aes=list(fill=cond_cols[c(3,1)]), 
                                  title.theme=element_text(face='bold', size=11)),
               fill=guide_legend(title.theme=element_text(face='bold',size=11),
                                 label.theme=element_text(size=9))) +
        theme(axis.title.x=element_text(size=12, face='bold'), axis.text.y=element_text(size=8),
              plot.title=element_text(face='bold',hjust=0.5), legend.position='right'))
dev.off()

pdf(file=paste0(output_dir, '../PseudoDEG_layerRep_x_hrm_ssc.pdf'), width=9.5, height=4)
print(ggbarplot(sample_dat[!is.na(sample_dat$cc_v4_class_hrm),], 
                x='cc_v4_class_hrm', y='Proportion', add=c('mean_sd'), fill='LayerRep', 
                position=position_dodge(0.8), palette=layer_cols_5) + 
        geom_jitter(aes(cc_v4_class_hrm, Proportion, color=factor(condition, levels=c('HC','GERD','SSc')),
                        group=sample_dat[!is.na(sample_dat$cc_v4_class_hrm),'LayerRep']),
                    fill=cond_cols[c(3,2,1)][factor(sample_dat[!is.na(sample_dat$cc_v4_class_hrm),'condition'], levels=c('HC','GERD','SSc'))], 
                    position=position_dodge(0.8), pch=21,size=1.5) +
        labs(fill='Layer', color='Condition') + scale_color_manual(values=c(1,1,1)) + 
        guides(color=guide_legend(override.aes=list(fill=cond_cols[c(3,2,1)]), 
                                  title.theme=element_text(face='bold', size=11)),
               fill=guide_legend(title.theme=element_text(face='bold',size=11),
                                 label.theme=element_text(size=9))) +
        theme(axis.title.x=element_text(size=12, face='bold'), axis.text.y=element_text(size=8),
              plot.title=element_text(face='bold',hjust=0.5), legend.position='right'))
dev.off()


prop_dat <- read.delim(paste0(output_dir,'../epi_cells/epi_cluster_props.tsv'), na.strings=c('NA','NULL','null'))
prop_dat$group <- paste(prop_dat$Sample, tolower(prop_dat$Cluster), sep='_')
sample_dat <- pb_obj$samples
sample_dat$sample <- word(sample_dat$sample,1,sep='\\_')

sample_dat <- merge(prop_dat[,c('Sample','Location','Cluster','Proportion','group')],
                    sample_dat[sample_dat$layer=='basal',c('sample','condition','cc_v4_class_hrm','flip_contract_pattern_v2')],
                    by.y='sample', by.x='Sample')
sample_dat$Cluster <- factor(sample_dat$Cluster, levels=c(1:8))

sample_dat$flip_contract_pattern_v2 <- factor(sample_dat$flip_contract_pattern_v2, ordered=T, 
                                              exclude=NA,levels=flip_vals)
sample_dat$cc_v4_class_hrm <- factor(sample_dat$cc_v4_class_hrm, ordered=T, levels=hrm_vals[1:4])

#library(scales)
cluster_cols <- hue_pal()(8)

#library(ggpubr)
pdf(file=paste0(output_dir, '../epi_cells/PseudoDEG_epiCluster_x_flip_ssc.pdf'), width=8, height=4)
print(ggbarplot(sample_dat[!is.na(sample_dat$flip_contract_pattern_v2),], 
                x='flip_contract_pattern_v2', y='Proportion', add=c('mean_sd'), fill='Cluster', 
                position=position_dodge(0.8), palette=cluster_cols) + 
        geom_jitter(aes(flip_contract_pattern_v2, Proportion, color=condition,
                        group=sample_dat[!is.na(sample_dat$flip_contract_pattern_v2),'Cluster']),
                    fill=cond_cols[c(3,1)][as.factor(sample_dat[!is.na(sample_dat$flip_contract_pattern_v2),'condition'])], 
                    position=position_dodge(0.8), pch=21, size=1.5) +
        labs(fill='Cluster', color='Condition') + scale_color_manual(values=c(1,1)) + 
        guides(color=guide_legend(override.aes=list(fill=cond_cols[c(3,1)]), 
                                  title.theme=element_text(face='bold', size=11)),
               fill=guide_legend(title.theme=element_text(face='bold',size=11),
                                 label.theme=element_text(size=9))) +
        theme(axis.title.x=element_text(size=12, face='bold'), axis.text.y=element_text(size=8),
              plot.title=element_text(face='bold',hjust=0.5)))
dev.off()

pdf(file=paste0(output_dir, '../epi_cells/PseudoDEG_epiCluster_x_hrm_ssc.pdf'), width=8, height=4)
print(ggbarplot(sample_dat[!is.na(sample_dat$cc_v4_class_hrm),], 
                x='cc_v4_class_hrm', y='Proportion', add=c('mean_sd'), fill='Cluster', 
                position=position_dodge(0.8), palette=cluster_cols) + 
        geom_jitter(aes(cc_v4_class_hrm, Proportion, color=factor(condition, levels=c('HC','GERD','SSc')),
                        group=sample_dat[!is.na(sample_dat$cc_v4_class_hrm),'Cluster']),
                    fill=cond_cols[c(3,2,1)][factor(sample_dat[!is.na(sample_dat$cc_v4_class_hrm),'condition'], levels=c('HC','GERD','SSc'))], 
                    position=position_dodge(0.8), pch=21,size=1.5) +
        labs(fill='Cluster', color='Condition') + scale_color_manual(values=c(1,1,1)) + 
        guides(color=guide_legend(override.aes=list(fill=cond_cols[c(3,2,1)]), 
                                  title.theme=element_text(face='bold', size=11)),
               fill=guide_legend(title.theme=element_text(face='bold',size=11),
                                 label.theme=element_text(size=9))) +
        theme(axis.title.x=element_text(size=12, face='bold'), axis.text.y=element_text(size=8),
              plot.title=element_text(face='bold',hjust=0.5)))
dev.off()


### Test quant trait associations with cell type proportions


library(hrbrthemes)
library(gridExtra)

get_prop_stats <- function(prop_file, pb_data, traits, grouping, colors) {
  
  #prop_file <- paste0(output_dir,'../epi_cells/epi_cluster_props.tsv') 
  #pb_data <- pb_obj$samples
  #traits <- quant_clin_traits
  #grouping <- 'Cluster'
  #colors <- hue_pal()(8)
  
  prop_dat <- read.delim(prop_file, na.strings=c('NA','NULL','null'))
  prop_dat$group <- paste(prop_dat$Sample, tolower(prop_dat[,grouping]), sep='_')
  pb_data$sample <- ifelse(pb_data$location=='Proximal',
                           paste0(pb_data$indv, '-','P'),
                           paste0(pb_data$indv,'-','D'))
  sample_dat <- merge(prop_dat, pb_data[pb_data$layer=='basal', 
                                        c('sample',quant_clin_traits)], by.x='Sample',by.y='sample')
  data <- sample_dat
  types <- as.factor(unique(data[,grouping]))
  names(colors) <- types
  prop.stats <- data.frame()
  for (type in types) {
    for (trait in traits) {
      for (loc in c('Proximal','Distal')) {
        lc <- ifelse(loc=='Proximal','p','d')
        #print(c(type, trait, loc))
        #print(head(data[data[,grouping]==type & data$Location==loc,c('Sample',grouping,'Location',trait,'Proportion')]))
        prop.test <- lm(data[data[,grouping]==type & data$Location==loc,trait] ~ 
                          data[data[,grouping]==type & data$Location==loc,'Proportion'])
        prop.test.spearman <- cor.test(data[data[,grouping]==type & data$Location==loc,trait], 
                                       data[data[,grouping]==type & data$Location==loc,'Proportion'],
                                       method='spearman', exact=F)
        test_sum <- summary(prop.test)$coefficients
        b <- test_sum[2,1]
        e <- test_sum[2,2]
        p <- test_sum[2,4]
        p_s <- prop.test.spearman$p.value
        cor <- cor(data[data[,grouping]==type & data$Location==loc,c(trait,'Proportion')])[1,2]
        cor_s <- cor(data[data[,grouping]==type & data$Location==loc,c(trait,'Proportion')],
                     method='spearman')[1,2]
        prop.stats.trait <- data.frame(cor, p, cor_s, p_s, stringsAsFactors = F)
        res <- c(paste('cor', lc, trait, sep='.'), paste('p',lc,trait, sep='.'),
                 paste('cor_s', lc, trait, sep='.'), paste('p_s',lc,trait, sep='.'))
        names(prop.stats.trait) <- res
        rownames(prop.stats.trait) <- rownames(count_norm)[i]
        prop.stats[type,res] <- prop.stats.trait
      }
    }
  }
  
  plot_list <- list()
  for (i in 1:length(quant_clin_traits)) {
    trait <- quant_clin_traits[[i]]
    trait_name <- names(quant_clin_traits[i])
    for (loc in c('Proximal','Distal')) {
      data <- sample_dat[sample_dat$Location==loc,]
      data[,grouping] <- as.factor(data[,grouping])
      plot_list[[paste(trait,loc,sep='_')]] <- ggplot(data, aes(x=.data[[trait]], y=Proportion, group=.data[[grouping]])) + 
        geom_point(aes(col=.data[[grouping]])) + theme_ipsum() + scale_color_manual(name='Epithelial cells', values=colors) +
        geom_smooth(aes(col=.data[[grouping]], fill=.data[[grouping]]), method=lm,se=T) + ggtitle(loc) +
        scale_fill_manual(name='Epithelial cells', values=alpha(colors,0.3)) + 
        theme(plot.title=element_text(size=12, face='italic', hjust=0.5),
              axis.title.x=element_text(size=11, face='bold'),
              axis.title.y=element_text(size=10.5, face='bold')) + labs(x=trait_name)
    }
  }
  
  for (i in seq(1, length(plot_list),2)) {
    trait <- quant_clin_traits[(i+1)/2]
    cairo_pdf(file=paste0(output_dir, paste('/QuantClinTrait',grouping,trait,'regression.pdf', sep='_')), width=11, height=5)
    print(ggarrange(plot_list[[i]], plot_list[[i+1]], ncol=2, common.legend=T, legend='right'))
    dev.off()
  }
  
  return(prop.stats)
}


names(quant_clin_traits) <- c('Median IRP','Mean DCI','Pressure (60ml)', 'EGD DI (60ml)', 
                              'End-expiratory EGJ pressure', 'EGJ contractile index')

prop.stats.layerRep <- get_prop_stats(paste0(output_dir,'../epi_cells/epi_layerRep_props.tsv'), 
                                      pb_obj$samples, quant_clin_traits,'LayerRep', layer_cols_5)

prop.stats.layers <- get_prop_stats(paste0(output_dir,'../epi_cells/epi_layer_props.tsv'),
                                    pb_obj$samples, quant_clin_traits,'Layer', layer_cols_3)

prop.stats.clusters <- get_prop_stats(paste0(output_dir,'../epi_cells/epi_cluster_props.tsv'),
                                      pb_obj$samples, quant_clin_traits,'Cluster', hue_pal()(8))


grouping <- 'Cluster'
prop_dat <- read.delim(paste0(output_dir,'../epi_cells/epi_cluster_props.tsv'), na.strings=c('NA','NULL','null'))
prop_dat$group <- paste(prop_dat$Sample, tolower(prop_dat[,grouping]), sep='_')
pb_data$sample <- ifelse(pb_data$location=='Proximal',
                         paste0(pb_data$indv, '-','P'),
                         paste0(pb_data$indv,'-','D'))
sample_dat <- merge(prop_dat, pb_data[pb_data$layer=='basal', 
                                      c('sample',quant_clin_traits)], by.x='Sample',by.y='sample')

for (i in unique(sample_dat[,grouping])) {
  sample_dat[nrow(sample_dat)+1,] <- sample_dat[sample_dat$Sample=='P0080-D' & sample_dat[,grouping]==i,]
  sample_dat[nrow(sample_dat),c('Sample','Location','nCells','totCells','Proportion')] <- c('P0080-P','Proximal',0,0,0)
}
#sample_dat <- head(sample_dat,152)
sample_dat$nCells <- as.numeric(sample_dat$nCells)
sample_dat[,grouping] <- factor(sample_dat[,grouping], levels=rev(levels(sample_dat[,grouping])))

pdf(file=paste0(output_dir, paste('../epi_cells/Epi_sampleCellCounts',grouping,'.pdf', sep='_')), width=8, height=5)
print(ggarrange(ggplot(sample_dat[sample_dat$Location=='Proximal',], aes(x=indv, y=nCells, fill=Cluster)) + 
                  geom_bar(stat='identity',position='stack') + theme_minimal() + ggtitle('Proximal') + 
                  scale_fill_manual(values=rev(cluster_cols)) + scale_y_continuous(limits=c(0,9000)),
                ggplot(sample_dat[sample_dat$Location=='Distal',], aes(x=indv, y=nCells, fill=Cluster)) + 
                  geom_bar(stat='identity',position='stack') + theme_minimal() + ggtitle('Distal') + 
                  scale_fill_manual(values=rev(cluster_cols)) + scale_y_continuous(limits=c(0,9000))
                , ncol=1, common.legend=T, legend='right'))
dev.off()

pdf(file=paste0(output_dir, paste('../epi_cells/Epi_sampleProportions',grouping,'.pdf', sep='_')), width=8, height=5)
print(ggarrange(ggplot(sample_dat[sample_dat$Location=='Proximal',], aes(x=indv, y=nCells, fill=Cluster)) + 
                  geom_bar(stat='identity',position='fill') + theme_minimal() + ggtitle('Proximal') + 
                  scale_fill_manual(values=rev(cluster_cols)),
                ggplot(sample_dat[sample_dat$Location=='Distal',], aes(x=indv, y=nCells, fill=Cluster)) + 
                  geom_bar(stat='identity',position='fill') + theme_minimal() + ggtitle('Distal') + 
                  scale_fill_manual(values=rev(cluster_cols)), ncol=1, common.legend=T, legend='right'))
dev.off()


for (i in 1:8) {
  dat <- merge(sample_dat[sample_dat$Cluster==i & sample_dat$Location=='Proximal',],
               sample_dat[sample_dat$Cluster==i & sample_dat$Location=='Distal',], 
               by='indv')[,c('Proportion.x','Proportion.y')]
  dat <- tail(dat,9)
  #print(dat)
  print(c(i,round(cor(cbind(as.numeric(dat$Proportion.x),as.numeric(dat$Proportion.y)))[1,2],2)))}

#library(scales)
cluster_cols <- hue_pal()(8)


# "#3288BD" "#66C2A5" "#E6F598" "#FDAE61" "#D53E4F"


grouping <- 'Layer'
cluster_prop_dat <- prop_dat <- read.delim(paste0(output_dir,'../epi_cells/epi_cluster_props.tsv'), na.strings=c('NA','NULL','null'))

prop_dat <- read.delim(paste0(output_dir,'../epi_cells/epi_layer_props.tsv'), na.strings=c('NA','NULL','null'))
prop_dat <- merge(prop_dat, cluster_prop_dat[cluster_prop_dat$Cluster==1,c('Sample','totCells')], by='Sample')
pb_data$sample <- ifelse(pb_data$location=='Proximal',
                         paste0(pb_data$indv, '-','P'),
                         paste0(pb_data$indv,'-','D'))
sample_dat <- merge(prop_dat, pb_data[pb_data$layer=='basal', 
                                      c('sample','indv',quant_clin_traits)], by.x='Sample',by.y='sample')
sample_dat$nCells <- sample_dat$totCells*sample_dat$Proportion

for (i in unique(sample_dat[,grouping])) {
  sample_dat[nrow(sample_dat)+1,] <- sample_dat[sample_dat$Sample=='P0080-D' & sample_dat[,grouping]==i,]
  sample_dat[nrow(sample_dat),c('Sample','Location','nCells','totCells','Proportion')] <- c('P0080-P','Proximal',0,0,0)
}
#sample_dat <- head(sample_dat,152)
sample_dat$nCells <- as.numeric(sample_dat$nCells)
sample_dat[,grouping] <- factor(sample_dat[,grouping], levels=c('Superficial','Suprabasal','Basal'))

pdf(file=paste0(output_dir, paste('../epi_cells/Epi_sampleCellCounts',grouping,'.pdf', sep='_')), width=8, height=5)
print(ggarrange(ggplot(sample_dat[sample_dat$Location=='Proximal',], aes(x=indv, y=nCells, fill=Layer)) + 
                  geom_bar(stat='identity',position='stack') + theme_minimal() + ggtitle('Proximal') + 
                  scale_fill_manual(values=rev(layer_cols_3)) + scale_y_continuous(limits=c(0,9000)),
                ggplot(sample_dat[sample_dat$Location=='Distal',], aes(x=indv, y=nCells, fill=Layer)) + 
                  geom_bar(stat='identity',position='stack') + theme_minimal() + ggtitle('Distal') + 
                  scale_fill_manual(values=rev(layer_cols_3)) + scale_y_continuous(limits=c(0,9000))
                , ncol=1, common.legend=T, legend='right'))
dev.off()

pdf(file=paste0(output_dir, paste('../epi_cells/Epi_sampleProportions',grouping,'.pdf', sep='_')), width=8, height=5)
print(ggarrange(ggplot(sample_dat[sample_dat$Location=='Proximal',], aes(x=indv, y=nCells, fill=Layer)) + 
                  geom_bar(stat='identity',position='fill') + theme_minimal() + ggtitle('Proximal') + 
                  scale_fill_manual(values=rev(layer_cols_3)),
                ggplot(sample_dat[sample_dat$Location=='Distal',], aes(x=indv, y=nCells, fill=Layer)) + 
                  geom_bar(stat='identity',position='fill') + theme_minimal() + ggtitle('Distal') + 
                  scale_fill_manual(values=rev(layer_cols_3)), ncol=1, common.legend=T, legend='right'))
dev.off()





grouping <- 'LayerRep'
cluster_prop_dat <- read.delim(paste0(output_dir,'../epi_cells/epi_cluster_props.tsv'), na.strings=c('NA','NULL','null'))

prop_dat <- read.delim(paste0(output_dir,'../epi_cells/epi_layerRep_props.tsv'), na.strings=c('NA','NULL','null'))
prop_dat <- merge(prop_dat, cluster_prop_dat[cluster_prop_dat$Cluster==1,c('Sample','totCells')], by='Sample')
pb_data$sample <- ifelse(pb_data$location=='Proximal',
                         paste0(pb_data$indv, '-','P'),
                         paste0(pb_data$indv,'-','D'))
sample_dat <- merge(prop_dat, pb_data[pb_data$layer=='basal', 
                                      c('sample','indv',quant_clin_traits)], by.x='Sample',by.y='sample')
sample_dat$nCells <- sample_dat$totCells*sample_dat$Proportion

for (i in unique(sample_dat[,grouping])) {
  sample_dat[nrow(sample_dat)+1,] <- sample_dat[sample_dat$Sample=='P0080-D' & sample_dat[,grouping]==i,]
  sample_dat[nrow(sample_dat),c('Sample','Location','nCells','totCells','Proportion')] <- c('P0080-P','Proximal',0,0,0)
}
#sample_dat <- head(sample_dat,152)
sample_dat$nCells <- as.numeric(sample_dat$nCells)
sample_dat[,grouping] <- factor(sample_dat[,grouping], levels=c('Superficial','Suprabasal','Replicating_suprabasal',
                                                                'Replicating_basal','Basal'))

pdf(file=paste0(output_dir, paste('../epi_cells/Epi_sampleCellCounts',grouping,'.pdf', sep='_')), width=8, height=5)
print(ggarrange(ggplot(sample_dat[sample_dat$Location=='Proximal',], aes(x=indv, y=nCells, fill=LayerRep)) + 
                  geom_bar(stat='identity',position='stack') + theme_minimal() + ggtitle('Proximal') + 
                  scale_fill_manual(values=rev(layer_cols_5)) + scale_y_continuous(limits=c(0,9000)),
                ggplot(sample_dat[sample_dat$Location=='Distal',], aes(x=indv, y=nCells, fill=LayerRep)) + 
                  geom_bar(stat='identity',position='stack') + theme_minimal() + ggtitle('Distal') + 
                  scale_fill_manual(values=rev(layer_cols_5)) + scale_y_continuous(limits=c(0,9000))
                , ncol=1, common.legend=T, legend='right'))
dev.off()

pdf(file=paste0(output_dir, paste('../epi_cells/Epi_sampleProportions',grouping,'.pdf', sep='_')), width=8, height=5)
print(ggarrange(ggplot(sample_dat[sample_dat$Location=='Proximal',], aes(x=indv, y=nCells, fill=LayerRep)) + 
                  geom_bar(stat='identity',position='fill') + theme_minimal() + ggtitle('Proximal') + 
                  scale_fill_manual(values=rev(layer_cols_5)),
                ggplot(sample_dat[sample_dat$Location=='Distal',], aes(x=indv, y=nCells, fill=LayerRep)) + 
                  geom_bar(stat='identity',position='fill') + theme_minimal() + ggtitle('Distal') + 
                  scale_fill_manual(values=rev(layer_cols_5)), ncol=1, common.legend=T, legend='right'))
dev.off()







for (i in 1:8) {
  dat <- merge(sample_dat[sample_dat$Cluster==i & sample_dat$Location=='Proximal',],
               sample_dat[sample_dat$Cluster==i & sample_dat$Location=='Distal',], 
               by='indv')[,c('Proportion.x','Proportion.y')]
  dat <- tail(dat,9)
  #print(dat)
  print(c(i,round(cor(cbind(as.numeric(dat$Proportion.x),as.numeric(dat$Proportion.y)))[1,2],2)))}




# Compute meta p-values for correlated variables
for (l in layers) {
  df <- lm.res[[l]]
  for (trait in quant_clin_traits) {
    df[,paste0('p.acat.', trait)] <- ACAT(rbind(df[,paste('p.p',trait, sep='.')], df[,paste('p.d',trait, sep='.')]))
  }
  lm.res[[l]] <- df
}





#### Pathway enrichment Analysis ####
background_genes <- genes

test_genes <- rownames(qlf.res[[paste0('SSc_HC')]][(qlf.res[[paste0('SSc_HC')]]$PValue<0.05) &
                                                     (qlf.res[[paste0('SSc_HC')]]$logFC>-40),])
ego <- enricher(gene          = test_genes,
                universe      = background_genes,
                pAdjustMethod = "fdr",
                pvalueCutoff  = 1,
                TERM2GENE     = annot_tf)

#View(as.data.frame(ego@result))
if(!is.null(ego)) {
  ego_df <- as.data.frame(ego@result)
  ego_df <- ego_df[ego_df$Count>2,]
  ego_df$p.adjust <- p.adjust(ego_df$pvalue)
  ego_df <- ego_df[,c(1,3:ncol(ego_df))]
  rownames(ego_df) <- NULL
}
View(ego_df)


kegg_enrich.bins.res <- lapply(setNames(levels(diff.bin), diff.bin_labs), function(b) {
  print(b)
  background_genes <- b.genes[[b]]
  test_genes <- rownames(qlf.res[[paste0('SSc_HC_',b)]][(qlf.res[[paste0('SSc_HC_',b)]]$p_adj<0.05) &
                                                          (qlf.res[[paste0('SSc_GERD_',b)]]$p_adj<0.05),])
  ego <- enricher(gene          = test_genes,
                  universe      = background_genes,
                  pAdjustMethod = "fdr",
                  pvalueCutoff  = 1,
                  TERM2GENE     = annot_kegg)
  
  #View(as.data.frame(ego@result))
  if(!is.null(ego)) {
    ego_df <- as.data.frame(ego@result)
    ego_df <- ego_df[ego_df$Count>2,]
    ego_df$p.adjust <- p.adjust(ego_df$pvalue)
    ego_df <- ego_df[,c(1,3:ncol(ego_df))]
    rownames(ego_df) <- NULL
    if(nrow(ego_df)>0) {
      return(ego_df)
    }
  }
})

fc <- gene.df[gene.df$SYMBOL %in% test_genes, 'logFC']
names(fc) <- gene.df[gene.df$SYMBOL %in% test_genes, 'ENTREZID']
if (nrow(as.data.frame(ego)) > 200) {
  cnetplot(ego, showCategory=sum(ego@result$p.adjust<0.05), layout='kk',
           color.params=list(foldChange=fc,category=alpha("#E5C494",0.8)), 
           cex.params=list(category_label=0.4, gene_label=0.55)) + scale_size(name='N genes') +
    theme(legend.title=element_text(size=8, face='bold'), plot.title=element_text(hjust=0.5, size=10,face='bold')) + 
    ggtitle(b) + scale_colour_gradient2(name='LogFC',low='blue3',high='red3')
}
return(ego)
})

for (b in layers) {
  pdf(file=paste0(output_dir, '/PseudoDEG_enrich_GO_',b,'_cnet.pdf'), width=6, height=6)
  print(attributes(go_enrich.res[[b]])$plot)
  dev.off()
  n_sig_pathways <- nrow(as.data.frame(go_enrich.res[[b]]))
  if (n_sig_pathways>0){
    pdf(file=paste0(output_dir, '/PseudoDEG_enrich_GO_',b,'_dotplot.pdf'), width=6, height=8)
    print(dotplot(go_enrich.res[[b]], showCategory=min(n_sig_pathways,20)) + ggtitle(b) + 
            theme(axis.text.y=element_text(size=8), legend.title=element_text(size=8, face='bold'),
                  plot.title=element_text(hjust=0.5, size=10,face='bold')) + 
            scale_color_gradient(name='q-value',low='rosybrown2',high='red3'))
    dev.off()
  }
  write.table(as.data.frame(go_enrich.res[[b]]), 
              paste0(output_dir, '/PseudoDEG_enrich_GO_',b,'.tsv'), quote=F, row.names=F,col.names=T, sep="\t")
}




