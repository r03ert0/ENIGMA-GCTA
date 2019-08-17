library(GenomicRanges)
suppressMessages(library(data.table))
suppressMessages(library(foreach))
suppressMessages(library(ggplot2))
library(ggbio)
suppressMessages(library(qqman))

plotManhattan <- function(gr, pval_column = 'P', label_column = 'SNP', pval_threshold = 1e-4, phenotype = '') {
  
  ## manhattan plot
  dd <-   gr
  dd$pval <-  -log10(mcols(gr)[[pval_column]])
  #dd <- keepSeqlevels(dd, as.character(c(1:19,'X','Y')), pruning.mode = 'coarse')
  #seqlengths(dd)
  dd$label <-  mcols(dd)[[label_column]]
  dd[is.na(dd$pval) | dd$pval < -log10(pval_threshold)]$label <- NA
  
  gg <- plotGrandLinear(dd, 
                        aes(y = pval), 
                        cutoff.color = "blue", cutoff.size = 0., cex = 0.7) + 
    ylab('-log10( p-value )') + #+ cutoff = 3, 
    #ggforce::facet_zoom (x = seqnames == '15', zoom.size = 1) +  
    geom_text(aes(label = label), size = 3, check_overlap = F, vjust = 0, nudge_y = 0.05) +
    ggtitle(phenotype)
  
  return(gg)
  
}


plotGWAS <- function(res_path) {

  gwasfiles <- list.files(res_path, 
                        pattern = '.linear', 
                        full.names = T)

  names(gwasfiles) <- gsub('all.|.assoc.linear', '', basename(gwasfiles))
  
  
  
  allgwas <- foreach( gwasfile = gwasfiles,
                      phenotype = names(gwasfiles)) %do% {
                        
                        gwasres<- fread(gwasfile)
                        ## Filter out results for covariables
                        gwasres <- gwasres[TEST == 'ADD',]
                        return(gwasres)
                      }
  names(allgwas) <- names(gwasfiles)
  
  tutu <- foreach( gwasres = allgwas,
                   phenotype = names(allgwas)) %do% {
                     
                     if (!all(is.na(gwasres$P))) {
                       pdf(paste0(res_path,'',gsub('log10|_freesurfer', '', phenotype),'.pdf'), 12, 8)
                       ## transform into GRanges object
                       gwasresgr <- copy(gwasres)
                       gwasresgr[, start := BP]
                       gwasresgr[, end := BP]
                       gwasresgr <- GenomicRanges::makeGRangesFromDataFrame(gwasresgr, keep.extra.columns = T)
                       ## manhattan plot
                       print(plotManhattan(gwasresgr, phenotype = phenotype))
                       ## qq-plot
                       qq(gwasres$P, col = alpha('black',0.5), cex = 0.5, main = phenotype)
                       dev.off()
                       
                     }
                   }

}


#res_path <- "/Volumes/cinq/rto/ukb-annex/data/derived/ukb-22110//05.gwas/gwas-pruned/"
#res_path <- "/Volumes/cinq-1/rto/ukb-annex/data/derived/IMAGEN123//05.gwas/gwas-pruned/"
#plotGWAS(res_path)
