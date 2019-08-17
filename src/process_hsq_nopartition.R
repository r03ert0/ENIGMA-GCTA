# Title     : Read and plot partitioned genetic variance outputs 
# Objective : Read and plot partitioned genetic variance outputs 
# Created by: ntraut, abiton
# Created on: 07/09/2018

library(tools)
library(data.table)
library(foreach)
library(ggplot2)
library(cowplot)

readHsq <- function(path_res, prefix_res = '.hsq') {
  
  fres <- list.files(path_res, pattern = prefix_res , full.names = T)
  
  if (length(fres) == 0) {
    # if no file found, try to suppress margin from pattern
    fres <- list.files(path_res, pattern = gsub('-margin..|-margin.|-margin..', '*', prefix_res) , full.names = T)
    print(fres)
    if (length(fres) == 0) stop(paste("Error: pattern",  prefix_res, "not found in", path_res))
  }
  
  # load files
  res <- lapply(fres, function(f) {
    res1 <- read.delim(f, sep = '\t') # stops at empty line before LRT information
    phenosel <- sub('.*\\.(.*?).hsq', '\\1', basename(f)) # to check
    res2 <- melt(res1, id.vars=1)
    row.names(res2) = sub('_variance', '', paste0(res2$Source,'_',tolower(res2$variable)))
    res3 <- as.data.table(t(res2[, "value", drop=F]))
    res3$pheno = phenosel
    return(res3)
  })
  
  
  # add missing names to results without max number of rows
  if (length(which(unlist(lapply(res, length)) != max(unlist(lapply(res, length))))) > 0) {
    coln <- names(res[which(unlist(lapply(res, length)) == max(unlist(lapply(res, length))))[1]][[1]])
    
    res[which(unlist(lapply(res, length)) != max(unlist(lapply(res, length))))] <-
      lapply(res[which(unlist(lapply(res, length)) != max(unlist(lapply(res, length))))],
             function(x, coln) {
               x <- cbind(x, as.data.table(t(structure(.Names = setdiff(coln,colnames(x)),
                                                       .Data = rep(NA, length(setdiff(coln,colnames(x))))))))
               x[, coln, with = F]
             }, coln = coln)
  }
  res <- as.data.table(do.call(rbind, res))
  res[, file := basename(fres)]
  return(res)
}



plot_hsq <- function(hsq, normByPropSNPs = FALSE, nbSNPs = NULL, colorby = NULL, 
                     outpdf = NULL, pdfwidth = 15, pdfheight = 8, width = 0.5,
                     pattern_pheno_ignore = 'left|right'
                     ) {
  

  if (is.character(hsq)) {
    hsq <- fread(hsq, sep = '\t', stringsAsFactors = F)
  }
  
  hsq <- hsq[grep(pattern_pheno_ignore, pheno, invert = T),]
  #hsq <- hsq[, setdiff(names(hsq),'n_ind'), with = F]
  hsq <- hsq[,c("pheno", "V(G)/Vp",  "V(G)/Vp_se", "Pval", "n", colorby), with = F]
  setnames(hsq, "V(G)/Vp", "vgvp")
  setnames(hsq, "V(G)/Vp_se", "vgvp_se")
  
  if (!is.factor(hsq$pheno))
    hsq[, Phenotype := gsub('log10|_freesurfer', '', pheno)]
  else 
    hsq[, Phenotype := pheno]
  hsq[, FDR_enrich := p.adjust(Pval, method = 'BH'), by = c('pheno')]
  

  hsq_sub <- hsq
  hsq_sub[, Penrich := as.numeric(Pval)]
  hsq_sub[, star := ""]
  # hsq_sub[Penrich <= .05, star :=  "*"]
  # hsq_sub[Penrich <= .01, star :=  "**"]
  # hsq_sub[Penrich <= .001, star :=  "***"]
  hsq_sub[FDR_enrich <= .05, star :=  "*"]
  hsq_sub[FDR_enrich <= .01, star :=  "**"]
  hsq_sub[FDR_enrich <= .001, star :=  "***"]
  
 
  if (!is.null(colorby)) {
    gg <- ggplot(hsq_sub, aes(x=Phenotype, y=vgvp, fill = get(colorby))) + #, group=variable)) + 
      geom_bar(position=position_dodge2(1), stat="identity", 
               width = width, 
               colour="black",  # Use black outlines,
               size=.3) +
      ylab('VG/VP  ± s.e' ) +
      xlab('Phenotype') + 
      scale_fill_brewer(palette = 'Set3') +
      #scale_fill_hue(name="Supplement type", # Legend label, use darker colors
      #               breaks=c("OJ", "VC"),
      #               labels=c("Orange juice", "Ascorbic acid")) +
      #facet_wrap(~margin, nrow = 3)+ 
      theme_bw( ) +
      theme( axis.text.x = element_text(angle = 90, size = 11, vjust = 0.5, hjust = 1), 
             axis.title.y = element_text(size = 10),
             axis.text.y = element_text(size = 10),
             axis.text.x.top = element_text(size = 10),
             strip.text.x = element_text(size=11), 
             legend.text = element_text(size = 10), legend.title = element_blank())
    
    
    gg <- gg + geom_errorbar(aes(ymin=vgvp-vgvp_se, ymax=vgvp+vgvp_se),
                             size=.2,    # Thinner lines
                             width=width, 
                             position=position_dodge2(1))
    
    #gg <- gg + geom_text(aes(label=star), colour="black", vjust=0, size=5)
    
  } else {
    gg <- ggplot(hsq_sub, aes(x=Phenotype, y=vgvp)) + #, group=variable)) + 
      geom_bar(position=position_dodge2(1), stat="identity", 
               width = width, fill = "aquamarine4", 
               colour="black",  # Use black outlines,
               size=.3) +
      ylab('VG/VP  ± s.e' ) +
      xlab('Phenotype') + 
      #scale_fill_brewer(palette = 'Set3') +
      #scale_fill_hue(name="Supplement type", # Legend label, use darker colors
      #               breaks=c("OJ", "VC"),
      #               labels=c("Orange juice", "Ascorbic acid")) +
      #facet_wrap(~margin, nrow = 3)+ 
      theme_bw( ) +
      theme( axis.text.x = element_text(angle = 90, size = 11, vjust = 0.5, hjust = 1), 
             axis.title.y = element_text(size = 10),
             axis.text.y = element_text(size = 10),
             axis.text.x.top = element_text(size = 10),
             strip.text.x = element_text(size=11), 
             legend.text = element_text(size = 10))
    
    
    gg <- gg + geom_errorbar(aes(ymin=vgvp-vgvp_se, ymax=vgvp+vgvp_se),
                             size=.2,    # Thinner lines
                             width=width, 
                             position=position_dodge2(1))
    
    #gg <- gg + geom_text(aes(label=star), colour="black", vjust=0, size=5)
  }
  
  if (!is.null(outpdf)) {
    pdf(outpdf, width = pdfwidth, height = pdfheight)
    print(gg)
    dev.off()
  }  
  
  return(list(hsq, gg))
  
}




