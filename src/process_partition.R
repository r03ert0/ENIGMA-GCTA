# Title     : Read and plot partitioned genetic variance outputs 
# Objective : Read and plot partitioned genetic variance outputs 
# Created by: ntraut, abiton
# Created on: 07/09/2018

library(tools)
library(data.table)
library(foreach)
library(ggplot2)
library(cowplot)
library(mgsub)


readHsqForPart <- function(path_res, prefix_res = '.hsq') {
  
  
  if (file.exists(path_res)) {
  
    fres <- list.files(path_res, pattern = prefix_res , full.names = T)
  
    if (length(fres == 0)) {
      # if no file found, try to suppress margin from pattern
      fres <- list.files(path_res, pattern = prefix_res , full.names = T)
      if (length(fres) == 0) stop(paste("Error: pattern",  prefix_res, "not found in", path_res))
    }
  
    # load files
    res <- sapply(fres, function(f) {
      # read log file instead of hsq file
      f2 <- paste0(file_path_sans_ext(f), '.log')
      lines <- readLines(f2)
      line_summary <- grep("Summary result of REML analysis:", lines)
      line_variance <- grep("Sampling variance/covariance", lines)
      line_hsq <- grep("Summary result of REML analysis has been saved", lines)
      line_n <- grep("individuals are in common in these files", lines)
      
      summary_table <- read.delim(text=lines[(line_summary+1):(line_variance-1)], stringsAsFactors=F)
      nbsamples <-  gsub('individuals are in common in these files.', '', lines[line_n])
  
      pval <- fread(paste0(file_path_sans_ext(f), '.hsq'), sep = '\t', stringsAsFactors = F, fill = T)[Source == 'Pval', ]$Variance
      
      cov_matrix <- read.delim(text=lines[(line_variance+1):(line_hsq-1)], header=F)
      cov_matrix <- as.matrix(cov_matrix[,1:nrow(cov_matrix)])
      
      data_dir <- sub("/data/.*", "/data", f)
      mgrm_file <- strsplit(grep("^--mgrm-bin", lines, perl=T, value = T), " ")[[1]][2]
      mgrm_file2 <- paste(data_dir, sub(".*/data/", "", mgrm_file), sep="/")
      grm_files <- readLines(mgrm_file2, warn=F)
      grm_files2 <- paste(data_dir, sub(".*/data/", "", grm_files), sep="/")
      nbsnp <- sapply(grm_files2, function(f) {
        N <- readBin(paste0(f, ".grm.N.bin"), n=1, what=numeric(0), size=4)
      })
      
      
      partition=sub("^.*/([^/]+)\\.[^/.]+$", "\\1", f)
      
      return(list(summary_table=summary_table, cov_matrix=cov_matrix, nbsnp=nbsnp, partition=partition, nbsamples = nbsamples, Pval = pval))
      
    }, simplify = F)
    return(res)
  } else {
    message(sprintf("path %s could not be found.", path_res ))
    return(NULL)
  }
}

process_partition <- function(path_res, prefix_res, output_file = NULL) {

  #### ==== Read original GCTA results  ==== ###
  res <- readHsqForPart(path_res = path_res, prefix_res = prefix_res)

  res.zscores <- sapply(res, function(r) {
    
    #print(r$partition)
    nbsamples <- r$nbsamples
    pvallrt <- r$Pval
    
    n <- r$nbsnp    # number of SNPs per group
    #n <- n[-length(n)]
    ng <- length(n)    # number of groups
    names(n) <- paste0("n", 1:ng)
    # variance explained by partitions
    v <- sapply(1:ng, function(i) subset(r$summary_table, Source == sprintf("V(G%d)", i))$Variance)
    w <- r$cov_matrix[1:ng, 1:ng] # covariance matrix of errors
    
    # expected v from htot without enrichment
    Ev <- sum(v)*n/sum(n)
    
    # theoretical standard errors of the estimations
    se.tot.th <- sqrt(sum(w))
    se.z.th <- sqrt((se.tot.th * n/sum(n))^2 -2*n/sum(n)*colSums(w) + diag(w))
    
    # theoretical p-values
    z.obs.th <- (v - Ev) / se.z.th
    names(z.obs.th) <- paste0("z", 1:ng)
    #p.z.obs.th <- 1 - abs(0.5 - pnorm(z.obs.th)) * 2 #two-sided
    p.z.obs.th <- pnorm(z.obs.th, lower.tail = FALSE)  #one-sided
    names(p.z.obs.th) <- paste0("p", 1:ng)
    
    # heritability values
    sum.melt=melt(r$summary_table, id.vars=1)
    h <- sum.melt$value
    names(h) <- sub("_Variance", "", paste(sum.melt$Source, sum.melt$variable, sep="_"))
    h <- h[grep("/Vp", names(h))]
    
    # compute V(Gi)/Vg
    ratios <- h[1:ng] / sum(h[1:ng])
    names(ratios) <- sprintf("V(G%d)/Vg", 1:ng)
    se.ratios <- abs(ratios) * sqrt(t(diag(w)/t((v)^2)) + (se.tot.th/sum(v))^2 - 2*colSums(w)/(v*sum(v)))
    names(se.ratios) <- paste0(names(ratios), "_SE")
    # if we have an a priori that E(v)=sum(v)*n/sum(n) (not recommended for standard errors)
    # derived z-scores differ from z.obs.th because of different precision in values reported by GCTA
    se.ratios2 <- n/sum(n) * sqrt(t(diag(w)/t((Ev)^2)) + (se.tot.th/sum(v))^2 - 2*colSums(w)/(Ev*sum(v)))
    
    return(c(as.list(h), ratios, se.ratios, n, z.obs.th, p.z.obs.th, nind = nbsamples, Plrt = pvallrt))
  }, simplify=F)
  
  allnames <- unique(unlist(lapply(res.zscores, names)))
  res.zscores <- lapply(res.zscores, 
                        function(x, allnames) {
                          xx <- x[allnames] 
                          xx[which(sapply(xx,length) ==0)] <- NA
                          names(xx) <- allnames
                          xx
                        }, 
                        allnames = allnames)
  pvals <- rbindlist(res.zscores)
  pvals <- cbind(partition=sapply(res, "[[", "partition"), pvals)
  
  i2id <- structure(.Names = as.character(1:length(res[[1]]$nbsnp)),
                    .Data = basename(names(res[[1]]$nbsnp)))
  
  setnames(pvals, names(pvals), mgsub::mgsub(string  = names(pvals), pattern = names(i2id), replacement = i2id))
  setnames(pvals, names(pvals), gsub('^z', 'Z_', names(pvals)))
  setnames(pvals, names(pvals), gsub('^p', 'Penrich_', names(pvals)))
  setnames(pvals, names(pvals), gsub('^n', 'n_', names(pvals)))
  setnames(pvals, names(pvals), gsub('^V\\(G', 'V', names(pvals)))
  setnames(pvals, names(pvals)[grep('Sum', names(pvals), invert = T)], gsub('\\)', '', names(pvals)[grep('Sum', names(pvals), invert = T)]))
  setnames(pvals, names(pvals), gsub('Penrich_artition', 'partition', names(pvals)))
  setnames(pvals, names(pvals), gsub('Penrich_lrt', 'Plrt', names(pvals)))
  
  pvals[, Phenotype := basename(gsub('\\.', '/', partition))]
  pvals[, snpgroup := gsub('/', '\\.', dirname(gsub('\\.', '/', partition)))]

    # flattens GCTA pvalues 
  if (!is.null(pvals$Plrt)) {
    pvalbygroup <- split(pvals$Plrt, pvals$snpgroup)
    names(pvalbygroup) <-  gsub('.*\\.', '',   names(pvalbygroup))
    names(pvalbygroup) <-  mgsub::mgsub(names(pvalbygroup), pattern = names(i2id), replacement = i2id)
    names(pvalbygroup) <- paste0('Plrt_',  names(pvalbygroup))
  } else {
    pvalbygroup <- NULL
  }
  pvals <- pvals[grep(paste('\\.',2:length(res[[1]]$nbsnp), collapse='|', sep = ''), snpgroup, invert = T),]
  
  for (i in names(pvalbygroup)) {
    pvals[[i]] <- pvalbygroup[[i]]
  }
  pvals[['Plrt']] <- NULL
    
  #setnames(pvals, 'Pval', 'Plrt')
  pvals <- pvals[, c('partition', 'snpgroup', 'Phenotype', setdiff(names(pvals),c('partition', 'snpgroup', 'Phenotype'))), with = F]
  
  if (!is.null(output_file))
    fwrite(pvals, file = output_file, sep = '\t')

  return(pvals)
}




plot_hsq_partition <- function(hsqsummary, normByPropSNPs = FALSE, nbSNPs = NULL, 
                                             outpdf = NULL, pdfwidth = 15, pdfheight = 8, width = 0.5,
                                             ignore_snpgroup = NULL) {
  
  
  if (is.character(hsqsummary)) {
    hsq <- fread(hsqsummary, header = T, sep = '\t')
  } else {
    hsq <- hsqsummary
  }
  
  #hsq <- hsq[, setdiff(names(hsq),'n_ind'), with = F]
  hsq <- melt(hsq, id.vars = names(hsq)[1:3])
  
  
  if (!is.null(ignore_snpgroup)){
    hsq <- hsq[!(snpgroup %in% ignore_snpgroup),]
  }
  
  hsq[, variable := as.character(variable)]
  ses <- hsq[grep('Vp_se', variable, ignore.case = T, invert = F), ]
  hsq[match(paste(gsub('_SE', '', ses$variable),ses$snpgroup,ses$Phenotype),paste(variable,snpgroup,Phenotype)), se := ses$value]
  hsq <- hsq[grep('Vp_se', variable, ignore.case = T, invert = T), ]
  plrt <- hsq[grep('Plrt', variable, ignore.case = T, invert = F), ]
  hsq[match(paste(paste0(gsub('Plrt_', 'V', plrt$variable),'/Vp'),plrt$snpgroup,plrt$Phenotype),paste(variable,snpgroup,Phenotype)), Plrt := plrt$value]
  hsq <- hsq[grep('Plrt', variable, ignore.case = T, invert = T), ]
  
  vgs <- hsq[grep('/Vg$', variable, ignore.case = T, invert = F), ]
  hsq[match(paste(gsub('/Vg', '/Vp', vgs$variable),vgs$snpgroup,vgs$Phenotype),
            paste(variable,snpgroup,Phenotype)), 
      vgvg := vgs$value]
  
  vgses <- hsq[grep('/Vg_SE$', variable, ignore.case = T, invert = F), ]
  hsq[match(paste(gsub('/Vg_SE', '/Vp', vgses$variable),vgses$snpgroup,vgses$Phenotype),
            paste(variable,snpgroup,Phenotype)), 
      vgvg_se := vgses$value]
  
  penrich <- hsq[grep('Penrich', variable, ignore.case = T, invert = F), ]
  hsq[match(paste(paste0(gsub('Penrich_', 'V', penrich$variable),'/Vp'),penrich$snpgroup,penrich$Phenotype),paste(variable,snpgroup,Phenotype)), Penrich := penrich$value]
  hsq <- hsq[grep('Penrich', variable, ignore.case = T, invert = T), ]
  
  zscores <- hsq[grep('Z_', variable, ignore.case = T, invert = F), ]
  hsq[match(paste(paste0(gsub('Z_', 'V', zscores$variable),'/Vp'),zscores$snpgroup,zscores$Phenotype),
            paste(variable,snpgroup,Phenotype)), Zscores := zscores$value]
  hsq <- hsq[grep('Z_', variable, ignore.case = T, invert = T), ]
  
  hsq <- hsq[grep('/Vg$', variable, ignore.case = T, invert = T), ]
  hsq <- hsq[grep('/Vg_SE$', variable, ignore.case = T, invert = T), ]
  
  nbind <- structure(.Data = hsq[variable == 'n_ind',]$value, .Names = hsq[variable == 'n_ind',]$Phenotype)
  hsq <- hsq[variable != 'n_ind', ]
  
  nbsnps <- structure(.Data = hsq[grep('n_', variable),]$value, .Names = gsub('n_', '', hsq[grep('n_', variable),]$variable))
  nbsnps <- nbsnps[unique(names(nbsnps))]
  levvar <- names(nbsnps)
  
  hsq <- hsq[grep('n_', variable, invert = T), ]
  hsq <- hsq[grep('Z_', variable, invert = T), ]
  hsq[, variable := gsub(' ', '_', as.character(variable))]
  

  hsqsave <- hsq
  hsq <- hsq[grep('^V', variable),]
  
  hsq[grep('/Vp', variable), vgvp := as.numeric(gsub(' .*', '',value))]
  hsq[grep('/Vp', variable), vgvp_se := as.numeric(se)]
  
  hsq <- hsq[, partition := gsub('V|/Vp|/Vg', '', variable)]
  hsq[, nbsnp := nbsnps[partition]]
  
  hsq <- hsq[variable != 'Sum_of_V(G)/Vp',]
  
  if (normByPropSNPs) {
    hsq[, f := as.numeric(nbsnp)/sum(as.numeric(nbsnps))]
    hsq[, vgvgdivf := as.numeric(vgvg)/f]
    hsq[, vgvgdivf_se := as.numeric(vgvg_se)/f]
  }
  
  # add nb SNPs to group name (such that it will show in the barplot legend with ggplot2)
  hsq[, variable := paste0(variable, ' (', nbsnp, ' SNPs)')]
  hsq[, variable := factor(variable, levels = sapply(levvar,  function(x) hsq$variable[grep(pattern = x, x = hsq$variable)[1]]))]
  
  hsq[, Phenotype := gsub('log10', '', Phenotype)]
  
  hsq <- hsq[grep ("\\.1", snpgroup),]
  hsq[, margin := gsub('.*-margin','margin', gsub("\\.1", "", snpgroup))]
  hsq[, partitions := partition]
  
  
  hsq[, FDR_enrich := p.adjust(Penrich, method = 'BH'), by = c('variable')]
  
  hsq_sub <- hsq
  hsq_sub[, Penrich := as.numeric(Penrich)]
  hsq_sub[, star := ""]
  hsq_sub[FDR_enrich <= .05, star :=  "*"]
  hsq_sub[FDR_enrich <= .01, star :=  "**"]
  hsq_sub[FDR_enrich <= .001, star :=  "***"]
  
 
  hsq_sub[, variable := factor(variable, levels = sapply(levvar,  function(x) hsq_sub$variable[grep(pattern = x, x = hsq_sub$variable)[1]]))]
  
  if (!normByPropSNPs) {
    gg <- ggplot(hsq_sub, aes(x=Phenotype, y=vgvp, fill=variable)) + #, group=variable)) + 
      geom_bar(position=position_dodge2(1), stat="identity", 
               width = width,
               colour="black",  # Use black outlines,
               size=.3) +
      ylab(if (!normByPropSNPs) 'VG/VP  ± s.e' ) +
      xlab('Phenotype') + 
      scale_fill_brewer(palette = 'Set3') +
      theme_bw( ) +
      theme( axis.text.x = element_text(angle = 90, size = 11, vjust = 0.5, hjust = 1), 
             axis.title.y = element_text(size = 10),
             axis.text.y = element_text(size = 10),
             axis.text.x.top = element_text(size = 10),
             strip.text.x = element_text(size=11), 
             legend.text = element_text(size = 10))
    
    
    gg <- gg + geom_errorbar(aes(ymin=vgvp-vgvp_se, ymax=vgvp+vgvp_se, group=variable),
                             size=.2,    # Thinner lines
                             width=width, 
                             position=position_dodge2(1))
    
    #gg <- gg + geom_text(aes(label=star), colour="black", vjust=0, size=5)
    
  } else {
    gg <- ggplot(hsq_sub, aes(x=Phenotype, y=vgvgdivf, fill=variable, group = variable)) + #, group=variable)) + 
      geom_bar(position=position_dodge2(1), stat="identity", 
               width = width,
               colour="black",  # Use black outlines,
               size=.3) +
      ylab( '%VG  divided by %SNPs') +
      xlab('Phenotype') + 
      scale_fill_brewer(palette = 'Set3') +
      coord_cartesian(ylim = c(0, max(hsq_sub$vgvgdivf+hsq_sub$vgvgdivf_se) + 0.5)) + 
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
    
    
    gg <- gg + geom_errorbar(aes(ymin=vgvgdivf-vgvgdivf_se, ymax=vgvgdivf+vgvgdivf_se, group=variable),
                             size=.2,    # Thinner lines
                             width=width, 
                             position=position_dodge2(1))
    gg <- gg  + geom_text(aes(x = Phenotype, y=vgvgdivf+vgvgdivf_se, label=star, group = variable), colour="black", size=3, vjust = -0.5, 
                           position= position_dodge2(width = ((length(unique(hsq_sub$variable))-1)/(length(unique(hsq_sub$variable))))), hjust = 0.5, #label.r = unit(0, "lines"), 
                           #label.size = 0, 
                          show.legend = FALSE)
    
  }
  

  
  if (!is.null(outpdf)) {
    pdf(outpdf, width = pdfwidth, height = pdfheight)
    print(gg)
    dev.off()
  }
  
  return(list(hsq, gg))
  
}



write_plot_hsq_partition <- function(hsq_outdir, output_dir, type_partition, group_partition) {
  
    
    partition <- type_partition
    group <- group_partition
    path_res <- file.path(hsq_outdir, paste0('hsq-', partition))
    
    output_file <- file.path(output_dir, paste0('enrichment_', partition, if (group != partition) paste0('_',group) else '', '.txt'))
    res <- process_partition(path_res = path_res,  prefix_res = paste0(group,'.*.hsq'),
                             output_file = output_file)

    if (!is.null(res)) {
      res[, Phenotype:= gsub('_freesurfer|log10', '', Phenotype)]

      #pdf(paste0(output_dir,'hsq-', partition, '-', gsub('\\.', '', group), '.pdf'), 11, 5)
      gg1 <- plot_hsq_partition(hsqsummary = res[grep('left|right', Phenotype, invert = T),],
                                normByPropSNPs = FALSE,
                                width = 0.7)
      gg2 <- plot_hsq_partition(hsqsummary = res[grep('left|right', Phenotype, invert = T),],
                                normByPropSNPs = TRUE,
                                width = 0.7)
      # gg3 <- plot_hsq_partition(hsqsummary = res[grep('left|right', Phenotype, invert = F),],
      #                           normByPropSNPs = FALSE,
      #                           width = 0.7)
      # gg4 <- plot_hsq_partition(hsqsummary = res[grep('left|right', Phenotype, invert = F),],
      #                           normByPropSNPs = TRUE,
      #                           width = 0.7)
      #dev.off()

      if (group == partition)
      group <- ''
      pdf(file.path(output_dir, paste0('enrichment_', partition, if (group != '') '_' else '', gsub('\\.', '', group), '.pdf')), 11, 8)
      print(plot_grid(gg1[[2]] + xlab('') +  theme(axis.text.x =element_blank()) ,
                      gg2[[2]],
                      labels = c("A", "B"), align = "v",  ncol = 1, rel_heights = c(1,1.3)))
      dev.off()
    }
  
}

# output_dir <- '/Volumes/cinq/rto/ukb-annex/data/derived/ukb-22110_maf001/06.hsq_FreeSurfer/hsq-summary/'
# hsq_outdir <- '/Volumes/cinq/rto/ukb-annex/data/derived/ukb-22110_maf001/06.hsq_FreeSurfer/'
# 
# 
# allres <- foreach (partition=c('genic',
#                                'cnsexpression', 
#                                'neurodev', 'maf', 
#                                rep('genic', 10)),
#                    group=c("updown-margin20-50", 
#                            'cnsexpression', 
#                            'neurodev', 'maf', 
#                            "genic-exon1intron1",  
#                            'genic-margin10', 'genic-margin20', 'genic-margin30', 'genic-margin40', 'genic-margin50',
#                            "updown-margin20\\.", "updown-margin30", "updown-margin40", "updown-margin50")
# ) %do% {
#   
#   path_res <- paste0(hsq_outdir,'hsq-', partition)
#   
#   output_file <- paste0(output_dir, 'enrichment_', paste0(partition, if (group != partition) paste0('_',group) else ''), '.txt')
#   res <- process_partition(path_res = path_res,  prefix_res = paste0(group,'.*.hsq'), 
#                            output_file = output_file)
#   
#   if (!is.null(res)) {
#     res[, Phenotype:= gsub('_freesurfer|log10', '', Phenotype)]
#     
#     #pdf(paste0(output_dir,'hsq-', partition, '-', gsub('\\.', '', group), '.pdf'), 11, 5)
#     gg1 <- plot_hsq_partition(hsqsummary = res[grep('left|right', Phenotype, invert = T),],
#                               normByPropSNPs = FALSE,
#                               width = 0.7)
#     gg2 <- plot_hsq_partition(hsqsummary = res[grep('left|right', Phenotype, invert = T),],
#                               normByPropSNPs = TRUE,
#                               width = 0.7)
#     gg3 <- plot_hsq_partition(hsqsummary = res[grep('left|right', Phenotype, invert = F),],
#                               normByPropSNPs = FALSE,
#                               width = 0.7)
#     gg4 <- plot_hsq_partition(hsqsummary = res[grep('left|right', Phenotype, invert = F),],
#                               normByPropSNPs = TRUE,
#                               width = 0.7)
#     #dev.off()
#     
#     pdf(paste0(output_dir,'hsq_', partition, '-', gsub('\\.', '', group), '.pdf'), 11, 8)
#     print(plot_grid(gg1[[2]] + xlab('') +  theme(axis.text.x =element_blank()) ,
#                     gg2[[2]],
#                     labels = c("A", "B"), align = "v",  ncol = 1, rel_heights = c(1,1.3)))
#     dev.off()
#     
#   }
#   
#   return(res)
# }
# 
# names(allres) <- c("updown-margin20-50", #"updown-margin20-40", 
#                    'cnsexpression', 'neurodev', 'maf', 
#                    "genic-exon1intron1",  
#                    'genic-margin10', 'genic-margin20', 'genic-margin30', 'genic-margin40', 'genic-margin50',
#                    "updown-margin20", "updown-margin30", "updown-margin40", "updown-margin50")
# 
# WriteXLS::WriteXLS(x = lapply(allres[!is.null(allres)], function(x) x[, -c(1:2)]), 
#                    ExcelFileName = paste0(output_dir,'hsq_partitions.xls'), row.names = F)

