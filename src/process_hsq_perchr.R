# Writes, and plots genetic variance partitioned per chromosome.

library(tools)
library(data.table)
library(foreach)
library(ggplot2)
library(cowplot)

readHsq_perchr <- function(path_res, prefix_res = '.hsq') {
  
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
    
    summary_table <- read.delim(text=lines[c((line_summary+1):(line_variance-1), line_n)], stringsAsFactors=F)
    nbsamples <-  gsub('individuals are in common in these files.', '', read.delim(text=lines[line_n], stringsAsFactors=F, sep = ' ', header = F, check.names = F)[1,1])
    
    pval <- fread(paste0(file_path_sans_ext(f), '.hsq'), sep = '\t', stringsAsFactors = F, fill = T)[Source == 'Pval', ]$Variance
    
    cov_matrix <- read.delim(text=lines[(line_variance+1):(line_hsq-1)], header=F)
    cov_matrix <- as.matrix(cov_matrix[,1:nrow(cov_matrix)])
    
    data_dir <- sub("/data/.*", "/data", f)
    grm_file <- strsplit(grep("Reading the GRM from \\[", lines, perl=F, value = T), " ")[[1]][5]
    grm_file <- gsub("\\].|\\[", "", grm_file)
    grm_file <- paste(data_dir, gsub(".*/data/|.grm.bin", "", grm_file), sep="/")
    
    nbsnp <- sapply(grm_file, function(f) {
      N <- readBin(paste0(f, ".grm.N.bin"), n=1, what=numeric(0), size=4)
    })
    
    partition=sub("^.*/([^/]+)\\.[^/.]+$", "\\1", f)
    
    return(list(summary_table=summary_table, cov_matrix=cov_matrix, 
                nbsnp=nbsnp, partition=partition, nbsamples = nbsamples, Pval = pval))
    
  }, simplify = F)
  return(res)
}

process_perchr <- function(path_res, prefix_res, make_plot = FALSE) {
  
  #### ==== Read original GCTA results  ==== ###
  res <- readHsq_perchr(path_res = path_res, prefix_res = prefix_res)
  
  res.zscores <- sapply(res, function(r) {
    
    #print(r$partition)
    nbsamples <- r$nbsamples
    pvallrt <- r$Pval
    
    n <- r$nbsnp    # number of SNPs per group
    #n <- n[-length(n)]
    ng <- length(n)    # number of groups
    names(n) <- paste0("n", 1:ng)
    
    # heritability values
    sum.melt=melt(r$summary_table, id.vars=1)
    h <- sum.melt$value
    names(h) <- sub("_Variance", "", paste(sum.melt$Source, sum.melt$variable, sep="_"))
    h <- h[grep("/Vp", names(h))]
    
    return(c(as.list(h), n = n, nind = nbsamples, Plrt = pvallrt))
  }, simplify=F)
  
  allnames <- unique(unlist(lapply(res.zscores, names)))
  #allnames[allnames == 'n'] <- 'f'
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
  
  #setnames(pvals, 'n.n1', 'n')
  i2id <- structure(.Names = as.character(1:length(res[[1]]$nbsnp)),
                    .Data = basename(names(res[[1]]$nbsnp)))
  
  setnames(pvals, names(pvals), mgsub::mgsub(string  = names(pvals), pattern = names(i2id), replacement = i2id))
  #setnames(pvals, names(pvals), gsub('^V\\(G', 'V', names(pvals)))
  
  pvals[, Phenotype := basename(gsub('\\.', '/', partition))]
  pvals[, snpgroup := gsub('/', '\\.', dirname(gsub('\\.', '/', partition)))]
  
  
  # pvals <- pvals[grep(paste('\\.',2:length(res[[1]]$nbsnp), collapse='|', sep = ''), snpgroup, invert = T),]
  # 
  # for (i in names(pvalbygroup)) {
  #   pvals[[i]] <- pvalbygroup[[i]]
  # }
  # pvals[['Plrt']] <- NULL
  
  #setnames(pvals, 'Pval', 'Plrt')
  pvals <- pvals[, c('partition', 'snpgroup', 'Phenotype', setdiff(names(pvals),c('partition', 'snpgroup', 'Phenotype'))), with = F]
  
  #fwrite(pvals, file = output_file, sep = '\t')
  
  return(pvals)
}



plot_perchr <- function() {
  
  plotperchr <- foreach (pheno = unique(res$Phenotype)) %do% {
    
    print(pheno)
    
    hsqsub <- res[Phenotype == pheno,]
    setnames(hsqsub, c('V(G)/Vp',  'V(G)/Vp_SE'), c('vgvp', 'vgvp_se'))
    gg <- 
      ggplot(hsqsub, 
             aes(x = nbsnpsperchr[hsqsub$snpgroup],
                 y = vgvp)) + #, colour = snpgroup)) + 
      geom_point(size = 2) + geom_smooth(method=lm, alpha = 0.2, color = alpha('blue', 0.3)) + 
      ggrepel::geom_text_repel(aes(x = nbsnpsperchr[hsqsub$snpgroup],
                                   y = vgvp, 
                                   label = snpgroup)) + 
      geom_errorbar(aes(ymin=vgvp-vgvp_se, ymax=vgvp+vgvp_se),
                    size=.3,    # Thinner lines
                    width=.2,
                    position=position_dodge(.9)) +
      scale_x_continuous(breaks = round(seq(0, max(nbsnpsperchr), by = 10^3),1)) +
      ylab('VG/VP  ± s.e.') + xlab('#SNPs') + 
      # geom_text(x = 12000, y = max(hsqsub$vgvp)-0.005, 
      #           label = lm_eqn(data.frame(x =nbsnpsperchr[hsqsub$snpgroup],y = hsqsub$vgvp)), 
      #           parse = TRUE, size = 4) +
      ggtitle(pheno, subtitle = paste0('(r = ',signif(cor(nbsnpsperchr[hsqsub$snpgroup],hsqsub$vgvp), 3),
                                       ', p-value = ',signif(cor.test(nbsnpsperchr[hsqsub$snpgroup],hsqsub$vgvp)$p.value, 3), ')')) + 
      theme_bw( ) +
      theme( axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5, hjust = 1), 
             axis.text.y = element_text(size = 10),
             axis.title.y = element_text(size = 15),
             legend.position = 'none',
             axis.text.x.top = element_text(size = 15),
             legend.text = element_text(size = 15),
             strip.text.x = element_text(size=20))
    
    return(gg)
    
  }
  
}




write_plot_hsq_perchr <- function(hsq_outdir, output_dir) {
  
  path_res <- file.path(hsq_outdir, 'hsq-perchr')
  print(path_res)
  res <- process_perchr(path_res = path_res,  prefix_res = '^chr.*.hsq', 
                        make_plot = FALSE)
  
  res[, Phenotype:= gsub('_freesurfer|log10', '', Phenotype)]
  res[, partition:= gsub('_freesurfer|log10', '', partition)]
  res[, dataset := gsub('06.hsq.*|.*derived', '', hsq_outdir)]
  setnames(res, "n.nall-0.025-chr1", "n")
  
  output_file <- file.path(output_dir, 'hsq-perchr.txt')
  fwrite(res, file = output_file, sep = '\t')

  nbsnpsperchr <- structure(.Data = unique(res[, c('snpgroup', 'n')])$n,
                            .Names = unique(res[, c('snpgroup', 'n')])$snpgroup)
  
  plotperchr <- foreach (pheno = unique(res$Phenotype)) %do% {
    hsqsub <- res[Phenotype == pheno,]
    setnames(hsqsub, c('V(G)/Vp',  'V(G)/Vp_SE'), c('vgvp', 'vgvp_se'))
    gg <- 
      ggplot(hsqsub, 
             aes(x = nbsnpsperchr[hsqsub$snpgroup],
                 y = vgvp)) + #, colour = snpgroup)) + 
      geom_point(size = 2) + geom_smooth(method=lm, alpha = 0.2, color = alpha('blue', 0.3)) + 
      ggrepel::geom_text_repel(aes(x = nbsnpsperchr[hsqsub$snpgroup],
                                   y = vgvp, 
                                   label = snpgroup)) + 
      geom_errorbar(aes(ymin=vgvp-vgvp_se, ymax=vgvp+vgvp_se),
                    size=.3,    # Thinner lines
                    width=.2,
                    position=position_dodge(.9)) +
      scale_x_continuous(breaks = round(seq(0, max(nbsnpsperchr), by = 10^3),1)) +
      ylab('VG/VP  ± s.e.') + xlab('#SNPs') + 
      # geom_text(x = 12000, y = max(hsqsub$vgvp)-0.005, 
      #           label = lm_eqn(data.frame(x =nbsnpsperchr[hsqsub$snpgroup],y = hsqsub$vgvp)), 
      #           parse = TRUE, size = 4) +
      ggtitle(pheno, subtitle = paste0('(r = ',signif(cor(nbsnpsperchr[hsqsub$snpgroup],hsqsub$vgvp), 3),
                                       ', p-value = ',signif(cor.test(nbsnpsperchr[hsqsub$snpgroup],hsqsub$vgvp)$p.value, 3), ')')) + 
      theme_bw( ) +
      theme( axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5, hjust = 1), 
             axis.text.y = element_text(size = 10),
             axis.title.y = element_text(size = 15),
             legend.position = 'none',
             axis.text.x.top = element_text(size = 15),
             legend.text = element_text(size = 15),
             strip.text.x = element_text(size=20))
    
    return(gg)
  }
  
  names(plotperchr) <- unique(res$Phenotype)
  
  pdf(file.path(output_dir, 'hsq-perchr.pdf'), 16, 14)
  print(plot_grid(plotlist = plotperchr[grep('left|right', names(plotperchr), invert = T)], align = "v",  ncol = 3, nrow = 4))
  dev.off()
  
  fwrite(res[,-1], file = file.path(output_dir, 'hsq-perchr.txt'), sep = '\t')

}
