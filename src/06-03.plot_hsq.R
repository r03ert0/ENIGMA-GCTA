library(ggplot2)
suppressMessages(library(data.table))
suppressMessages(library(foreach))
library(corrplot)
library(ggrepel)

## hsqsummary: data.table containing the hsq summary output 
##             or a string giving the path of the file containing the hsq output summary

## TO DO: voir normalisation pour per chr
plot_hsq_univariate <- function(hsqsummary, outpdf = NULL, normByPropSNPs = FALSE, 
                                nbSNPs = NULL, pdfwidth = 15, pdfheight = 8, 
                                ignore_snpgroup = NULL) {
  
  if (!is.null(outpdf))
    pdf(outpdf, width = pdfwidth, height = pdfheight)
  
  
  #directory       snpgroup        Phenotype       VG/VP ± s.e.    LRT     df      P-value N
  if (is.character(hsqsummary)) {
    hsq <- fread(hsqsummary, header = T, sep = '\t')
  } else {
    hsq <- hsqsummary
    if (!(all(c('snpgroup', 'Phenotype', 'VG/VP ± s.e.') %in% colnames(hsq))))
        message('the hsq summary data.table must contain the colnames snpgroup, Phenotype, and VG/VP ± s.e.')
  }
  
  hsq$vgvp <- as.numeric(gsub(' .*', '',hsq$`VG/VP ± s.e.`))
  hsq$vgvp_se <- as.numeric(gsub('.* ', '',hsq$`VG/VP ± s.e.`))
  
  if (!is.null(ignore_snpgroup)){
      hsq <- hsq[!(snpgroup %in% ignore_snpgroup),]
  }
  
  if (!is.null(nbSNPs)) {
    # extract SNP groups used to partition variance and get nb of SNPs
    snpgroups <- hsq$snpgroup
    snpgroups[snpgroups %in% c('allchr','nopca')] <- 'all'
    nbsnpbygroup <- nbSNPs[snpgroups]
    names(nbsnpbygroup) <- snpgroups
    hsq$nbsnp <- nbsnpbygroup

    hsq$snpgroup <- as.character(hsq$snpgroup)
    
    if (normByPropSNPs) {
      # extract VG/VP coputed using allchr to normalize per-chr VG/VP
      vgvptotal <- hsq[snpgroup == 'allchr',]$vgvp
      names(vgvptotal) <- hsq[snpgroup == 'allchr',]$Phenotype
      hsq <- hsq[snpgroup != 'allchr',]
      hsq$vgvp <- hsq$vgvp/vgvptotal[hsq$Phenotype]/(hsq$nbsnp/nbSNPs['all'])
    }
    
    # add nb SNPs to group name (such that it will show in the barplot legend with ggplot2)
    hsq$snpgroup <- paste0(hsq$snpgroup, ' (', hsq$nbsnp, ' SNPs)')
  }  else {
    if (normByPropSNPs)
      message('The number of SNPs are missing and the h^2 estimates won\'t be divided')
  }


  ## for each directory, plot barplots for each subset of SNPs tested
  foreach (dir = unique(hsq$directory)) %do% {
    
    hsq_sub <- hsq[directory == dir, ]
    hsq_sub <- hsq_sub[snpgroup != 'allchr',]
    
    gg <- ggplot(hsq_sub, aes(x=Phenotype, y=vgvp, fill=Phenotype)) + 
      geom_bar(position=position_dodge(), stat="identity",
               colour="black", # Use black outlines,
               size=.3) +
       ylab(if (!normByPropSNPs) 'VG/VP  ± s.e.' else '%VG/VP divided by %SNPs') +
      xlab('Phenotype') + 
      scale_fill_hue(name="Supplement type", # Legend label, use darker colors
                     breaks=c("OJ", "VC"),
                     labels=c("Orange juice", "Ascorbic acid")) +
      facet_wrap(~snpgroup)+ 
      theme_bw( ) +
      theme( axis.text.x = element_text(angle = 90, size = 15, vjust = 0.5, hjust = 1), 
             axis.text.y = element_text(size = 10),
             axis.title.y = element_text(size = 15),
             axis.text.x.top = element_text(size = 15),
             strip.text.x = element_text(size=15), 
             legend.text = element_text(size = 15))
    
    if (!normByPropSNPs) {
      gg <- gg+geom_errorbar(aes(ymin=vgvp-vgvp_se, ymax=vgvp+vgvp_se),
                             size=.3,    # Thinner lines
                             width=0,
                             position=position_dodge()) 
    }
    
    print(gg)
    
    
    return(NULL)   
  }
  
  if (!is.null(outpdf))
    dev.off()
  
  return(hsq)
  
}


plot_hsq_univariate_partition <- function(hsqsummary, normByPropSNPs = FALSE, nbSNPs = NULL, 
                                          outpdf = NULL, pdfwidth = 15, pdfheight = 8,
                                          ignore_snpgroup = NULL) {
  
  if (!is.null(outpdf))
    pdf(outpdf, width = pdfwidth, height = pdfheight)
  
  #directory       snpgroup        Phenotype       VG/VP ± s.e.    LRT     df      P-value N
  if (is.character(hsqsummary)) {
    hsq <- fread(hsqsummary, header = T, sep = '\t')
  } else {
    hsq <- hsqsummary
    if (!(all(c('snpgroup', 'Phenotype', 'VG/VP ± s.e.') %in% colnames(hsq))))
      message('the hsq summary data.table must contain the colnames snpgroup, Phenotype, and VG/VP ± s.e.')
  }
  
  if (nrow(hsq) == 0) {
    stop(paste("Error: void table in", hsqsummary))
  }
  
  hsq <- hsq[, setdiff(names(hsq),'nb_samples'), with = F]
  hsq <- melt(hsq, id.vars = names(hsq)[1:3])
  
  
  if (!is.null(ignore_snpgroup)){
    hsq <- hsq[!(snpgroup %in% ignore_snpgroup),]
  }
  
  ## Work-around for now to take into account the number of SNPs when hsq for genic-marginxx: 
  # if snpgroup contains genic-marginxx, the margin information is added to the variable id for each line
  # by replacing 'genic' of the variable column by 'genic-marginxx' in the snpgroup column
  if (length(grep('genic-margin', hsq$snpgroup)) == nrow(hsq)) {
    hsq[, variable := as.character(variable)]
    hsq[, variable := gsub('genic',gsub('\\.1$|\\.2$|\\.3$||\\.4$','',snpgroup),variable), by = seq_len(nrow(hsq))]
  } else {
    if (length(grep('updown-margin.*-.*', hsq$snpgroup)) == nrow(hsq)) {
      hsq[, variable := as.character(variable)]
      hsq[grep('-bin2', variable), variable := gsub('-bin2', '', gsub('updown-margin',gsub('\\.1$|\\.2$|\\.3$||\\.4$','',snpgroup),variable)), by = seq_len(nrow(hsq[grep('-bin2', variable),]))]
      hsq[grep('-bin1', variable), variable := gsub('-..-bin1', '', gsub('updown-margin',gsub('\\.1$|\\.2$|\\.3$||\\.4$','',snpgroup),variable)), by = seq_len(nrow(hsq[grep('-bin1', variable),]))]
      hsq[grep('nongenic', variable), 
          variable := gsub('nongenic', paste0('nongenic-margin',gsub('.*-', '', gsub('\\.1$|\\.2$|\\.3$||\\.4$','',snpgroup))), variable), 
          by = seq_len(nrow(hsq[grep('nongenic', variable),]))]
    } else {
      if (length(grep('updown-margin', hsq$snpgroup)) == nrow(hsq)) {
          hsq[, variable := as.character(variable)]
          hsq[, variable := gsub('updown-margin',gsub('\\.1$|\\.2$|\\.3$||\\.4$','',snpgroup),variable), by = seq_len(nrow(hsq))]
          hsq[grep('nongenic', variable), 
              variable := gsub('nongenic', paste0('nongenic-margin',gsub('.*margin', '', gsub('\\.1$|\\.2$|\\.3$||\\.4$','',snpgroup))), variable), 
              by = seq_len(nrow(hsq[grep('nongenic', variable),]))]
          
      }
    }
  }
  
  
  
  hsq$vgvp <- as.numeric(gsub(' .*', '',hsq$value))
  hsq$vgvp_se <- as.numeric(gsub('.* ', '',hsq$value))
  
  if (!is.null(nbSNPs)) {
    # extract SNP groups used to partition variance and get nb of SNPs
    snpgroups <- gsub('^V|/VP ± s.e.','',hsq$variable)
    nbsnpbygroup <- nbSNPs[snpgroups]
    names(nbsnpbygroup) <- snpgroups
    nbsnpbygroup[names(nbsnpbygroup) == "Sum_of_V(G)/Vp"] <- sum(nbsnpbygroup[unique(snpgroups)], na.rm = TRUE)
    hsq$nbsnp <- nbsnpbygroup
    hsq$variable <- as.character(hsq$variable)
    # extract VG/VP results for normalization and remove it from plots
    vgvptotal <- hsq[variable == 'Sum_of_V(G)/Vp',]$vgvp
    names(vgvptotal) <- hsq[variable == 'Sum_of_V(G)/Vp',]$Phenotype
    vgvptotal_se <- hsq[variable == 'Sum_of_V(G)/Vp',]$vgvp_se
    names(vgvptotal_se) <- hsq[variable == 'Sum_of_V(G)/Vp',]$Phenotype
    hsq <- hsq[variable != 'Sum_of_V(G)/Vp',]
    #nbsnpbygroup <- nbsnpbygroup[names(nbsnpbygroup) != 'Sum_of_V(G)/Vp']

    if (normByPropSNPs) {
      if (!is.null(nbSNPs)) {
        hsq$vgvp <- (hsq$vgvp/vgvptotal[hsq$Phenotype])/(hsq$nbsnp/nbsnpbygroup["Sum_of_V(G)/Vp"])
        #hsq$vgvp_se <- rep(0,length(hsq$vgvp_se))  #need to use permutations to estimate SE(hsq$vgvp_se/vgvptotal_se[hsq$Phenotype])/(nbsnpbygroup/nbsnpbygroup["Sum_of_V(G)/Vp"])
      } else {
        message('The number of SNPs are missing and the h^2 estimates won\'t be normalized')
      }
        }
    
    # add nb SNPs to group name (such that it will show in the barplot legend with ggplot2)
    hsq$variable <- paste0(hsq$variable, ' (', hsq$nbsnp, ' SNPs)')
  }
    
  
  ## for now only .1 bc no p-values extracted from results, see summary_hsq.sh
  hsq <- hsq[grep ("\\.1", hsq$snpgroup),]
  
  ## for each directory, plot barplots for each subset of SNPs tested
  foreach (dir = unique(hsq$directory)) %do% {
    
    hsq_sub <- hsq[directory == dir, ]
    
    gg <- ggplot(hsq_sub, aes(x=Phenotype, y=vgvp, fill=variable)) + 
      geom_bar(position=position_dodge(1), stat="identity",
               colour="black", # Use black outlines,
               size=.3) +
      ylab(if (!normByPropSNPs) 'VG/VP  ± s.e.' else '%VG/VP  divided by %SNPs') +
      xlab('Phenotype') + 
      #scale_fill_hue(name="Supplement type", # Legend label, use darker colors
      #               breaks=c("OJ", "VC"),
      #               labels=c("Orange juice", "Ascorbic acid")) +
      facet_wrap(~snpgroup, nrow = 3, scales = 'free_y')+ 
      theme_bw( ) +
      theme( axis.text.x = element_text(angle = 90, size = 15, vjust = 0.5, hjust = 1), 
             axis.title.y = element_text(size = 15),
             axis.text.y = element_text(size = 15),
             axis.text.x.top = element_text(size = 15),
             strip.text.x = element_text(size=20), 
             legend.text = element_text(size = 15))
    
    
    if (!normByPropSNPs) {
      gg <- gg+geom_errorbar(aes(ymin=vgvp-vgvp_se, ymax=vgvp+vgvp_se),
                  size=.3,    # Thinner lines
                  width=0,#.2,
                  position=position_dodge(1)) 
    }
      
    print(gg)
    
    return(NULL)   
  }
  
  if (!is.null(outpdf))
    dev.off()
  
  return(hsq)
  
}




plot_hsq_bivariate <- function(hsqsummary, outpdf = NULL, pdfwidth = 15, pdfheight = 10, ignore_snpgroup = NULL) {
  
  if (!is.null(outpdf))
    pdf(outpdf, width = pdfwidth, height = pdfheight)
  
  if (is.character(hsqsummary)) {
    hsq <- fread(hsqsummary, header = T, sep = '\t')
  } else {
    hsq <- hsqsummary
    if (!(all(c('snpgroup', 'phenotype1', 'phenotype2', 
                'VG/VP_phenotype1 ± s.e.', 'VG/VP_phenotype2 ± s.e.', 
                'rG ± s.e.', 'pvalue') %in% colnames(hsq))))
      message('the hsq summary data.table must contain the colnames snpgroup, phenotype1, phenotype2, VG/VP_phenotype1 ± s.e., VG/VP_phenotype2 ± s.e., 
              and rG ± s.e., and pvalue.')
  }
  
  if (!is.null(ignore_snpgroup)){
    hsq <- hsq[!(snpgroup %in% ignore_snpgroup),]
  }
  
  hsq$rg <- as.numeric(gsub(' .*', '',hsq$`rG ± s.e.`))
  hsq$rg_se <- as.numeric(gsub('.* ', '',hsq$`rG ± s.e.`))
  hsq[, rg_sebis := rg_se]
  hsq[rg < 0, rg_sebis := -rg_se]
  hsq$pvalue <- as.numeric(gsub('.* ', '',hsq$pvalue))
  hsq[, snpgroup:= gsub("\\\\",'',snpgroup)]
  hsq[, pvaluechar:= '']
  hsq[pvalue < 0.05, pvaluechar:= '*']
  hsq[pvalue < 0.01, pvaluechar:= '**']
  hsq[pvalue < 0.01, pvaluechar:= '***']
  hsq[, phenotype1 := gsub("all.rg=1.|all.rg=0.", '', phenotype1)]
  hsq[, phenotype2 := gsub("all.rg=1.|all.rg=0.", '', phenotype2)]
  hsq[, phenotype1 := gsub("rg=1.|rg=0.", '', phenotype1)]
  hsq[, phenotype2 := gsub("rg=1.|rg=0.", '', phenotype2)]
  
  ## for now only .1 bc no p-values extracted from results, see summary_hsq.sh
  
  hsq[, snpgroup := gsub("\\", '', snpgroup, fixed = T),]
  hsq[, pheno_pair := paste(sort(c(phenotype1, phenotype2)), collapse = '_'), by = seq_len(nrow(hsq))]
  
  # ## assigning NA to aberrant values
  # hsq[rg < -1 | rg > 1, rg_se := NA ]
  # hsq[rg < -1 | rg > 1, rg_sebis := NA ]
  # hsq[rg < -1 | rg > 1, rg := NA ]
  # hsq[rg < -1 | rg > 1, pvalue := NA ]
  # hsq[rg < -1 | rg > 1, pvaluechar := NA ]
  
  ## 1. Plot rG  ± s.e. as barplots
  gg <- ggplot(hsq[which(!duplicated(pheno_pair)),], 
               aes(x=pheno_pair, y=rg, fill=pheno_pair)) + 
    geom_bar(position=position_dodge(), stat="identity",
             colour="black", # Use black outlines,
             size=.3) + ggplot2::coord_cartesian(ylim = c(-1.2,1.2)) +
    geom_errorbar(aes(ymin=rg-rg_se, ymax=rg+rg_se),
                  size=.3,    # Thinner lines
                  width=.2,
                  position=position_dodge(.9)) +
    ylab('rG  ± s.e.') +
    xlab('Phenotype') + 
    facet_wrap(~snpgroup)+ 
    theme_bw( ) +
    theme( axis.text.x = element_text(angle = 90, size = 15, vjust = 0.5, hjust = 1), 
           axis.text.y = element_text(size = 15),
           axis.title.y = element_text(size = 15),
           legend.position = 'none',
           axis.text.x.top = element_text(size = 15),
           legend.text = element_text(size = 15),
           strip.text.x = element_text(size=15))
  
  
  print(gg)

  ## 2. Plot rG and indicate p-values for allrg=1 and allrg=0 as barplots
  for (snpgroupsel in c('all.rg=0','all.rg=1')) {
    hsqsub <- hsq[snpgroup %in% snpgroupsel,][which(!duplicated(pheno_pair)),]#[grep('ICV', pheno_pair, invert = TRUE)]
    hsqsub[, vjust := -0.5]
    hsqsub[rg_sebis < 0, vjust := 1]
    if (nrow(hsqsub) == 0) {
      warning(sprintf("No bivariate value found for group %s, skipping plot.", snpgroupsel))
      next
    }
    gg <- ggplot(hsqsub, 
                 aes(x=pheno_pair, y=rg, fill=pheno_pair)) + 
      geom_bar(position=position_dodge(), stat="identity",
               colour="black", # Use black outlines,
               size=.3) + ggplot2::coord_cartesian(ylim = c(-1.2,1.2)) +
      geom_text(aes(x=pheno_pair, y=rg+rg_sebis,  label = pvaluechar, vjust = vjust),
        position = position_dodge(width = 1))+
       # vjust = hsqsub$vjust) +
      geom_errorbar(aes(ymin=rg-rg_se, ymax=rg+rg_se),
                    size=.3,    # Thinner lines
                    width=.2,
                    position=position_dodge(.9)) +
      ylab('rG  ± s.e.') +
      xlab('Phenotype') + 
      facet_wrap(~snpgroup)+ 
      theme_bw( ) +
      theme( axis.text.x = element_text(angle = 90, size = 15, vjust = 0.5, hjust = 1), 
             axis.text.y = element_text(size = 15),
             axis.title.y = element_text(size = 15),
             legend.position = 'none',
             axis.text.x.top = element_text(size = 15),
             legend.text = element_text(size = 15),
             strip.text.x = element_text(size=15))
    
       print(gg)
  }
  
  ## 2bis. Plot pvalue for allrg=1 and allrg=0 as barplots
  for (snpgroupsel in c('all.rg=0','all.rg=1')) {
      hsqsub <- hsq[snpgroup %in% snpgroupsel,][which(!duplicated(pheno_pair)),]
      if (nrow(hsqsub) == 0) {
        warning(sprintf("No bivariate value found for group %s, skipping plot.", snpgroupsel))
        next
      }
      gg <- ggplot(hsqsub, 
                 aes(x=pheno_pair, y=-log10(pvalue), fill=pheno_pair)) + 
      geom_bar(position=position_dodge(), stat="identity",
               colour="black", # Use black outlines,
               size=.3) + #ggplot2::coord_cartesian(ylim = c(-1.5,1.5)) +
      ylab('-log10(p-value)') +
      xlab('Phenotype') + geom_abline(slope = 0,intercept = -log10(0.05), lty = 3) + 
      facet_wrap(~snpgroup)+ 
      theme_bw( ) +
      theme( axis.text.x = element_text(angle = 90, size = 15, vjust = 0.5, hjust = 1), 
             axis.text.y = element_text(size = 15),
             axis.title.y = element_text(size = 15),
             legend.position = 'none',
             axis.text.x.top = element_text(size = 15),
             legend.text = element_text(size = 15),
             strip.text.x = element_text(size=15))
    
    
      print(gg)
  }
  
  ## 3. Plot rG  as corrplot (find a way to add ± s.e.)
  # hsqsub = hsq[snpgroup == 'all',]
  hsqsub = hsq
  avpheno <- unique(c(hsqsub$phenotype1,hsqsub$phenotype2))
  npheno <- length(avpheno)
  # correlation matrix
  mcor = matrix(NA, nrow = npheno, ncol = npheno, dimnames = list(avpheno, avpheno))
  for (i in 1:nrow(hsqsub)) {
    mcor[hsqsub$phenotype1[i], hsqsub$phenotype2[i]] <- 
      mcor[hsqsub$phenotype2[i], hsqsub$phenotype1[i]] <- hsqsub$rg[i]
    
  }
  #mcor <- mcor[grep('ICV',rownames(mcor), invert = T),grep('ICV',colnames(mcor), invert = T)]
  diag(mcor) <- 1
  # standard error matrix
  mse = matrix(NA, nrow = npheno, ncol = npheno, dimnames = list(avpheno, avpheno))
  for (i in 1:nrow(hsqsub)) {
    mse[hsqsub$phenotype1[i], hsqsub$phenotype2[i]] <- 
      mse[hsqsub$phenotype2[i], hsqsub$phenotype1[i]] <- hsqsub$rg_se[i]
    
  }
  #mse <- mse[grep('ICV',rownames(mse), invert = T),grep('ICV',colnames(mse), invert = T)]
  diag(mse) <- 0
  #mcor <- mcor[which(rowSums(mse, na.rm = T) < 10000),]
  #mse <- mse[which(rowSums(mse, na.rm = T) < 10000),]
  mcor <- mcor[,rownames(mcor)]
  mse <- mse[,rownames(mse)]
  
  ## assigning NA to aberrant values
  mcor[mcor > 1] <- NA
  mcor[mcor < -1] <- NA
  mse[mse+mcor > 1] <- NA
  mse[mse-mcor < -1] <- NA
  
  # plot both
  # plot corr values 
  # corrplot.mixed(mcor,order = "AOE", lower.col = "black", 
  #                upper.col =  col3(200), cl.lim = c(-1, 1), 
  #                na.label = 'square')
  corrplot(mcor, type = "upper", method = 'circle', diag = F, 
           tl.pos = 'td', tl.cex = 0.8, 
           #plotCI = plotCI, #col = upper.col, 
           addCoef.col = "black", number.cex = .7,  #col = upper.col, 
           lowCI.mat = mcor-mse, uppCI.mat = mcor+mse, 
           mar = c(5,5,5,5)) 
  
  # plot corr values and confidence intervals
  tl.pos = c("d", "lt", "n")[1]
  diag = c("n", "l", "u")[1]
  plotCI = 'circle' #c("n", "square",   "circle", "rect")
  col3 <- colorRampPalette(c("red", "white", "blue"))
  corrplot(mcor, type = "lower", method = 'color', diag = TRUE, order = 'alphabet',
           tl.pos = 'n', tl.cex = 0.8, col =  col3(200),
           #plotCI = plotCI, #col = upper.col,
           addCoef.col = "black", number.cex = .8,  #col = upper.col,
           mar = c(0,0,0,0), lowCI.mat = mcor-mse, uppCI.mat = mcor+mse)

  if (sum((mcor+mse) > 1 | (mcor-mse) < -1, na.rm = T) > ncol(mcor)) { 
    corrplot(mcor, type = "upper", method = 'circle', diag = T, add = TRUE,order = 'alphabet',
             col =  col3(200),cl.offset = -1, cl.pos = 'n',
             tl.pos = 'd', tl.cex = 0.8, plotCI = plotCI, #col = upper.col,
             mar = c(0,0,0,0), lowCI.mat = mcor-mse, uppCI.mat = mcor+mse)
  }
  
 # symbols(1:npheno, npheno:1, add = TRUE, bg = 'white', fg = 'grey', inches = FALSE, squares = rep(1, npheno))

  ## 4 Plot rG as corrplot for gcta parameters rg=0 and rg=1
  for (snpgroupsel in c('all.rg=0','all.rg=1')) {
  hsqsub = hsq[snpgroup == snpgroupsel,]
  if (nrow(hsqsub) == 0) {
    warning(sprintf("No bivariate value found for group %s, skipping plot.", snpgroupsel))
    next
  }
  avpheno <- unique(c(hsqsub$phenotype1,hsqsub$phenotype2))
  npheno <- length(avpheno)
  mcor = matrix(NA, nrow = npheno, ncol = npheno, dimnames = list(avpheno, avpheno))
  for (i in 1:nrow(hsqsub)) {
    mcor[hsqsub$phenotype1[i], hsqsub$phenotype2[i]] <- 
      mcor[hsqsub$phenotype2[i], hsqsub$phenotype1[i]] <- hsqsub$rg[i]
    
  }
  #mcor <- mcor[grep('ICV',rownames(mcor), invert = T),grep('ICV',colnames(mcor), invert = T)]
  diag(mcor) <- 1
  mpval = matrix(NA, nrow = npheno, ncol = npheno, dimnames = list(avpheno, avpheno))
  for (i in 1:nrow(hsqsub)) {
    mpval[hsqsub$phenotype1[i], hsqsub$phenotype2[i]] <- 
      mpval[hsqsub$phenotype2[i], hsqsub$phenotype1[i]] <- hsqsub$pvalue[i]
    
  }
  #mpval <- mpval[grep('ICV',rownames(mpval), invert = T),grep('ICV',colnames(mpval), invert = T)]
  diag(mpval) <- 0
  
  mse = matrix(NA, nrow = npheno, ncol = npheno, dimnames = list(avpheno, avpheno))
  for (i in 1:nrow(hsqsub)) {
    mse[hsqsub$phenotype1[i], hsqsub$phenotype2[i]] <- 
      mse[hsqsub$phenotype2[i], hsqsub$phenotype1[i]] <- hsqsub$rg_se[i]
    
  }
  #mse <- mse[grep('ICV',rownames(mse), invert = T),grep('ICV',colnames(mse), invert = T)]
  diag(mse) <- 0
  
  ## assigning NA to aberrant values
  mcor[mcor > 1] <- NA
  mcor[mcor < -1] <- NA
  mse[mse+mcor > 1] <- NA
  mse[mse-mcor < -1] <- NA
  
  if (sum((mcor+mse) > 1 | (mcor-mse) < -1, na.rm = T) > ncol(mcor)) { 
    
    col3 <- colorRampPalette(c("red", "white", "blue")) 
    
    
    # plot corr values, confidence intervals with p-value information
    cp1 <- corrplot(mcor, col =  col3(200),order = 'alphabet',
                    type = "lower", method = 'shade',#'number',
                    diag = TRUE, tl.pos = "n", #cl.pos = "n",
                    mar = c(0,0,4,0),
                    addCoef.col = "black", number.cex = .8, title = snpgroupsel)
    # insig = 'label_sig', pch = 1, pch.cex = 1, sig.level = c(.001, .01, .05))#col = lower.col,
    cp2 <- corrplot(mcor, p.mat = mpval, order = 'alphabet',
                    type = "upper", method = 'circle', diag = T, col =  col3(200),
                    tl.pos = 'd', tl.cex = 0.8, plotCI = 'circle',  add = TRUE, cl.pos = 'n',
                    #addCoef.col = "black", number.cex = .7, #col = upper.col,
                    lowCI.mat = mcor-mse, uppCI.mat = mcor+mse, number.cex = .8,
                    insig = 'label_sig', pch = 1, pch.cex = 1, sig.level = c(.001, .01, .05))
    
  }
  # charsign <- apply(apply(mpval[rownames(cp1),colnames(cp1)],2, as.character), 2, 
  #                   function(x) { 
  #                     xx <- as.numeric(x)
  #                     x <- ''
  #                     x[xx < 0.001] <-  '***';  
  #                     x[xx < 0.01 & xx > 0.001] <-  '**';  
  #                     x[xx < 0.05 & xx > 0.01] <- '*' 
  #                     x
  #                   })
  # 
  # for (i in nrow(cp2):1) {
  #   for (j in 1:ncol(cp2)) {
  #      text(x=i, y = j, 
  #          labels = charsign[nrow(cp2)-i+1,j], #col = pch.col, cex = pch.cex, 
  #          lwd = 2)
  #   }
  # }
  }
  
  if (!is.null(outpdf))
    dev.off()
  
  return(hsq)
}

plotNbSNPs <- function(hsq_outdir) {
  
  #hsq_outdir <- "/Volumes/cinq/rto/ukb-annex/data/derived/ukb-9891/06.hsq/"
  nbsnpsfiles <- list.files(path = hsq_outdir, 
                       pattern = '.nbSNPs.txt',
                       full.names = TRUE,
                       recursive = TRUE)
  
  nbsnpsfiles <-  nbsnpsfiles[grep('hsq-summary/', nbsnpsfiles, invert = T)]
             
  names(nbsnpsfiles) <- gsub('.nbSNPs.txt','',basename(nbsnpsfiles))
  
  nbSNPs <- do.call(rbind,lapply(nbsnpsfiles, fread))
  setnames(nbSNPs, c('snp_group', 'nbsnps'))
  nbSNPsall <- nbSNPs
  
  nbSNPsbychr <- nbSNPs[grep('chr', snp_group),]
  nbSNPsbychr <- nbSNPsbychr[,snp_group := factor(snp_group, levels = paste0('chr',1:22))]
  nbSNPs <- nbSNPs[grep('chr', snp_group, invert = TRUE),]
  
  nbSNPs[, snp_group1 := sub('maf.*','maf',sub('-.*', '', snp_group))]
  
  #nbSNPs[, snp_group1 := factor(snp_group1, levels = c("all", "cnsexpression", "noncnsexpression", "genic", "nongenic", "maf", "neurodev", "nonneurodev"))]
  
  gg <- ggplot(nbSNPs, aes(x = as.factor(snp_group), y = as.numeric(nbsnps))) + ylab('#SNPs') + 
    geom_bar(stat = 'identity', position = 'dodge', colour="black",  size=.3) + #, fill = 'grey90') 
    #facet_wrap(~snp_group1, scales = 'free_x', ncol = 8) +
    xlab('SNP group') +
    scale_y_continuous(breaks = round(seq(0, max(nbSNPs$nbsnps), by = 10^4),1)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 15), strip.text.x = element_text(size=14))  #, panel.spacing.y = unit(4, "lines")) #
  print(gg)
 
  gg <- ggplot(nbSNPsbychr, aes(x = as.factor(snp_group), y = as.numeric(nbsnps))) + ylab('#SNPs') + 
    geom_bar(stat = 'identity', position = 'dodge', colour="black",  size=.3) + #, fill = 'grey90') 
    #facet_wrap(~snp_group1, scales = 'free_x', ncol = 8) +
    xlab('SNP group') +
    scale_y_continuous(breaks = round(seq(0, max(nbSNPs$nbsnps), by = 10^3),1)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 15), strip.text.x = element_text(size=14))  #, panel.spacing.y = unit(4, "lines")) #
  print(gg)
  
  return(nbSNPsall )
}

plotHSQ <- function(prefix_hsq_summary_file, hsq_outdir) {
  
  pdf(file = paste0(prefix_hsq_summary_file, 'nbSNPs.pdf'), width = 8, height = 6)
  nbsnps <- plotNbSNPs(hsq_outdir)
  dev.off()
  fwrite(nbsnps, file = paste0(prefix_hsq_summary_file, 'nbSNPs.txt'), sep = '\t', row.names = FALSE)
  
  
  nbsnpsperchr <- nbsnps[grep('chr',snp_group),]
  nbsnpsperchr <- nbsnpsperchr[,snp_group := factor(snp_group, levels = paste0('chr',1:22))]
  
  
  outfile_univariate <- c(all = paste0(prefix_hsq_summary_file,'hsq-all.txt'), 
                          nopca = paste0(prefix_hsq_summary_file,'hsq-nopca.txt'),
                          perchr = paste0(prefix_hsq_summary_file,'hsq-perchr.txt'))
  
  
  outfile_univariate_partition  <- c(cnsexpression = paste0(prefix_hsq_summary_file,'hsq-cnsexpression.txt'),
                                     maf = paste0(prefix_hsq_summary_file,'hsq-maf.txt'),
                                     neurodev = paste0(prefix_hsq_summary_file,'hsq-neurodev.txt'),
                                     genic = paste0(prefix_hsq_summary_file,'hsq-genic.txt'),
                                     genic_updown = paste0(prefix_hsq_summary_file,'hsq-genic_updown-margin.txt'),
                                     # genic_updown_intron1exon1 = paste0(prefix_hsq_summary_file,'hsq-genic_genic-exon1intron1.txt'),
                                     genic_updown_2bins = paste0(prefix_hsq_summary_file,'hsq-genic_updown-margin-2bins.txt'))
  
  
  outfile_bivariate <- paste0(prefix_hsq_summary_file,'hsq-biv.txt')
  
  ## ======== barplots   =========
  hsquniv <- mapply(hsqsummary = outfile_univariate, 
                    outpdf = gsub('.txt', '.pdf', outfile_univariate), 
                    plot_hsq_univariate, 
                    ignore_snpgroup = list(NULL,'allchr',NULL),
                    SIMPLIFY = FALSE)
  

  hsqunivpart <- mapply(hsqsummary = outfile_univariate_partition,
                            function(hsqsummary) {
                              plot_hsq_univariate_partition( 
                                hsqsummary = hsqsummary,
                                outpdf = gsub('.txt', '.pdf',  hsqsummary),
                                normByPropSNPs = FALSE, 
                                nbSNPs = structure(.Data = nbsnps$nbsnps, .Names = nbsnps$snp_group)
                              )
                            },SIMPLIFY = FALSE)

  
  ## ======== divide h^2 estimates by the proportion of SNPs per group  =========
  hsqunivperchrnorm <-  plot_hsq_univariate( 
    hsqsummary = outfile_univariate[3],
    outpdf = gsub('.txt', '_normByPropSNPs.pdf',  outfile_univariate[3]),
    normByPropSNPs = TRUE, 
    nbSNPs = structure(.Data = nbsnps$nbsnps, .Names = nbsnps$snp_group)
  )
  
  
  hsqunivpartnorm <- mapply(hsqsummary = outfile_univariate_partition,
                            function(hsqsummary) {
                              plot_hsq_univariate_partition( 
                                hsqsummary = hsqsummary,
                                outpdf = gsub('.txt', '_normByPropSNPs.pdf',  hsqsummary),
                                normByPropSNPs = TRUE, 
                                nbSNPs = structure(.Data = nbsnps$nbsnps, .Names = nbsnps$snp_group)
                              )
                            },SIMPLIFY = FALSE)
  

  ## ======== hsq per chr vs #SNPs per chr  =========
  lm_eqn <- function(df){ #source http://goo.gl/K4yh
    m <- lm(y ~ x, df);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                     list(a = format(coef(m)[1], digits = 2), 
                          b = format(coef(m)[2], digits = 2), 
                          r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));                 
  }
  
  
  pdf(file = paste0(prefix_hsq_summary_file, 'hsq-perchr_vs_nbSNPs.pdf'), width = 8, height = 6)
  for (pheno in unique(hsquniv$perchr$Phenotype)) {
    
    print(pheno)
    hsqsub <- hsquniv$perchr[Phenotype == pheno & snpgroup != 'allchr',]
    gg <- 
      ggplot(hsqsub, 
             aes(x = nbsnpsperchr[match(as.character(hsqsub$snpgroup),snp_group),]$nbsnps,
                 y = vgvp)) + #, colour = snpgroup)) + 
      geom_point(size = 2) + geom_smooth(method=lm) + 
      geom_text_repel(aes(x = nbsnpsperchr[match(as.character(hsqsub$snpgroup),snp_group),]$nbsnps,
                     y = vgvp, 
                     label = snpgroup)) + 
      geom_errorbar(aes(ymin=vgvp-vgvp_se, ymax=vgvp+vgvp_se),
                    size=.3,    # Thinner lines
                    width=.2,
                    position=position_dodge(.9)) +
      scale_x_continuous(breaks = round(seq(0, max(nbsnps$nbsnps), by = 10^3),1)) +
      ylab('VG/VP  ± s.e.') + xlab('#SNPs') + 
      geom_text(x = 8000, y = max(hsqsub$vgvp)-0.05, 
                label = lm_eqn(data.frame(x = nbsnpsperchr[match(as.character(hsqsub$snpgroup),snp_group),]$nbsnps,y = hsqsub$vgvp)), 
                parse = TRUE, size = 4) +
      ggtitle(paste0(pheno, ', Pearson cor =',signif(cor(nbsnpsperchr[match(as.character(hsqsub$snpgroup),snp_group),]$nbsnps,hsqsub$vgvp), 3),
                     ', p-value =',signif(cor.test(nbsnpsperchr[match(as.character(hsqsub$snpgroup),snp_group),]$nbsnps,hsqsub$vgvp)$p.value, 3))) + 
      theme_bw( ) +
      theme( axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5, hjust = 1), 
             axis.text.y = element_text(size = 10),
             axis.title.y = element_text(size = 15),
             legend.position = 'none',
             axis.text.x.top = element_text(size = 15),
             legend.text = element_text(size = 15),
             strip.text.x = element_text(size=20))
    
    print(gg)
    
  }
  dev.off()
  
  plot_hsq_bivariate(outfile_bivariate, outpdf = gsub('.txt', '.pdf',  outfile_bivariate))
  
}


#hsq_outdir <- "/Volumes/cinq-2/rto/ukb-annex/data/derived/ukb-22110/06.hsq/"
#prefix_hsq_summary_file <- "/Volumes/cinq-2/rto/ukb-annex/data/derived/ukb-22110/06.hsq/hsq-summary/summary_"
#hsq_outdir <- "/Volumes/cinq/rto/ukb-annex/data/derived/IMAGEN123/before_longrangeLD_filtering/06.hsq/"
#prefix_hsq_summary_file <- "/Volumes/cinq/rto/ukb-annex/data/derived/IMAGEN123/before_longrangeLD_filtering/06.hsq/hsq-summary/summary_"
#hsq_outdir <- "/Volumes/cinq/rto/ukb-annex/data/derived/LOTHIAN/06.hsq/"
#prefix_hsq_summary_file <- "/Volumes/cinq/rto/ukb-annex/data/derived/LOTHIAN/06.hsq/hsq-summary/summary_"
#plotHSQ(prefix_hsq_summary_file = prefix_hsq_summary_file, hsq_outdir = hsq_outdir)


