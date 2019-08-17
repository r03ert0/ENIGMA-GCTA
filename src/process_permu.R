
suppressMessages(library(data.table))
suppressMessages(library(foreach))
library(ggplot2)


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

process_permu <- function(path_permu, prefix_permu, path_res, prefix_res, output_file = NULL, make_plot = FALSE) {
  
  #### ==== Read original GCTA results  ==== ###
  res <- readHsq(path_res = path_res, prefix_res = prefix_res)
  
  #### ==== Read and process permutation results ==== ###
  fpermu <- list.files(path_permu, pattern = prefix_permu, full.names = T)
  
  if (length(fpermu) == 0)
    stop(sprintf('Error: pattern %s not found in %s\nNo permutation result file has been found.', prefix_permu, path_permu))
  
  # read only files which are not empty
  respermu <- lapply(fpermu, function(f) {if (file.info(f)$size > 1) data.table::fread(f) else NA})
  respermu <- respermu[!is.na(respermu)]
  respermu <- data.table::rbindlist(respermu)
  # number of groups
  ng <- length(grep(pattern = '^V\\(G.)/Vp_se', x = names(respermu)))
  # number of snps in each group
  nbsnps <- unlist(respermu[1,paste0('n',1:ng), with = F])
  print(paste('Number of SNPs per group:',paste(nbsnps, collapse = ',')))
  # phenotypes
  phenos <- unique(res$pheno)
  
  indtot_var <- match("Sum of V(G)/Vp", names(respermu))
  indtot_se <- match("Sum of V(G)/Vp_se", names(respermu))
  
  res_perm <- sapply(phenos, function(phenosel) {
    
    print(phenosel)
    respheno <- respermu[pheno == phenosel,]
    
    if (nrow(respheno) < 2) {
      stop(sprintf("Error: not enough permutation result for phenotype %s.", phenosel))
    }
    
    # remove outliers for V(G1)/Vp_se
    V_se <- respheno[["V(G1)/Vp_se"]]
    permusel <- abs(V_se - median(V_se)) < 3*mad(V_se)
    
    hsq <- respheno[permusel, paste0('V(G', 1:ng, ')/Vp'), with=F]
    v <- respheno[permusel, paste0('V(G', 1:ng, ')'), with=F]
    
    f <- nbsnps/sum(nbsnps)
    score <- hsq - as.matrix(respheno[permusel, "Sum of V(G)/Vp"]) %*% f
    score_v <- score2 <- v - as.matrix(rowSums(v)) %*% f
    
    se <- apply(score, 2, sd)
    se_v <-apply(score_v, 2, sd)
    
    zscore <- t(t(score) / se)
    zscore_v <- t(t(score_v) / se_v)
    
    hsq_res <- res[pheno == phenosel,][,paste0('V(G', 1:ng, ')/Vp'), with=F]
    score_res <- hsq_res - res[pheno == phenosel,][["Sum of V(G)/Vp"]] * f
    zscore_res <- t(t(score_res) / se)
    pval <- 1 - abs(0.5 - rowMeans(t(zscore) >= c(zscore_res))) * 2
    
    list(hsq=hsq, se=se, zscore=zscore, hsq_res=hsq_res, zscore_res=zscore_res, pval=pval)
    
  }, simplify=FALSE, USE.NAMES=T)
  

  if (make_plot) {
    pdf('densityVgVp_permu_ukb_genic_nongenic.pdf', 10, 10)
    par(mfrow = c(4,4))
    for (phenosel in phenos) {
      for (numgroup in 1:ng) {
        pheno_var <- paste0('V(G',numgroup,')/Vp')
        
        plot(density(res_perm[[phenosel]]$hsq[[numgroup]]),
             main = paste(phenosel, pheno_var),
             xlab = pheno_var)
        abline(v = res[pheno == phenosel,][[pheno_var]], col = 'blue')
      
      }
    }
    dev.off()
  
    pdf('densityZscores_permu_ukb_genic_nongenic.pdf', 10, 10)
    par(mfrow = c(3,3))
    for (phenosel in phenos) {
      for (numgroup in 1:ng) {
        phenovar <- paste0('V(G',numgroup,')/Vp')
        
        zs <- res_perm[[phenosel]]$zscore[,numgroup]
        hist(zs, main = paste(phenosel, phenovar),
             breaks = 50,
             xlab = 'z-scores')
        
        zres <- res_perm[[phenosel]]$zscore_res[numgroup]
        abline(v = zres, col = 'blue')
      }
    }
    dev.off()
  }
  
  pvals <- sapply(res_perm, '[[', 'pval')
  pvals <- as.data.frame(rbind(pvals, unlist(lapply(res_perm, function(x) nrow(x$zscore)))))
  pvals$nbSNPs <-  as.numeric(c(nbsnps,NA))
  rownames(pvals) <- c(rownames(pvals)[-3], 'nbPermu')
    
  if (!is.null(output_file))
    fwrite(pvals, file = output_file, sep = '\t', row.names = TRUE)

  return(pvals)
}
  
# process_permu(path_permu, prefix_permu, path_res, prefix_res, output_file, make_plot = FALSE)

  