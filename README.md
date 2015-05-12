# ENIGMA-GCTA

ENIGMA-GCTA is a project within the ENIGMA consortium to establish the role of frequent, strongly polygenic, genetic polymorphism on the determination of human neuroanatomical diversity.

This repository contains the bas scripts and C source code used to estimate the phenotypic variance captured by common genetic variation in the IMAGEN consortium, described in our article "Genomic architecture of human neuroanatomical diversity", published in Molecular Psychiatry, doi: 10.1038/mp.2014.99.

The scripts are within the directory 'bin', and the different lists of genes used for genomic partition are within the directory 'lists.

## Description of the files

The 'scripts' directory contains the following:

1. grm-sge.sh
This script launches the parallel computation (using SGE) of various genetic relationship matrices (GRMs) computed from different subsets of single nucleotide polymorphisms (SNPs). The script is divided in 5 sections:

1. Compute a GRM using all (filtered, R2 pruned) SNPs
2. Compute 1 GRM per chromosome
3. Compute GRMs from genic and nongenic SNPs
4. Make GRMs for low, medium and high MAF
5. Make GRMs for different functional gene partitions. Currently: genes preferentially expressed in CNS and genes involved in neurodevelopment.

The execution of each section is controlled by the boolean value in the if/then block. Set the value to true for execution, put no value to prevent execution.

This script relies on GCTA for the computation of GRMs, and on the C command genes2snps to compute list of SNPs from list of genes (see description bellow).

2. hsq-sge.sh
This script launches the estimations of phenotypic variance captured by the GRMs computed from the different genomic partitions. It is divided in 7 sections, whose execution can be controlled through the boolean value in their if/then blocks:

1. table 1: Using GRMs computed from all SNPs, including 10 PCs as covariates
2. table 2: Using GRMs computed from all SNPs, not including 10 PCs as covariates
3. table 3: Using GRMs computed from genic and nongenic SNP partitions
4. table 4: Using GRMs computed from genic SNPs preferentially expressed in the CNS, versus remaining genic SNPs, versus nongenic SNPs
5. table 5: Using GRMs computed from genic SNPs involved in neurodev, versus remaining genic SNPs, versus nongenic SNPs
6. table 6: Using partitions based on MAF
7. table 7: Using one GRM per chromosome

hsq-sge.sh relies on GCTA for the estimation of variance captured by SNPs, and on the C command 'mygcta' to simplify the syntax of statistical models (see description bellow).

* get-results.sh
This script read all the log and result files produced by hsq-sge.sh, collects relevant data, and produces tables with the results. It is divided in 7 sections, reflecting each of the sections of hsq-sge.sh.

* genes2snps/
genes2snps takes text files with lists of genes (like those inside the 'lists' directory) and produces lists of SNPs, to be used in the creation of GRMs. The directory contains the C code for the command, genes2snps.c, and a binary file compiled for Mac OS X. The compilation does not require any library, and can be performed with the command:

    gcc -Wall genes2snps.c -o genes2snps

* mygcta/
mygcta is a C wrapper to gcta which simplifies the definition of statistical models. gcta requires a single file with all quantitative covariates (using the --qcovar switch) and another single file with all categorial covariates (using the --covar switch). mygcta allow to have multiple --qcovar and --covar switches, each with a single covariate. It combines then all these individual covariates in a format suitable for gcta. The command does not require any library, and can be compiled using:

    gcc -Wall mygcta.c -o mygcta


## Description of the analysis

The scripts suppose that you have set up 3 directories: genotype (original genotypes, filtered for maf, hwe, etc., and pruned for R^2), lists (functional lists of cns genes, neurodev, and genic), and phenotype (in the case of our paper, these were neuroimaging phenotypes, plus age, centre, viq, piq, height).

In addition to the scripts in the bin directory, you need to have gcta and plink (we used versions gcta123 and plink-1.07). Calling grm-sge.sh will launch all the computations of GRM matrices, for the complete genome, chromosome by chromosome, and by functional sets. hsq-sge.sh will then compute the heritabilities for the different phenotypes. The analyses produce many log files. Running the script get-results.sh inside the results directory scans all the log files and makes pure text tables.

As we said previously, gcta only accepts one --covar switch and one --qcovar switch to add covariates. So, if you have several categorial or quantitative covariates you have to make a special text file with one covariate per column. I wrote mygcta because I wanted to have more flexibility to add covariates to the model. When you use mygcta you just add covars and qcovars to the command line and mygcta will make the special text files for you. Look at hsq-sge.sh to see how to use it.

hsq-sge.sh uses age, sex, centre and 10 first principal components as covariates. In our analyses, age does not have any effect but we used them nevertheless for completness.

The effect of the first 10 PCs is small, but it is important to add them to control for population structure.

On the contrary, the sex and scanning centre effects were very important, especially centre. In addition to methodological differences (different MRI machines, etc), there seem to be also differences in the populations recruited.

We used only autosomes in your analyses (using the --autosomes switch in gcta).