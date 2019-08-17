#!/usr/bin/env python3
"""Compute heritability for all SNPs"""

import os
import sys
import preprocessing
import config_dataset


def main(config_file):

    """Entry point if called as an executable"""

    config = config_dataset.config_dataset(config_file)
    # quantitative covariables
    # qcovar = [os.path.join(PHE_DIR, 'age.txt'),PCS]
    qcovar_par = []
    for qcov in config.quant_covar:
        qcovar_par.append('--qcovar')
        qcovar_par.append(qcov)

    # qualitative covariables
    # covar = [os.path.join(PHE_DIR, 'sex.txt'),os.path.join(PHE_DIR, 'centre.txt')]
    covar_par = []
    for cov in config.qual_covar:
        covar_par.append('--covar')
        covar_par.append(cov)

    var_par = qcovar_par + covar_par

    # ========= 1. All SNPS, 10 PCs =============

    in_file = os.path.join(config.grm_dir, 'grm-all', 'all')
    out_hsq = os.path.join(config.hsq_dir, 'hsq-all', 'all')

    for pheno in config.phe_list:

        out_file = out_hsq + '.' + pheno
        print('running h^2 gcta estimation for phenotype: ' + pheno)
        pheno = os.path.join(config.phe_dir, pheno+'.txt')
        pars = var_par + ['--pheno', pheno, str(config.reml_call)]

        preprocessing.gcta_hsq(in_file=in_file,
                               out_file=out_file,
                               gcta=config.gcta,
                               other_gcta_par=pars,
                               ncpus=config.nbproc,
                               mygcta=config.mygcta,
                               sbatch=config.use_sbatch,
                               sbatch_par_j="hsq-all")

    # extract number of SNPs per chromosome
    with open(out_hsq + '.nbSNPs.txt', 'w') as in_filenb:
        nsnp = preprocessing.read_grm_bin_n(in_file)
        in_filenb.write('all' + ' ' + str(nsnp) + '\n')


if __name__ == '__main__':
    main(sys.argv[1])
