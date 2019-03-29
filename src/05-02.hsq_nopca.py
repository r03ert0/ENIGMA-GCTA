#!/usr/bin/env python3
"""compute heritability without including PCs as covariates"""

# from future import print_function
import os
# from configall import GCTA, MYGCTA, NBPROC, USE_SBATCH, REML_CALL
import preprocessing
import config_dataset
# import ipdb
import sys


# ipdb.set_trace()
# C-c C-z : open a python shell
# C-c C-c : run the content of the buffer in the opened python shell
# C-c C-r : run the selected region in the python shell
# Simple command
# subprocess.call(['ls', '-1'], shell=True)
# http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html

def main(config_file):
    """Entry point if called as an executable"""

    config = config_dataset.config_dataset(config_file)
    # quantitative covariables
    qcovar_par = []
    for qcov in config.quant_covar:
        if qcov == config.pcs:
            continue
        qcovar_par.append('--qcovar')
        qcovar_par.append(qcov)

    # qualitative covariables
    covar_par = []
    for cov in config.qual_covar:
        covar_par.append('--covar')
        covar_par.append(cov)

    var_par = qcovar_par + covar_par

    # ========= 1bis. All SNPS, no PC =============

    in_file = os.path.join(config.grm_dir, 'grm-all/all')
    out_dir = os.path.join(config.hsq_dir, 'hsq-nopca', 'nopca')

    for pheno in config.phe_list:
        out_file = out_dir+'.'+pheno
        print('running h^2 gcta estimation (no PCs) for phenotype: ' + pheno)
        pheno = os.path.join(config.phe_dir, pheno + '.txt')
        pars = var_par + ['--pheno', pheno, str(config.reml_call)]

        preprocessing.gcta_hsq(in_file=in_file,
                               out_file=out_file,
                               gcta=config.gcta,
                               mygcta=config.mygcta,
                               ncpus=config.nbproc,
                               other_gcta_par=pars,
                               sbatch=config.use_sbatch,
                               sbatch_par_j="hsq-nopca")


if __name__ == '__main__':
    main(sys.argv[1])
