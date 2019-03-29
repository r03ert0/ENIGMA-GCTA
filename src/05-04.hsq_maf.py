#!/usr/bin/env python3
"""Compute heritability partitioned by maf"""

import os
# from configall import GCTA, MYGCTA, NBPROC, USE_SBATCH, REML_CALL
import preprocessing
import config_dataset
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
        qcovar_par.append('--qcovar')
        qcovar_par.append(qcov)

    # qualitative covariables
    covar_par = []
    for cov in config.qual_covar:
        covar_par.append('--covar')
        covar_par.append(cov)

    var_par = qcovar_par + covar_par

    # ========= 4. MAF  =============

    in_dir_grm_maf = os.path.join(config.grm_dir, 'grm-maf')
    out_hsq_maf = os.path.join(config.hsq_dir, 'hsq-maf')
    maf_intervals = config.maf_intervals

    if not os.path.exists(out_hsq_maf):
        os.makedirs(out_hsq_maf)

    # write the number of SNPs associated with each subgroup
    in_filenb = open(os.path.join(out_hsq_maf, 'maf.nbSNPs.txt'), 'w')
    # write the GRM used to partitionate h^2
    in_file = open(os.path.join(out_hsq_maf, 'maf.test.txt'), 'w')
    for maf_int in maf_intervals:
        maf_int_char = str(maf_int[0]) + '-' + str(maf_int[1])
        in_file_grm_mafint = os.path.join(in_dir_grm_maf, 'maf' + str(maf_int_char),
                                          'maf.' + str(maf_int_char))
        in_file.write(in_file_grm_mafint + '\n')

        # extract number of SNPs
        nsnp = preprocessing.read_grm_bin_n(in_file_grm_mafint)
        in_filenb.write('maf.' + str(maf_int_char) + ' ' + str(nsnp) + '\n')

    in_file.close()
    in_filenb.close()

    for pheno in config.phe_list:

        for lrt in [1,2,3,4]:

            out_file = os.path.join(out_hsq_maf, 'maf' + '.' + str(lrt) + '.' + pheno)

            phenopath = os.path.join(config.phe_dir, pheno + '.txt')
            pars = var_par + ['--pheno', phenopath, '--reml-lrt', str(lrt), str(config.reml_call)]

            preprocessing.gcta_hsq(in_file=in_file.name,
                                   out_file=out_file,
                                   gcta=config.gcta,
                                   mygcta=config.mygcta,
                                   ncpus=config.nbproc,
                                   other_gcta_par=pars,
                                   par_input='--mgrm-bin',
                                   sbatch=config.use_sbatch,
                                   sbatch_par_j="hsq-maf")


if __name__ == '__main__':
    main(sys.argv[1])
