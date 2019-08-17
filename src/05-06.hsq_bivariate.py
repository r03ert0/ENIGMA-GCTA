#!/usr/bin/env python3
"""bivariate heritability analysis"""

import os
# from configall import GCTA, MYGCTA, NBPROC, USE_SBATCH, REML_BIVAR_CALL
import preprocessing
# import ipdb
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
    # Run if .hsq file was already computed ?
    # If TRUE, will compute again the hsq,
    # otherwise won't compute the hsq if the .hsq file is present
    overwrite = False

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

    # ========= 2. All SNPS, 10 pcs, bivariate analysis =============

    out_biv = os.path.join(config.hsq_dir, 'hsq-biv')
    in_file = os.path.join(config.grm_dir, 'grm-all-0.025/all-0.025')
    treated = []

    for pheno1 in config.phe_list:
        for pheno2 in config.phe_list:

            treated.append(pheno1+'_'+pheno2)
            if (pheno1 != pheno2 and (pheno1+'.'+pheno2 not in treated) and
                    (pheno2+'.'+pheno1 not in treated)):

                pheno1path = os.path.join(config.phe_dir, pheno1 + '.txt')
                pheno2path = os.path.join(config.phe_dir, pheno2 + '.txt')
                phenopair = os.path.splitext(pheno1)[0] + '.' + os.path.splitext(pheno2)[0]

                out_biv_all = os.path.join(out_biv, phenopair)

                if config.reml_bivar_call == '--reml-bivar-no-constrain':
                    pars = var_par + ['--pheno', pheno1path,
                                      '--pheno', pheno2path,
                                      '--reml-maxit', str(200),
                                      '--reml-bivar', str(config.reml_bivar_call)]

                else:
                    pars = var_par + ['--pheno', pheno1path,
                                      '--pheno', pheno2path,
                                      '--reml-maxit', str(200),
                                      '--reml-bivar']  # , '--reml-bendV']

                # if not os.path.isfile(out_biv_all+'.hsq') or overwrite:
                #     preprocessing.gcta_hsq(in_file=in_file,
                #                            out_file=out_biv_all,
                #                            gcta=config.gcta,
                #                            mygcta=config.mygcta,
                #                            other_gcta_par=pars,
                #                            ncpus=config.nbproc,
                #                            sbatch=config.use_sbatch,
                #                            sbatch_par_j="hsq-biv",
                #                            sbatch_par_p="dedicated",  # "common",
                #                            sbatch_par_qos="fast")  # "normal")

                out_biv_rg0 = os.path.join(out_biv, 'all.rg=0.' + phenopair)
                pars_rg0 = pars + ['--reml-bivar-lrt-rg', str(0)]

                if not os.path.isfile(out_biv_rg0+'.hsq') or overwrite:
                    preprocessing.gcta_hsq(in_file=in_file,
                                           out_file=out_biv_rg0,
                                           gcta=config.gcta,
                                           mygcta=config.mygcta,
                                           ncpus=config.nbproc,
                                           other_gcta_par=pars_rg0,
                                           sbatch=config.use_sbatch,
                                           sbatch_par_j="hsq-biv",
                                           sbatch_par_p="dedicated",  # "common",
                                           sbatch_par_qos="fast")  # "normal")

                # out_biv_rg1 = os.path.join(out_biv, 'all.rg=1.' + phenopair)
                # pars_rg1 = pars + ['--reml-bivar-lrt-rg', str(1)]
                #
                # if not os.path.isfile(out_biv_rg1+'.hsq') or overwrite:
                #     preprocessing.gcta_hsq(in_file=in_file,
                #                            out_file=out_biv_rg1,
                #                            gcta=config.gcta,
                #                            mygcta=config.mygcta,
                #                            ncpus=config.nbproc,
                #                            other_gcta_par=pars_rg1,
                #                            sbatch=config.use_sbatch,
                #                            sbatch_par_j="hsq-biv",
                #                            sbatch_par_p="common",  # "common",
                #                            sbatch_par_qos="normal")  # "normal")


if __name__ == '__main__':
    main(sys.argv[1])
