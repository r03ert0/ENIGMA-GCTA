#!/usr/bin/env python3
"""Compute heritability partition for genic / non-genic"""

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
    # quantitative covariables

    config = config_dataset.config_dataset(config_file)

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

    margins = [0, 10, 20, 30, 40, 50]

    # ========= 3.1 Genic/non-genic  =============

    in_genic = os.path.join(config.grm_dir, 'grm-genic')  # + -0.025 ?
    in_nongenic = os.path.join(config.grm_dir, 'grm-genic')  # + -0.025 ?
    out_genic = os.path.join(config.hsq_dir, 'hsq-genic')

    if not os.path.exists(out_genic):
        os.makedirs(out_genic)

    in_filenb = open(os.path.join(out_genic, 'genic.nbSNPs.txt'), 'w')

    for margin in margins:

        in_genic_margin = os.path.join(in_genic, 'genic-margin' + str(margin),
                                       'genic-margin' + str(margin))
        in_nongenic_margin = os.path.join(in_nongenic, 'nongenic-margin' + str(margin),
                                          'nongenic-margin' + str(margin))

        out_genic_margin = os.path.join(out_genic, 'genic-margin' + str(margin))

        # input both genic and non-genic and genic grm for variance partitioning
        print(out_genic_margin + '.test.txt')
        with open(out_genic_margin + '.test.txt', 'w') as in_file:
            in_file.write(in_genic_margin + '\n')
            in_file.write(in_nongenic_margin)

        # save the number of SNPs associated with each subgroup
        lc_genic = preprocessing.read_grm_bin_n(in_genic_margin)
        lc_nongenic = preprocessing.read_grm_bin_n(in_nongenic_margin)
        in_filenb.write('{} {}\n'.format(os.path.basename(in_genic_margin), lc_genic))
        in_filenb.write('{} {}\n'.format(os.path.basename(in_nongenic_margin), lc_nongenic))

        for pheno in config.phe_list:
            for lrt in [1, 2]:
                # --reml-lrt  1
                # Calculate the log likelihood of a reduce model with one or multiple genetic
                # variance components dropped from the full model and calculate the LRT and p-value.
                # By default, GCTA will always calculate and report the LRT for the first genetic
                # variance component, i.e. --reml-lrt 1, unless you re-specify this option,
                # e.g. --reml-lrt 2 assuming there are a least two genetic variance components
                # included in the analysis. You can also test multiple components simultaneously,
                # e.g. --reml-lrt 1 2 4. See FAQ #1 for more details.
                out_file = out_genic_margin + '.' + str(lrt) + '.' + pheno

                print(pheno)
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
                                       sbatch_par_j="hsq-genic")

    # ========= 3.2 Genic / xxk upstream and downstream / non-genic  =============

    for margin in margins:

        if margin > 0:

            in_genic_margin = os.path.join(in_genic, 'genic-margin' + str(0),
                                           'genic-margin' + str(0))
            in_updown = os.path.join(in_genic, 'updown-margin' + str(margin),
                                     'updown-margin' + str(margin))
            in_nongenic_margin = os.path.join(in_genic, 'nongenic-margin' + str(margin),
                                              'nongenic-margin' + str(margin))

            out_updown_margin = os.path.join(out_genic, 'updown-margin' + str(margin))

            # input genic, genic +/- marginkb, and non-genic grm for variance partitioning
            print(out_updown_margin + '.test.txt')
            with open(out_updown_margin + '.test.txt', 'w') as in_file:
                in_file.write(in_genic_margin + '\n')
                in_file.write(in_updown + '\n')
                in_file.write(in_nongenic_margin)

            # save the number of SNPs associated with each subgroup
            lc_genic = preprocessing.read_grm_bin_n(in_genic_margin)
            lc_nongenic = preprocessing.read_grm_bin_n(in_nongenic_margin)
            lc_updown = preprocessing.read_grm_bin_n(in_updown)
            in_filenb.write('{} {}\n'.format(os.path.basename(in_genic_margin), lc_genic))
            in_filenb.write('{} {}\n'.format(os.path.basename(in_updown), lc_updown))
            in_filenb.write('{} {}\n'.format(os.path.basename(in_nongenic_margin), lc_nongenic))

            for pheno in config.phe_list:
                for lrt in [1, 2]:
                    # --reml-lrt  1
                    # Calculate the log likelihood of a reduce model with one or multiple genetic
                    # variance components dropped from the full model and calculate the LRT and p-value.
                    # By default, GCTA will always calculate and report the LRT for the first genetic
                    # variance component, i.e. --reml-lrt 1, unless you re-specify this option,
                    # e.g. --reml-lrt 2 assuming there are a least two genetic variance components
                    # included in the analysis. You can also test multiple components simultaneously,
                    # e.g. --reml-lrt 1 2 4. See FAQ #1 for more details.
                    out_file = out_updown_margin + '.' + str(lrt) + '.' + pheno

                    print(pheno)
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
                                           sbatch_par_j="hsq-genic")





    # ========= 3.2 Genic / 0-20k and 20-50k upstream and downstream / non-genic  =============


    if (20 in margins and 50 in margins):

        in_updown1 = os.path.join(in_genic, 'updown-margin' + str(20),
                                 'updown-margin' + str(20))
        in_updown2 = os.path.join(in_genic, 'updown-margin' + "20-50",
                                 'updown-margin' + "20-50")
        in_genic_margin = os.path.join(in_genic, 'genic-margin' + str(0),
                                       'genic-margin' + str(0))
        in_nongenic_margin = os.path.join(in_genic, 'nongenic-margin' + str(50),
                                          'nongenic-margin' + str(50))

        out_updown_margin = os.path.join(out_genic, 'updown-margin' + "20-50")

        # input genic, genic +/- marginkb, and non-genic grm for variance partitioning
        print(out_updown_margin + '.test.txt')
        in_file = open(out_updown_margin + '.test.txt', 'w')
        in_file.write(in_genic_margin + '\n')
        in_file.write(in_updown1 + '\n')
        in_file.write(in_updown2 + '\n')
        in_file.write(in_nongenic_margin)
        in_file.close()

        # save the number of SNPs associated with each subgroup
        lc_updown2 = preprocessing.read_grm_bin_n(in_updown2)
        in_filenb.write('updown-margin' + "20-50" + ' ' + str(lc_updown2) + '\n')
        print(lc_updown2)

        for pheno in config.phe_list:
            for lrt in [1, 2, 3, 4]:
                # --reml-lrt  1
                # Calculate the log likelihood of a reduce model with one or multiple genetic
                # variance components dropped from the full model and calculate the LRT and p-value.
                # By default, GCTA will always calculate and report the LRT for the first genetic
                # variance component, i.e. --reml-lrt 1, unless you re-specify this option,
                # e.g. --reml-lrt 2 assuming there are a least two genetic variance components
                # included in the analysis. You can also test multiple components simultaneously,
                # e.g. --reml-lrt 1 2 4. See FAQ #1 for more details.
                out_file = out_updown_margin + '.' + str(lrt) + '.' + pheno

                print(pheno)
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
                                       sbatch_par_j="hsq-genic")

    in_filenb.close()

if __name__ == '__main__':
    main(sys.argv[1])
