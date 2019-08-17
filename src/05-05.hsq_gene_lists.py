#!/usr/bin/env python3
"""Compute heritability partition by gene list"""

# from future import print_function
import os
# from configall import GCTA, MYGCTA, NBPROC, USE_SBATCH, REML_CALL
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


def main(dataset):
    """Entry point if called as an executable"""
    # quantitative covariables

    config = config_dataset.config_dataset(dataset)
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

    # ========= 4. gene lists =============

    margin = 50

    for namesel, grmsel in [('neurodev', os.path.join(config.grm_dir, 'grm-neurodev')),
                            ('cnsexpression', os.path.join(config.grm_dir, 'grm-cnsexpression'))
                            ]:

        print("Gene set name:", namesel)
        print("GRMs:", grmsel)

        in_genesel = os.path.join(config.grm_dir, 'grm-' + namesel)  # + -0.025 ?
        out_genesel = os.path.join(config.hsq_dir, 'hsq-' + namesel)  # + '-margin' + str(margin))

        in_genesel_margin = os.path.join(in_genesel,
                                         namesel + '-margin' + str(margin),
                                         namesel)
        in_nongenesel_margin = os.path.join(in_genesel,
                                            'non' + namesel + '-margin' + str(margin),
                                            'non' + namesel)
        in_nongenic_margin = os.path.join(config.grm_dir, 'grm-genic', 'nongenic-margin' + str(margin),
                                          'nongenic-margin' + str(margin))

        out_genesel_margin = os.path.join(out_genesel, namesel + '-margin' + str(margin))

        print("hsq outputs:", out_genesel)

        if not os.path.exists(out_genesel):
            os.makedirs(out_genesel)

        # input both snp in gene lists and snp not in gene lists grm for variance partitioning
        with open(out_genesel_margin + '.test.txt', 'w') as in_file:
            in_file.write(in_genesel_margin + '\n')
            in_file.write(in_nongenesel_margin + '\n')
            in_file.write(in_nongenic_margin + '\n')

        # save the number of SNPs associated with each subgroup
        lc_genesel = preprocessing.read_grm_bin_n(in_genesel_margin)
        lc_nongenesel = preprocessing.read_grm_bin_n(in_nongenesel_margin)
        lc_nongenic = preprocessing.read_grm_bin_n(in_nongenic_margin)
        in_filenb = open(os.path.join(out_genesel, namesel + '.nbSNPs.txt'), 'w')
        in_filenb.write('{} {}\n'.format(os.path.basename(in_genesel_margin), lc_genesel))
        in_filenb.write('{} {}\n'.format(os.path.basename(in_nongenesel_margin), lc_nongenesel))
        in_filenb.write('{} {}\n'.format(os.path.basename(in_nongenic_margin), lc_nongenic))
        in_filenb.close()

        for pheno in config.phe_list:

            print('running h^2 gcta estimation for phenotype: ' + pheno)

            for lrt in [1, 2, 3]:
                out_file = out_genesel_margin + '.' + str(lrt) + '.' + pheno

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
                                       sbatch_par_j="hsq-genesel")

            # Test partitioning versus no partitioning ?
            # out_file = os.path.join(out_genesel, namesel + '.' + '123' + '.' + pheno)

            # phenopath = os.path.join(PHE_DIR, pheno + '.txt')
            # pars = var_par + ['--pheno', phenopath, '--reml-lrt 1 2 3', '--reml']

            # preprocessing.gcta_hsq(in_file=in_file.name,
            #                        out_file=out_file,
            #                        gcta=GCTA,
            #                        mygcta=MYGCTA,
            #                        ncpus=NBPROC,
            #                        other_gcta_par=pars,
            #                        par_input='--mgrm-bin',
            #                        sbatch=USE_SBATCH)


if __name__ == '__main__':
    main(sys.argv[1])
