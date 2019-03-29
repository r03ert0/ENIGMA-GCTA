#!/usr/bin/env python3

"""GRMs using lists of genes involved in various processes"""

import os
import preprocessing
# from configall import GCTA, PLINK, NBPROC, USE_SBATCH
import config_dataset
import sys


def main(config_file):
    """Entry point if called as an executable"""
    in_margin = 50

    config = config_dataset.config_dataset(config_file)

    for namesel, filesel in [
            ('neurodev', config.selgenes_neurodev),
            ('cnsexpression', config.selgenes_cnsexpression)
    ]:
        print("Gene set name:", namesel)
        print("Gene set file:", filesel)

        in_file_pruned = os.path.join(config.pru_dir, 'all')
        in_genicregions = os.path.join(config.external_dir, 'hg19genes.txt')
        out_dir_grm_genesel = os.path.join(config.grm_dir, 'grm-' + namesel)
        out_dir_grm_genic = os.path.join(out_dir_grm_genesel, namesel + '-margin' + str(in_margin))
        out_dir_grm_nongenic = os.path.join(out_dir_grm_genesel, 'non' + namesel + '-margin' +
                                            str(in_margin))

        if not os.path.exists(out_dir_grm_genesel):
            os.makedirs(out_dir_grm_genesel)

        listsnp_genesel_margin = os.path.join(out_dir_grm_genesel, '{}-margin{}.snplist'.format(namesel, in_margin))

        # extract SNPs located in the gene set
        preprocessing.plink_make_set(bfile=in_file_pruned,
                                     border=in_margin,
                                     subset=filesel,
                                     gene_set=in_genicregions,
                                     out_file=os.path.splitext(listsnp_genesel_margin)[0],
                                     plink=config.plink
                                     )

        listsnp_nongenesel_margin = os.path.join(out_dir_grm_genesel,
                                                 'non{}-margin{}.snplist'.format(namesel, in_margin))

        # compute list of SNPs not involved in the given process
        listsnp_genic_margin = os.path.join(config.grm_dir, 'grm-genic', 'genic-margin' + str(in_margin) + '.snplist')
        snp_genic_margin = set(open(listsnp_genic_margin).readlines())
        snp_genesel_margin = set(open(listsnp_genesel_margin).readlines())
        snp_nongenesel_margin = snp_genic_margin.difference(snp_genesel_margin)
        open(listsnp_nongenesel_margin, 'w').writelines(snp_nongenesel_margin)

        # Compute GRM using SNPs in those genes using unrelated individuals
        preprocessing.gcta_grm(in_file=in_file_pruned,
                               par_input='--bfile',
                               out_file=os.path.join(out_dir_grm_genic, namesel),
                               gcta=config.gcta,
                               per_chr=False,
                               ncpus=config.nbproc,
                               other_gcta_par=["--extract", listsnp_genesel_margin,
                                               "--keep", config.keep_ind],
                               sbatch=config.use_sbatch)

        # Compute GRM using SNPs not in those genes using unrelated individuals
        preprocessing.gcta_grm(in_file=in_file_pruned,
                               par_input='--bfile',
                               out_file=os.path.join(out_dir_grm_nongenic, 'non'+namesel),
                               gcta=config.gcta,
                               per_chr=False,
                               ncpus=config.nbproc,
                               other_gcta_par=["--extract", listsnp_nongenesel_margin,
                                               "--keep", config.keep_ind],
                               sbatch=config.use_sbatch)


if __name__ == '__main__':
    main(sys.argv[1])
