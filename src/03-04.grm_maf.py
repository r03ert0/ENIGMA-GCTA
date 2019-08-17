#!/usr/bin/env python3

"""GRMs using SNPs in 3 ranges of MAF"""

import os
import preprocessing
# from configall import GCTA, NBPROC, USE_SBATCH
import config_dataset
import sys


def main(config_file):
    """Entry point if called as an executable"""

    config = config_dataset.config_dataset(config_file)

    in_file_pruned = os.path.join(config.pru_dir, 'all')
    out_dir_grm_maf = os.path.join(config.grm_dir, 'grm-maf')

    if not os.path.isdir(out_dir_grm_maf):
        os.makedirs(out_dir_grm_maf)

    maf_intervals = config.maf_intervals

    for maf_int in maf_intervals:
        print('Generating GRM for MAF interval:' + str(maf_int))
        maf_int_char = str(maf_int[0]) + '-' + str(maf_int[1])

        # Extract SNPs in the maf range
        listsnp_maf = os.path.join(out_dir_grm_maf, 'maf' + maf_int_char + '.snplist')
        preprocessing.run([config.plink,
                           '--bfile', in_file_pruned,
                           '--maf', str(maf_int[0])] +
                          (['--max-maf', str(maf_int[1]-1e-10)] if maf_int[1] != .5 else []) +
                          ['--write-snplist',
                           '--out', os.path.splitext(listsnp_maf)[0]])

        out_file_grm_mafint = os.path.join(out_dir_grm_maf, 'maf' + maf_int_char,
                                           'maf.' + maf_int_char)

        # Compute GRMs using unrelated individuals
        preprocessing.gcta_grm(in_file=in_file_pruned,
                               par_input='--bfile',
                               out_file=out_file_grm_mafint,
                               gcta=config.gcta,
                               per_chr=False,
                               ncpus=config.nbproc,
                               other_gcta_par=["--extract", listsnp_maf,
                                               "--keep", config.keep_ind],
                               sbatch=config.use_sbatch)


if __name__ == '__main__':
    main(sys.argv[1])
