#!/usr/bin/env python3

'''Compute GRMs'''

import os
# from configall import GCTA, NBPROC, USE_SBATCH  # PRU_DIR, GRM_DIR, GRM_CUTOFF,
import preprocessing
import config_dataset
import sys


def main(config_file):
    """Entry point if called as an executable"""
    config = config_dataset.config_dataset(config_file)

    in_file_pruned = os.path.join(config.pru_dir, 'all')
    out_file_grm = os.path.join(config.grm_dir, 'grm-all', 'all')
    out_file_grm_filtered = os.path.join(config.grm_dir, 'grm-all-' + str(config.grm_cutoff),
                                         'all-' + str(config.grm_cutoff))

    # create output directory if it doesn't exist
    if not os.path.exists(config.grm_dir):
        os.makedirs(config.grm_dir)

    # Using All SNPs
    preprocessing.gcta_grm(in_file=in_file_pruned,
                           par_input='--bfile',
                           out_file=out_file_grm,
                           gcta=config.gcta,
                           per_chr=True,  # If True, splits by chromosome merges them
                           ncpus=config.nbproc,
                           sbatch=config.use_sbatch)

    # filtering related individuals
    preprocessing.gcta_grm_filter(in_file=out_file_grm,
                                  out_file=out_file_grm_filtered,
                                  grm_cutoff=config.grm_cutoff,
                                  gcta=config.gcta,
                                  per_chr=True,  # If True, splits by chromosome merges them
                                  ncpus=config.nbproc,
                                  srun=config.use_sbatch)


if __name__ == '__main__':
    main(sys.argv[1])  # ,sys.argv[1:][1])
