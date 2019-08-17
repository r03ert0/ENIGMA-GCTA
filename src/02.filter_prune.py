#!/usr/bin/env python3

"""Filtering/pruning of the SNP data"""

import os
import sys
# from configall import PLINK  # , DERIVED_DIR, DATASET, EXCLUDE_IND, EXCLUDE_RANGE
import config_dataset
import preprocessing


def main(config_file):
    """Entry point if called as an executable"""
    config = config_dataset.config_dataset(config_file)

    print(config.plink)

    in_dir = os.path.join(config.derived_dir, config.dataset)
    in_prefix = 'all'

    if config.exclude_ind is not None:
        other_plink_par = ["--remove", config.exclude_ind]
    else:
        other_plink_par = []

    if config.exclude_range is not None:
        other_plink_par = ["--exclude", "range", config.exclude_range] + other_plink_par

    preprocessing.filter_prune(
        in_dir=in_dir,
        in_prefix=in_prefix,
        plink=config.plink,
        other_plink_par=other_plink_par,
        maf=config.maf)


if __name__ == '__main__':
    main(sys.argv[1])
