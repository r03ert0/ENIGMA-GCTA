#!/usr/bin/env python3
"""Compute permutation z-scores"""


import preprocessing
import os
import config_dataset
import sys
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages


def process_permu_r(config, partition, out_dir):

    r_source = robjects.r['source']

    if os.path.exists("process_permu.R"):
        r_source("process_permu.R")
    else:
        r_source(os.path.join(os.path.dirname(os.path.abspath(__file__)), "process_permu.R"))

    process_permu = robjects.globalenv['process_permu']

    path_permu = os.path.join(out_dir, os.path.basename(partition), 'hsq')

    prefix_permu = '_res_hsq.csv'
    path_res = os.path.join(config.hsq_dir, os.path.dirname(partition))
    prefix_res = os.path.basename(partition) + '.1.*hsq'

    hsq_sum = os.path.join(config.hsq_dir, 'hsq-summary')

    if not os.path.exists(hsq_sum):
        os.makedirs(hsq_sum)

    output_file = os.path.join(hsq_sum, 'permu_pval_' + os.path.basename(partition) + '.txt')

    process_permu(path_permu, prefix_permu, path_res, prefix_res, output_file, 'FALSE')


def main(config_file, partition, out_dir=None):
    config = config_dataset.config_dataset(config_file)
    if out_dir is None:
        out_dir = os.path.join(config.hsq_dir, 'permutations')
    process_permu_r(config, partition, out_dir)


if __name__ == '__main__':

    main(*sys.argv[1:])
