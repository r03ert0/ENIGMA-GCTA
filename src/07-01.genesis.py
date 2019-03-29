#!/usr/bin/env python3
"""Compute heritability for all SNPs"""

import os
import sys
import preprocessing
import config_dataset


def main(config_file):

    """Entry point if called as an executable"""

    config = config_dataset.config_dataset(config_file)

    # ========= 1. All SNPS =============

    in_dir_gwas_allsnps = os.path.join(config.gwa_dir, 'gwas-all')
    out_dir = os.path.join(config.hsq_dir, 'genesis')
    out_dir_log = os.path.join(out_dir, 'log')

    os.makedirs(out_dir_log, exist_ok=True)

    # slurm configuration
    if config.use_sbatch:
        mode = "sbatch"
    else:
        mode = "direct"

    for pheno in config.phe_list:
        assoc_file = os.path.join(in_dir_gwas_allsnps, "all." + pheno + ".assoc.linear")
        res_file = os.path.join(out_dir, "all." + pheno + ".RData")

        if not os.path.exists(assoc_file):
            print("Warning: {} not found.".format(assoc_file))
            continue

        cmd = ["Rscript",
               os.path.join(os.path.dirname(os.path.abspath(__file__)), 'genesis.R'),
               assoc_file,
               res_file,
               str(config.nbproc)]
        slurm_par = ["-J", "genesis",
                     "--qos", "ghfc",
                     "-p", "ghfc",
                     # "-p", "common",
                     "-D", out_dir_log,
                     "-o", "all." + pheno + "-%j.out",
                     "-e", "all." + pheno + "-%j.out",
                     "-c", str(config.nbproc),
                     "--mem", "4G"]
        preprocessing.run(cmd, mode=mode, slurm_par=slurm_par)


if __name__ == '__main__':
    main(sys.argv[1])
