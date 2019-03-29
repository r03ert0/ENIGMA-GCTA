#!/usr/bin/env python3
"""Run permutations to test for enrichment of partitions """

import preprocessing
import os
import config_dataset
import sys
import math


def main(config_file):
    config = config_dataset.config_dataset(config_file)

    nbit = config.nbit
    nbproc = config.nbproc

    out_dir = config.permu_dir
    out_dir_log = os.path.join(out_dir, 'log')
    if not os.path.exists(out_dir_log):
        os.makedirs(out_dir_log)

    if config.use_sbatch:
        mode = "sbatch"
    else:
        mode = "direct"

    # Run permutations
    def run_permutations():

        in_file_pruned = os.path.join(config.pru_dir, 'all')
        
        # Estimate required memory in kilobytes
        nind = preprocessing.linecount(config.keep_ind)
        nsnp_tot = preprocessing.linecount(in_file_pruned + '.bim')
        nsnps = [preprocessing.linecount(f) for f in snplist_files]
        if nsnp_tot == sum(nsnps):
            ncomp = 1 + len(snplist_files)
        else:
            ncomp = 2 + len(snplist_files)
        mem = math.ceil(500000 + (0.025 + 0.009 * ncomp) * nind**2)
        sbatch_par_mem = str(mem) + 'K'
        
        array_lim = config.array_lim

        job_id = preprocessing.run_array(
                ['python3',
                 os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              'permutations.py'),
                 '\\j' + '_' + '\\i',
                 in_file_pruned,
                 ",".join(snplist_files),
                 permu_path,
                 os.path.abspath(config_file)],
                mode=mode,
                slurm_par=["-J", "permu_pheno",
                           # "--qos", "fast",
                           "-p", "common",
                           "--mem", sbatch_par_mem,
                           "--cpus-per-task", str(nbproc),
                           "-D", out_dir_log],
                array=range(1, nbit+1),
                array_limit=array_lim)
        return job_id

    part2jid = {}

    # margins = [0, 10, 20, 30, 40, 50]
    margins = [0, 20, 50]

    # genic non-genic
    for margin in margins:
        snplist_files = [os.path.join(config.grm_dir, 'grm-genic', 'genic-margin' + str(margin) + '.snplist')]
        permu_path = os.path.join(out_dir, 'genic-margin' + str(margin))
        hsq_prefix = os.path.join("hsq-genic", "genic-margin" + str(margin))

        part2jid[hsq_prefix] = run_permutations()

        if margin > 0:
            snplist_files = [os.path.join(config.grm_dir, 'grm-genic', 'genic-margin0.snplist'),
                             os.path.join(config.grm_dir, 'grm-genic', 'updown-margin' + str(margin) + '.snplist')]
            permu_path = os.path.join(out_dir, 'updown-margin' + str(margin))
            hsq_prefix = os.path.join("hsq-genic", "updown-margin" + str(margin))

            part2jid[hsq_prefix] = run_permutations()

    # cnsexpression non-cnsexpression non-genic
    snplist_files = [os.path.join(config.grm_dir, 'grm-cnsexpression', 'cnsexpression-margin50.snplist'),
                     os.path.join(config.grm_dir, 'grm-cnsexpression', 'noncnsexpression-margin50.snplist')]
    permu_path = os.path.join(out_dir, "cnsexpression-margin50")
    hsq_prefix = os.path.join("hsq-cnsexpression", "cnsexpression-margin50")

    part2jid[hsq_prefix] = run_permutations()

    # neurodev non-neurodev non-genic
    snplist_files = [os.path.join(config.grm_dir, 'grm-neurodev', 'neurodev-margin50.snplist'),
                     os.path.join(config.grm_dir, 'grm-neurodev', 'nonneurodev-margin50.snplist')]
    permu_path = os.path.join(out_dir, "neurodev-margin50")
    hsq_prefix = os.path.join("hsq-neurodev", "neurodev-margin50")

    part2jid[hsq_prefix] = run_permutations()

    # maf
    maf_intervals = config.maf_intervals
    snplist_files = [os.path.join(config.grm_dir, 'grm-maf', 'maf{}-{}.snplist'.format(*maf_int)) for maf_int
                     in maf_intervals]
    permu_path = os.path.join(out_dir, "maf")
    hsq_prefix = os.path.join("hsq-maf", "maf")

    part2jid[hsq_prefix] = run_permutations()

    # == Once permutations are done, compute z-scores and p-value for the different partitions == #

    # to compute the p-values with sbatch
    for partition, jid in part2jid.items():
        cmd = ["python3",
               os.path.join(os.path.dirname(os.path.abspath(__file__)), 'permutations_zscores.py'),
               os.path.join(os.path.dirname(os.path.abspath(__file__)), config_file),
               partition,
               out_dir]
        slurm_par = ["-J", "zscores_" + partition,
                     "--qos", "fast",
                     "-p", "dedicated",
                     # "-p", "common",
                     "-D", out_dir_log,
                     "--mem", "2G",
                     "--dependency", "afterany:" + jid]
        preprocessing.run(cmd, mode=mode, slurm_par=slurm_par)

    # # to compute the p-values without sbatch
    # partitions = ["hsq-genic/genic-margin0",
    #               "hsq-genic/genic-margin20",
    #               "hsq-genic/genic-margin50",
    #               "hsq-cnsexpression/cnsexpression-margin50",
    #               "hsq-neurodev/neurodev-margin50"]
    #
    # for partition in partitions:
    #     permutations_zscores.process_permu_r(configs, partition)


if __name__ == '__main__':
    main(sys.argv[1])
