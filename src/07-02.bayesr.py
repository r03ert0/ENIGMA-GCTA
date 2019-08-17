#!/usr/bin/env python3
"""Compute heritability for all SNPs"""

import os
import sys
import math
import pandas as pd
import statsmodels.formula.api as smf
import preprocessing
import config_dataset


def main(config_file):

    """Entry point if called as an executable"""

    config = config_dataset.config_dataset(config_file)

    # numit = 10000
    # burnin = 5000
    seed = 333
    ndist = 4
    gpin = [0., 0.00001, 0.0001, 0.001]

    # slurm configuration
    if config.use_sbatch:
        mode = "sbatch"
    else:
        mode = "direct"

    ncpus = config.nbproc

    out_dir = os.path.join(config.hsq_dir, 'bayesR')
    log_dir = os.path.join(out_dir, 'log')
    tmp_dir = os.path.join(out_dir, 'plink')

    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(tmp_dir, exist_ok=True)

    in_files = {"all": os.path.join(config.fil_dir, 'all'), "pruned": os.path.join(config.pru_dir, 'all')}

    for key in in_files:

        in_file = in_files[key]

        for pheno in config.phe_list:
            pheno_file = os.path.join(config.phe_dir, pheno+'.txt')
            tmp_in_file = os.path.join(tmp_dir, key + "." + pheno)
            out_file = os.path.join(out_dir, key + "." + pheno)

            data = pd.read_table(pheno_file)
            data.rename(columns={data.columns[2]: "pheno"}, inplace=True)

            for qcov in config.quant_covar:
                data = data.merge(pd.read_table(qcov, sep='\s+'))
            for cov in config.qual_covar:
                data = data.merge(pd.read_table(cov, sep='\s+'))

            data.set_index(["FID", "IID"], inplace=True)

            model = smf.ols(formula='pheno~' + "+".join(data.columns.difference(["pheno"])), data=data).fit()
            resid_file = tmp_in_file + ".resid.txt"
            model.resid.to_csv(resid_file, sep=" ")

            preprocessing.run([config.plink,
                               "--bfile", in_file,
                               "--keep", config.keep_ind,
                               "--pheno", resid_file,
                               "--make-bed",
                               "--out", tmp_in_file
                               ])

            # Estimate required memory in megabytes
            nind = preprocessing.linecount(tmp_in_file + '.fam')
            nsnp = preprocessing.linecount(tmp_in_file + '.bim')
            mem = math.ceil((1 + 2e-6 * nind * nsnp) * 1.1)
            sbatch_par_mem = str(mem) + 'M'

            cmd = [config.bayesrv2,
                   "-bfile", tmp_in_file,
                   "-nthreads", str(ncpus),
                   "-ndist", str(ndist),
                   "-gpin", ",".join(map(str, gpin)),
                   "-out", out_file,
                   # "-numit", str(numit),
                   # "-burnin", str(burnin),
                   "-seed", str(seed)]
            slurm_par = ["-J", "bayesR_" + key,
                         "--qos", "ghfc",
                         "-p", "ghfc",
                         "-c", str(ncpus),
                         "-D", log_dir,
                         "-o", key + "." + pheno + "-%j.out",
                         "-e", key + "." + pheno + "-%j.out",
                         "--mem", sbatch_par_mem]
            jid = preprocessing.run(cmd, mode=mode, slurm_par=slurm_par)

            if config.clean_permu:
                cmd = ["rm", tmp_in_file + ".*"]
                slurm_par = ["-J", "clean_bayesR",
                             "-p", "common,dedicated",
                             "--qos", "fast",
                             "-D", log_dir,
                             "-o", "clean." + key + "." + pheno + "-%j.out",
                             "-e", "clean." + key + "." + pheno + "-%j.out",
                             "--mem", "500M",
                             "--dependency", "afterany:" + jid]
                preprocessing.run(cmd, mode=mode, slurm_par=slurm_par)


if __name__ == '__main__':
    main(sys.argv[1])
