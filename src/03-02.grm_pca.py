#!/usr/bin/env python3

"""PCA of GRM computed on filtered data"""

import os
# from configall import GCTA, MYGCTA, NBPROC, USE_SBATCH
import preprocessing
import config_dataset
import sys
import pandas as pd
import matplotlib as mpl
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
import matplotlib.pyplot as plt


def plot_pca(pca_file):
    png_file = os.path.splitext(pca_file)[0] + ".png"
    pca = pd.read_csv(pca_file+'.eigenvec', sep=' ', index_col=[0, 1])
    pca.plot.scatter("PC1", "PC2", alpha=0.1)
    # plt.show()
    plt.savefig(png_file)


def main(config_file):
    """Entry point if called as an executable"""
    config = config_dataset.config_dataset(config_file)

    in_grm = os.path.join(config.grm_dir, 'grm-all', 'all')
    out_file_pca = os.path.join(config.grm_dir, 'grm-all', 'all.pca')

    in_grm_filtered = os.path.join(config.grm_dir, 'grm-all-' + str(config.grm_cutoff),
                                   'all-' + str(config.grm_cutoff))
    out_file_pca_filtered = os.path.join(config.grm_dir, 'grm-all-' + str(config.grm_cutoff),
                                         'all-' + str(config.grm_cutoff) + '.pca')

    nbpcs = 10

    if config.use_sbatch:
        rmode = "srun"
    else:
        rmode = "direct"

    # compute PCA for all individuals
    preprocessing.run([config.mygcta, config.gcta,
                       "--grm-bin", in_grm,
                       "--pca", str(nbpcs),
                       "--out", out_file_pca,
                       "--thread-num", str(config.nbproc)],
                      mode=rmode,
                      slurm_par=["-J", "gcta_pca",
                                 "-p", "common,dedicated",
                                 "--qos", "fast",
                                 "-c", str(config.nbproc)])
    plot_pca(out_file_pca)

    # compute PCA for unrelated individuals
    preprocessing.run([config.mygcta, config.gcta,
                       "--grm-bin", in_grm_filtered,
                       "--pca", str(nbpcs),
                       "--out", out_file_pca_filtered,
                       "--thread-num", str(config.nbproc)],
                      mode=rmode,
                      slurm_par=["-J", "gcta_pca",
                                 "-p", "common,dedicated",
                                 "--qos", "fast",
                                 "-c", str(config.nbproc)])
    plot_pca(out_file_pca_filtered)


if __name__ == '__main__':
    main(sys.argv[1])
