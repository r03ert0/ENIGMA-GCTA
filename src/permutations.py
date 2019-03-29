import subprocess
import os
import sys
import struct
import preprocessing
import pandas as pd
import numpy as np
import glob
import config_dataset


def random_grm(in_file,
               out_dir,
               snplist_files,
               permu_type='simple',
               out_file_prefix='',
               # level='variant',
               # plink='plink',
               gcta='gcta',
               # mygcta='mygcta',
               other_gcta_par_grm=None,
               # other_gcta_par_hsq=None,
               ncpus=1,
               sbatch=False):
    """
    randomly sample sets of SNPs from a GCTA-format binary file and compute the associated GRMs
    Args:
        in_file (str): bfile fileset referencing plink files bfile.bed + bfile.bim + bfile.fam
        out_dir (str): output directory for grm and hsq files
        snplist_files (list): list of the SNPs of each partition to be randomly resampled without overlap
        permu_type (str): permutation method. Can be "simple" or "bins".
        out_file_prefix (str): prefix for out file
        # level (str): either 'variant' to sample variants, or 'sample' to sample samples
        # plink (str): path of plink executable
        gcta (str): path of gcta executable
        # mygcta (str): path of mygcta executable
        other_gcta_par_grm (list): list of additional gcta parameters to be used for the GRM analysis
        ncpus (int): number of CPUs to be used
        sbatch (bool): whether to run into sbatch mode

    Returns:
        The path of a file containing the paths of the GRMs computed by GCTA on the random sets of SNPs

    """

    out_dir_grm = os.path.join(out_dir, 'grm')

    if not os.path.exists(out_dir_grm):
        os.makedirs(out_dir_grm, exist_ok=True)

    # read SNP list
    all_snps = pd.read_csv(in_file + '.bim', sep='\t', header=None).iloc[:, 1].values
    partition = np.zeros(all_snps.size, dtype=int)
    for i, snplist_file in enumerate(snplist_files):
        part_snps = pd.read_csv(snplist_file, sep='\t', header=None).iloc[:, 0].values
        isin_part = np.in1d(all_snps, part_snps) * (i + 1)
        if (partition & isin_part).any():
            raise Exception('Error: overlapping SNP lists')
        partition = partition + isin_part
    partition[np.argwhere(partition == 0)] = len(snplist_files) + 1
    nb_part = partition.max()

    if permu_type == 'simple':
        partition2 = np.random.permutation(partition)

    elif permu_type == 'bins':
        bloc_ends = np.where(np.concatenate((partition[:-1] != partition[1:], [True])))[0]
        bloc_lengths = np.diff(np.concatenate(([0], bloc_ends + 1)))
        bloc_part = partition[bloc_ends]
        bloc_part_lengths = [bloc_lengths[bloc_part == i] for i in range(1, nb_part + 1)]

        nb_snps = all_snps.size
        nb_blocs = bloc_part.size
        print("{} blocks for {} SNPs.".format(nb_blocs, nb_snps))

        # shuffle blocs
        bloc_part2 = np.roll(bloc_part, np.random.randint(nb_blocs))
        bloc_part_lengths2 = [np.random.permutation(b) for b in bloc_part_lengths]
        bloc_lengths2 = np.zeros(nb_blocs, dtype=int)
        for i in range(1, nb_part + 1):
            bloc_lengths2[bloc_part2 == i] = bloc_part_lengths2[i-1]
        partition2 = np.repeat(bloc_part2, bloc_lengths2)

    else:
        raise ValueError("invalid value for permu_type parameter")

    out_grm_files = []

    # compute grm for each random partition
    for i in range(1, nb_part + 1):

        out_file_grm = os.path.join(out_dir_grm, out_file_prefix + '_random' + str(i))
        snplist_file = out_file_grm + '.snplist'

        # print(nb)
        snps = all_snps[partition2 == i]

        np.savetxt(X=snps,
                   fname=snplist_file,
                   fmt='%s')

        # compute grm using this set of sampled SNPs only
        preprocessing.gcta_grm(in_file=in_file,
                               out_file=out_file_grm,
                               gcta=gcta,
                               per_chr=False,
                               ncpus=ncpus,
                               other_gcta_par=["--extract",
                                               snplist_file] + other_gcta_par_grm,
                               sbatch=sbatch)

        out_grm_files.append(out_file_grm)

    return out_grm_files


def main(it, plink_file, snplist_files, permu_path, config_file):
    """ Called  as an executable when file is called"""

    # === configurations of paths === #
    config = config_dataset.config_dataset(config_file)

    snplist_files = snplist_files.split(',')

    if not os.path.exists(permu_path):
        os.makedirs(permu_path, exist_ok=True)

    print('snplist_files')
    print(snplist_files)

    other_gcta_par_grm = ["--keep", config.keep_ind]

    qcovar_par = []
    for qc in config.quant_covar:
        qcovar_par.append('--qcovar')
        qcovar_par.append(qc)

    covar_par = []
    for cv in config.qual_covar:
        covar_par.append('--covar')
        covar_par.append(cv)

    var_par = qcovar_par + covar_par

    save_permut_file = os.path.join(permu_path, 'hsq', 'it' + str(it) + '_res_hsq' + '.csv')
    out_file_prefix = 'it' + str(it)
    run_one_permutation(in_file=plink_file,
                        out_dir=permu_path,
                        out_file_prefix=out_file_prefix,
                        snplist_files=snplist_files,
                        # plink=config.plink,
                        gcta=config.gcta,
                        mygcta=config.mygcta,
                        other_gcta_par_grm=other_gcta_par_grm,
                        other_gcta_par_hsq=var_par,
                        ncpus=config.nbproc,
                        # sbatch=False,
                        phe_list=config.phe_list,
                        phe_dir=config.phe_dir,
                        save_permut_file=save_permut_file,
                        permu_type=config.permu_type,
                        cleanup=config.clean_permu)


def run_one_permutation(in_file,
                        out_dir,
                        snplist_files,
                        phe_list,
                        out_file_prefix='',
                        # plink='plink',
                        gcta='gcta',
                        mygcta='mygcta',
                        other_gcta_par_grm=None,
                        other_gcta_par_hsq=None,
                        ncpus=1,
                        # sbatch=False,
                        phe_dir='',
                        save_permut_file=None,
                        permu_type="simple",
                        cleanup=True):
    """Computes heritability on SNPs sampled randomly.

    Randomly sample n SNPs from a GCTA-format binary file,
    compute the GRM based on those SNPs and the GRM based on all other SNPs,
    compute the heritability estimates with GCTA.

    Args:
        in_file (str): bfile fileset referencing plink files bfile.bed + bfile.bim + bfile.fam
        out_dir (str): output directory for grm and hsq files
        snplist_files (list): list of the SNPs of each partition to be randomly resampled without overlap
        out_file_prefix (str): prefix for out file
        # plink (str): path of plink executable
        gcta (str): path of gcta executable
        mygcta (str): path of mygcta script
        other_gcta_par_grm (list): list of additional gcta parameters to be used for the GRM analysis
        other_gcta_par_hsq (list): list of additional gcta parameters to be used for the heritability analysis
        ncpus (int): number of CPUs to be used
        # sbatch (bool): TRUE to run the GRM in sbatch mode (FALSE is default)
        phe_dir (str): path of the phenotype data
        phe_list (list): list of phenotypes available in phe_dir that should be considered for heritability analyses
        save_permut_file (str): file in which all permutation result will be saved,
            if None out_dir/hsq/res_hsq_permutations.csv will be chosen.
        permu_type (str): permutation method. Can be "simple" or "bins".
        cleanup (bool): erase intermediate files at the end.
    Returns:
        A grm file name computed by GCTA on the randomly sampled SNPs or samples

    """

    # # total number of available variants
    # nb_tot = preprocessing.linecount(in_file + '.bim')
    # # number of SNP partitions
    # nb_part = len(snplist_files) + 1

    h2_res = pd.DataFrame()
    # print(snplist_files)

    # === Compute GRM for a random set of SNPs ===
    grm_files = random_grm(in_file=in_file,
                           out_dir=out_dir,
                           out_file_prefix=out_file_prefix,
                           snplist_files=snplist_files,
                           permu_type=permu_type,
                           # plink=plink,
                           gcta=gcta,
                           ncpus=ncpus,
                           other_gcta_par_grm=other_gcta_par_grm)

    # === Compute heritability for this random GRM for each phenotype of interest  ===

    # file where permutation results will be saved
    dir_hsq = os.path.join(out_dir, 'hsq')

    if not os.path.exists(dir_hsq):
        os.makedirs(dir_hsq, exist_ok=True)

    if save_permut_file is None:
        save_permut_file = dir_hsq + 'res_hsq_permutations' + '.csv'

    if os.path.exists(save_permut_file):
        header = False
    else:
        header = True

    out_dir_hsq_filename = os.path.join(dir_hsq, out_file_prefix)

    # input random GRMs for variance partitioning
    mgrm_file = out_dir_hsq_filename + '.mgrm.txt'
    np.savetxt(X=grm_files,
               fname=mgrm_file,
               fmt='%s')

    # extract numbers of SNPs for each GRM:
    n = [int(struct.unpack('f', open(g + '.grm.N.bin', 'rb').read(4))[0]) for g in grm_files]

    for pheno in phe_list:

        out_file = out_dir_hsq_filename + '.' + pheno
        print('running h^2 gcta estimation for phenotype: ' + pheno)
        phenofile = os.path.join(phe_dir, pheno + '.txt')
        pars = other_gcta_par_hsq + ['--pheno', phenofile, '--reml-no-constrain']  # '--reml-lrt', '1',

        # compute h^2 for the snp partitions
        try:
            subprocess.check_call([mygcta, gcta,
                                   '--mgrm-bin', mgrm_file,
                                   "--out", out_file,
                                   "--thread-num", str(ncpus)] + pars)
        # if not converging, try next phenotype
        except subprocess.CalledProcessError:
            # remove hsq log file as log file of slurm already contains error info
            os.remove(out_file + '.log')
            continue

        # extract h^2 and standard error estimates for each partition
        res = pd.read_csv(out_file + '.hsq', sep='\t')
        res_var = res.iloc[:, 0:2]
        res_se = res.iloc[:, [2]]
        res_se = res_se.assign(Source=[i + '_se' for i in res_var.Source.values])
        res_se = res_se.iloc[:, [1, 0]]
        res_var.columns = ['variable', 'value']
        res_se.columns = ['variable', 'value']
        res = res_var.append(res_se)
        res = res.transpose()
        res['pheno'] = pheno
        res.loc['variable', 'pheno'] = 'pheno'
        for i, nb in enumerate(n):
            res['n' + str(i + 1)] = ['n' + str(i + 1), str(nb)]
        res.columns = res.loc['variable', :].values
        res = res.iloc[[1], :]
        h2_res = h2_res.append(res)

    if cleanup:
        # remove hsq  and log files
        for f in glob.glob('{}.*.hsq'.format(out_dir_hsq_filename)) \
                 + glob.glob('{}.*.log'.format(out_dir_hsq_filename)):
            os.remove(f)

        # remove grm files
        os.remove(mgrm_file)
        for f in grm_files:
            for f2 in glob.glob('{0}*'.format(f)):
                os.remove(f2)

    with open(save_permut_file, 'a') as f:
        h2_res.to_csv(path_or_buf=f, sep=';',
                      header=header,
                      index=False,
                      na_rep='NA')

    return h2_res


if __name__ == '__main__':
    main(*sys.argv[1:])
