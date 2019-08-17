#!/usr/bin/env python3

"""Prepare UKB genotype
    The genotype files have to be decrypted and decompressed before.

    Note: we did not filter individuals base of their genotype.
    The file new_ukb_sqc_v2.txt contains columns such as Submitted.Gender, Inferred.Gender,
    sample.qc.missing.rate, heterozygosity, excluded.from.kinship.inference, excess.relatives,
    in.white.British.ancestry.subset, PCs which could be used for filtering
"""

import os
import subprocess
import pandas as pd
from configukb import PHE_DIR, GEN_DIR, RAW_DIR, PLINK


def main():
    """Entry point if called as an executable"""
    if not os.path.exists(GEN_DIR):
        os.makedirs(GEN_DIR)

    # load subjects to be analyzed
    subjects_file = os.path.join(PHE_DIR, "subjects.txt")

    # load SNP QC file to select SNPs
    snp_qc_file = os.path.join(RAW_DIR, "ukb-1226", "ukb_snp_qc.txt")
    snp_qc_table = pd.read_table(snp_qc_file, sep=' ')
    # keep only SNPs which pass QC for all batches
    keep_list = snp_qc_table.loc[snp_qc_table.filter(regex='_qc$').all(axis=1), 'rs_id']
    keep_file = os.path.join(GEN_DIR, "snp_keep.txt")
    keep_list.to_csv(keep_file, index=False)

    # extract genotypes for those individuals
    print("extracting genotypes")
    for i in range(1, 23):
        print("chr" + str(i))
        bed_file = os.path.join(RAW_DIR, "ukb-1226", "ukb_cal_chr" + str(i) + "_v2.bed")
        fam_file = os.path.join(RAW_DIR, "ukb-1226", "ukb1858_cal_v2_s488374.fam")
        bim_file = os.path.join(RAW_DIR, "ukb-1226", "ukb_snp_chr" + str(i) + "_v2.bim")
        out_file = os.path.join(GEN_DIR, "chr" + str(i))
        subprocess.run([PLINK,
                        "--bed", bed_file,
                        "--fam", fam_file,
                        "--bim", bim_file,
                        "--keep", subjects_file,
                        "--extract", keep_file,
                        "--make-bed",
                        "--out", out_file], check=True)

    # join plink files by chromosomes in one file called "all"
    print("merging chromosomes")
    out_file = os.path.join(GEN_DIR, "all")
    chr_file = os.path.join(GEN_DIR, "chr")
    with open(out_file + ".txt", "w") as thefile:
        for i in range(1, 23):
            thefile.write(chr_file + str(i) + "\n")
    subprocess.run([PLINK, "--bfile", chr_file + "1", "--merge-list", out_file + ".txt",
                    "--out", out_file], check=True)


if __name__ == '__main__':
    main()
