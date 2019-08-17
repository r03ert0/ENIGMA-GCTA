#!/usr/bin/env python3

"""Import data from 1000 genomes"""

import os
import subprocess
import shutil
from configsimu import (GEN_DIR, FIL_DIR, PRU_DIR, GRM_DIR, GCTA, PLINK, DERIVED_DIR,
                        GRM_CUTOFF, PHE_DIR, GRM, PCS, IN_DATASET)


def main():
    """Entry point if called as an executable"""

    # genotype

    if os.path.exists(GEN_DIR):
        shutil.rmtree(GEN_DIR)

    if os.path.exists(FIL_DIR):
        shutil.rmtree(FIL_DIR)

    if os.path.exists(PRU_DIR):
        shutil.rmtree(PRU_DIR)
    os.makedirs(PRU_DIR)

    # import genotype file
    in_file = os.path.join(DERIVED_DIR, IN_DATASET, "03.pruned", "all")
    out_file = os.path.join(PRU_DIR, "all")

    subprocess.run([PLINK,
                    "--bfile", in_file,
                    "--make-bed",
                    "--out", out_file], check=True)

    # grm

    # remove output directory and create a new one
    if os.path.exists(GRM_DIR):
        shutil.rmtree(GRM_DIR)
    os.makedirs(GRM_DIR)
    os.makedirs(os.path.dirname(GRM))

    # copy grm and PCs from ukb
    in_file_grm = os.path.join(DERIVED_DIR, IN_DATASET, '04.grm', 'grm-all-' + str(GRM_CUTOFF),
                               'all-' + str(GRM_CUTOFF))
    out_file_grm = GRM

    subprocess.run([GCTA,
                    "--grm", in_file_grm,
                    "--grm-cutoff", str(GRM_CUTOFF),
                    "--make-grm",
                    "--out", out_file_grm], check=True)

    in_file_pca = os.path.join(DERIVED_DIR, IN_DATASET, '04.grm', 'grm-all-' + str(GRM_CUTOFF),
                               'all-' + str(GRM_CUTOFF) + '.pca.eigenvec')
    out_file_pca = PCS

    shutil.copyfile(in_file_pca, out_file_pca)

    # phenotype

    if not os.path.exists(PHE_DIR):
        os.makedirs(PHE_DIR)

    # copy age, sex and centre phenotypes from original dataset
    phe_dir_ukb = os.path.join(DERIVED_DIR, IN_DATASET, "00.phenotype")
    shutil.copyfile(os.path.join(phe_dir_ukb, "sex.txt"), os.path.join(PHE_DIR, "sex.txt"))
    shutil.copyfile(os.path.join(phe_dir_ukb, "centre.txt"), os.path.join(PHE_DIR, "centre.txt"))
    shutil.copyfile(os.path.join(phe_dir_ukb, "age.txt"), os.path.join(PHE_DIR, "age.txt"))


if __name__ == '__main__':
    main()
