#!/usr/bin/env python3

"""Prepare UKB phenotype"""

import os
import subprocess
import pandas as pd
from configtest import PHE_DIR, GEN_DIR, GCTA, PHENOTYPES, RAW_DIR

def main():
    '''Entry point if called as an executable'''
    if not os.path.exists(PHE_DIR):
        os.makedirs(PHE_DIR)
    # simulate heritable phenotypes with GCTA (all SNPs associated)
    gen_file = os.path.join(GEN_DIR, "all")
    for key, hsq in PHENOTYPES.items():
        phe_file = os.path.join(PHE_DIR, key)
        snp_file = os.path.join(GEN_DIR, "mysnps.txt")
        subprocess.run([GCTA, "--bfile", gen_file, "--simu-qt", "--simu-hsq", str(hsq),
                        "--simu-causal-loci", snp_file, "--out", phe_file], check=True)
        # add header and change extension to .txt
        phe_table = pd.read_table(phe_file + ".phen", sep=' ', names=["FID", "IID", key],
                                  index_col=False)
        phe_table.to_csv(phe_file + ".txt", sep='\t', index=False)

    # read pedigree file to obtain sex and centre phenotypes
    ped_file = os.path.join(RAW_DIR, "1000genomes",
                            "integrated_call_samples_v2.20130502.ALL.ped")
    ped_table = pd.read_table(ped_file)
    sex_table = ped_table.loc[:, ["Family ID", "Individual ID", "Gender"]]
    sex_table.columns = ["FID", "IID", "sex"]
    sex_file = os.path.join(PHE_DIR, "sex.txt")
    sex_table.to_csv(sex_file, sep='\t', index=False, header=True)
    # use the population as centre
    centre_table = ped_table.loc[:, ["Family ID", "Individual ID", "Population"]]
    centre_table.columns = ["FID", "IID", "centre"]
    centre_file = os.path.join(PHE_DIR, "centre.txt")
    centre_table.to_csv(centre_file, sep='\t', index=False, header=True)

if __name__ == '__main__':
    main()
