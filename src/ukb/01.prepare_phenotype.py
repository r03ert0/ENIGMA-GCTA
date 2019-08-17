#!/usr/bin/env python3

"""Prepare UKB phenotype"""

import os
import sys
import subprocess
from glob import glob
import numpy as np
import pandas as pd
import statsmodels.api as sm
from configukb import BIN_DIR, UKB_DATA, PHE_DIR, REGIONS, REMOVE_IND, LOG10


def ukb_run(cmd):
    """Choose the rigth way to run ukb tools depending to the system"""
    if sys.platform == 'darwin':
        # use wine
        # could use noah but didn't work last time tested
        cmd[0] = cmd[0] + ".exe"
        cmd = ["wine"] + [arg.replace('/', '\\') for arg in cmd]

    print(' '.join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def main():
    """Entry point if called as an executable"""
    # extract variables for regions
    codes = [item[2:] for sublist in REGIONS.values() for item in sublist]
    ukbconv = os.path.join(BIN_DIR, "ukb_conv")
    ukb_unpack = os.path.join(BIN_DIR, "ukbunpack")
    encoding = os.path.join(BIN_DIR, "encoding.ukb")
    codesfile = os.path.join(PHE_DIR, "varcodes.txt")
    tablefile = os.path.join(PHE_DIR, "ukbtable")
    datafile = UKB_DATA + "_ukb"

    if not os.path.exists(PHE_DIR):
        os.makedirs(PHE_DIR)

    with open(codesfile, 'w') as thefile:
        for code in codes:
            thefile.write(code + "\n")

    # use relpath to try to make the paths no longer than 64 characters
    codesfile = os.path.relpath(codesfile, os.getcwd())
    tablefile = os.path.relpath(tablefile, os.getcwd())
    encoding = os.path.relpath(encoding, os.getcwd())
    ukb_data = os.path.relpath(UKB_DATA, os.getcwd())

    if not os.path.isfile(datafile):
        # try to find the key file
        key_files = glob(os.path.join(os.path.dirname(datafile), "*.key"))
        key_file = os.path.relpath(key_files[0])
        if not key_files:
            print("unable to find the key file to decode {}.".format(UKB_DATA))
            sys.exit(1)
        ukb_run([ukb_unpack, ukb_data, key_file])

    ukb_run([ukbconv, datafile, "r", "-o"+tablefile, "-i"+codesfile, "-e"+encoding])

    table = pd.read_csv(tablefile+".tab", parse_dates=[REGIONS["date.mri"], REGIONS["date.init"]], sep='\t', dtype=str)
    table = table.set_index([table.iloc[:, 0], table.iloc[:, 0]])
    table.index.names = ['FID', 'IID']
    for key, value in REGIONS.items():
        if type(value) == list:
            daf = pd.DataFrame(table[value].sum(axis=1, skipna=False)).dropna()
            if (len(value) == 2):
                dafleft = pd.DataFrame(table.loc[:, value[0]]).dropna()
                dafright = pd.DataFrame(table.loc[:, value[1]]).dropna()
                dafleftlog10 = dafleft.apply(np.log10)
                dafrightlog10 = dafright.apply(np.log10)
                dafleftlog10.columns = [key + 'leftlog10']
                dafrightlog10.columns = [key + 'rightlog10']
                dafleftlog10.to_csv(os.path.join(PHE_DIR, key + "leftlog10.txt"), sep='\t', index=True)
                dafrightlog10.to_csv(os.path.join(PHE_DIR, key + "rightlog10.txt"), sep='\t', index=True)
        else:
            daf = pd.DataFrame(table[value]).dropna()
        daf.columns = [key]
        if key == 'ICV':
            daf = daf.astype(np.float64).rtruediv(1)
        daf.to_csv(os.path.join(PHE_DIR, key+".txt"), sep='\t', index=True)
        if key in LOG10:
            daflog10 = daf.astype(np.float64).apply(np.log10)
            daflog10.columns = [key+'log10']
            daflog10.to_csv(os.path.join(PHE_DIR, key+"log10.txt"), sep='\t', index=True)

    # Compute age from age of intial session when age of mri session is missing
    age = table[REGIONS["age.mri"]].copy()
    age.name = "age"
    diff_age = (table[REGIONS["date.mri"]] - table[REGIONS["date.init"]]
                ).dropna().map(lambda x: round(x.days/365.25, 1))
    age2 = table[REGIONS["age.init"]].astype(np.float64) + diff_age
    age = age.fillna(age2)
    age = age.dropna()
    age.to_csv(os.path.join(PHE_DIR, "age.txt"), sep='\t', header=True, index=True)

    # Compute list of subject to keep
    # Remove individual in a low density region in the ICV vs brain volume plot
    data = table[[REGIONS["ICV"], REGIONS["brain"]]].astype(np.float64).dropna()
    dens = sm.nonparametric.KDEMultivariate(data=data, var_type='cc', bw='normal_reference')
    print("computing phenotypic density to detect ouliers...")
    pdf = dens.pdf(data)
    outliers = pdf < max(pdf) / 100
    outliers_file = os.path.join(PHE_DIR, "outliers.txt")
    data.loc[outliers, :].to_csv(outliers_file, columns=[], sep='\t')
    keep_table = data.loc[~outliers, []].reset_index()
    if REMOVE_IND is not None:
        remove_list = np.loadtxt(REMOVE_IND, dtype=str)
        keep_table = keep_table.loc[~keep_table["IID"].isin(remove_list), :]

    print("saving subject file...")
    subjects_file = os.path.join(PHE_DIR, "subjects.txt")
    keep_table.to_csv(subjects_file, sep='\t', index=False)


if __name__ == '__main__':
    main()
