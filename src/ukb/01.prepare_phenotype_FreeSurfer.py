#!/usr/bin/env python3

"""Prepare UKB phenotype"""

import os
import numpy as np
import pandas as pd
from configukb import RAW_DIR, PHE_DIR, REGIONS_FREESURFER, LOG10


def main():
    """Entry point if called as an executable"""
    # extract variables for regions

    fsdata_file = os.path.join(RAW_DIR, 'ukb-FreeSurfer', 'measures-fs-ukb.tsv')
    fsdata = pd.read_table(fsdata_file, sep=r'\s+')

    fsdata.Subject.replace('_20252_2_0', value='', regex=True, inplace=True)
    fsdata = fsdata.set_index([fsdata.iloc[:, 0], fsdata.iloc[:, 0]])
    fsdata.index.names = ['FID', 'IID']

    for key, value in REGIONS_FREESURFER.items():
        if len(value) == 2:
            dafleft = pd.DataFrame(fsdata.loc[:, value[0]]).dropna()
            dafright = pd.DataFrame(fsdata.loc[:, value[1]]).dropna()
            dafleftlog10 = dafleft.apply(np.log10)
            dafrightlog10 = dafright.apply(np.log10)
            dafleftlog10.columns = [key + 'leftlog10']
            dafrightlog10.columns = [key + 'rightlog10']
            dafleftlog10.to_csv(os.path.join(PHE_DIR, key + "leftlog10_freesurfer.txt"), sep='\t', index=True)
            dafrightlog10.to_csv(os.path.join(PHE_DIR, key + "rightlog10_freesurfer.txt"), sep='\t', index=True)

        daf = pd.DataFrame(fsdata.loc[:, value].sum(axis=1, skipna=False), columns=[key]).dropna()
        if key in LOG10:
                daflog10 = daf.apply(np.log10)
                daflog10.columns = [key+'log10']
                daflog10.to_csv(os.path.join(PHE_DIR, key+"log10_freesurfer.txt"), sep='\t', index=True)


if __name__ == '__main__':
    main()
