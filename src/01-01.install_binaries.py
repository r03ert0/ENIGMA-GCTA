#!/usr/bin/env python3

"""Download gcta and plink executables from the web"""

import os
import sys
import stat
import re
from io import BytesIO
from zipfile import ZipFile
from urllib.request import urlopen
import config_dataset


def main(config_file):
    """Entry point if called as an executable"""
    config = config_dataset.config_dataset(config_file)

    bin_dir = os.path.join(config.annex_dir, "bin")
    bin_subdirs = next(os.walk(bin_dir))[1]

    if sys.platform == "linux":
        plink_re = re.compile(r"plink_linux_x86_64")
        plink_dirs = sorted(filter(plink_re.match, bin_subdirs))
        if not plink_dirs:
            print("Installing plink in {}".format(bin_dir))
            plink_url = "http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20190304.zip"
            resp = urlopen(plink_url)
            zipfile = ZipFile(BytesIO(resp.read()))
            zipfile.extractall(os.path.join(bin_dir, os.path.splitext(os.path.basename(plink_url))[0]))

        gcta_re = re.compile(r"gcta(?!.*_[a-z]*$)")
        gcta_dirs = sorted(filter(gcta_re.match, bin_subdirs))
        if not gcta_dirs:
            print("Installing GCTA in {}".format(bin_dir))
            gcta_url = "https://cnsgenomics.com/software/gcta/gcta_1.92.0beta3.zip"
            resp = urlopen(gcta_url)
            zipfile = ZipFile(BytesIO(resp.read()))
            zipfile.extractall(bin_dir)
        

    elif sys.platform == "darwin":
        plink_re = re.compile(r"plink_mac")
        plink_dirs = sorted(filter(plink_re.match, bin_subdirs))
        if not plink_dirs:
            print("Installing plink in {}".format(bin_dir))
            plink_url = "http://s3.amazonaws.com/plink1-assets/plink_mac_20190304.zip"
            resp = urlopen(plink_url)
            zipfile = ZipFile(BytesIO(resp.read()))
            zipfile.extractall(os.path.join(bin_dir, os.path.splitext(os.path.basename(plink_url))[0]))

        gcta_re = re.compile(r"gcta.*_mac")
        gcta_dirs = sorted(filter(gcta_re.match, bin_subdirs))
        if not gcta_dirs:
            print("Installing GCTA in {}".format(bin_dir))
            gcta_url = "https://cnsgenomics.com/software/gcta/gcta_1.92.0beta3_mac.zip"
            resp = urlopen(gcta_url)
            zipfile = ZipFile(BytesIO(resp.read()))
            zipfile.extractall(bin_dir)

    else:
        raise Exception("Unknown platform: " + sys.platform)

    print("giving executable rights to downloaded executables")
    config = config_dataset.config_dataset(config_file)
    os.chmod(config.gcta, stat.S_IEXEC)
    os.chmod(config.plink, stat.S_IEXEC)


if __name__ == '__main__':
    main(sys.argv[1])
