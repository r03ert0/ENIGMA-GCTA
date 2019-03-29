#!/usr/bin/env python3

"""
    Master script that runs in the appropriate order all the scripts necessary to
    convert raw data to generic data for test dataset.
    @Author Nicolas Traut, 1 November 2017
"""

import subprocess

subprocess.run("./01.prepare_genotype.py", check=True)
subprocess.run("./02.prepare_phenotype.py", check=True)
