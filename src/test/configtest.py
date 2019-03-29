"""Make config"""

import subprocess
import os

def getannexdir():
    """call git executable to find annex dir"""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    return subprocess.run(["git", "rev-parse", "--show-toplevel"], stdout=subprocess.PIPE,
                          cwd=script_dir).stdout.decode().strip()


ANNEX_DIR = globals().get('ANNEX_DIR', getannexdir())
PLINK = 'plink'
GCTA = 'gcta64'

# indicate which dataset to use
DATASET = 'test'

# data subdirectories
RAW_DIR = os.path.join(ANNEX_DIR, 'data/raw')
DERIVED_DIR = os.path.join(ANNEX_DIR, 'data/derived')
EXTERNAL_DIR = os.path.join(ANNEX_DIR, 'data/external')

# directory containing all phenotypes
PHE_DIR = os.path.join(DERIVED_DIR, DATASET, '00.phenotype')

# derived data subdirectories
GEN_DIR = os.path.join(DERIVED_DIR, DATASET, '01.genotype')

# number of processors
NBPROC = 4

# plink input file
IN_FILE = os.path.join(RAW_DIR, "1000genomes",
                       "ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes")

# number of extracted SNPs
N_SNP = 20000

# simulated phenotypes and their heritability
PHENOTYPES = {
    "height": .6,
    "intelligence": .3,
    "age": 0,
    "ICV": .4,
    "brain": .5,
    "thalamus": .4,
    "caudate": .4,
    "putamen": .4,
    "pallidum": .4,
    "hippocampus": .4,
    "amygdala": .4,
    "accumbens": .4
}
