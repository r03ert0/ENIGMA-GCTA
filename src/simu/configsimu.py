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
GCTB = 'gctb'

# indicate which dataset to use
IN_DATASET = 'ukb-22110_maf001'
DATASET = 'simulations.ukb-22110_maf001'

# data subdirectories
RAW_DIR = os.path.join(ANNEX_DIR, 'data/raw')
DERIVED_DIR = os.path.join(ANNEX_DIR, 'data/derived')
EXTERNAL_DIR = os.path.join(ANNEX_DIR, 'data/external')

# derived data subdirectories
PHE_DIR = os.path.join(DERIVED_DIR, DATASET, '05.phenotype')
GEN_DIR = os.path.join(DERIVED_DIR, DATASET, '01.genotype')
FIL_DIR = os.path.join(DERIVED_DIR, DATASET, '02.filtered')
PRU_DIR = os.path.join(DERIVED_DIR, DATASET, '03.pruned')
GRM_DIR = os.path.join(DERIVED_DIR, DATASET, '04.grm')
GWA_DIR = os.path.join(DERIVED_DIR, DATASET, '05.gwas')
HSQ_DIR = os.path.join(DERIVED_DIR, DATASET, '06.hsq')

# GRM cutoff
GRM_CUTOFF = 0.025

# wrapper for gcta
MYGCTA = os.path.join(ANNEX_DIR, "bin", "mygcta", "mygcta.sh")

# wrapper for plink
MYPLINK = os.path.join(ANNEX_DIR, "bin", "myplink", "myplink.sh")

# number of processors
NBPROC = 4

# grm file
GRM = os.path.join(GRM_DIR, 'grm-all-' + str(GRM_CUTOFF), 'all-' + str(GRM_CUTOFF))

# Principal component of grm
PCS = os.path.join(GRM_DIR, 'grm-all-' + str(GRM_CUTOFF), 'all-' + str(GRM_CUTOFF) + '.eigenvec')

HSQ = [0.2, 0.5, 0.8]
N_SNP = [1000, 5000, 10000]  # number of causal SNPs
#N_IND = [100, 200, 400, 800, 1600]
N_IND = [10000]
N_ITER = 5

HSQ_RANGE = (0.5, 0.5)
N_SNP_RANGE = (10000, 10000)
N_IND_RANGE = (50, 9000)

if 'USE_SBATCH' not in globals():
    USE_SBATCH = False
