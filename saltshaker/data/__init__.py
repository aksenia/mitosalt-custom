"""
Default data files for SaltShaker
"""
import os

# Path to package data directory
DATA_DIR = os.path.dirname(os.path.abspath(__file__))

# Default gene annotations for human mitochondrial genome (hg38/rCRS)
DEFAULT_MT_GENES = os.path.join(DATA_DIR, 'gencode.v49.annotation.MT_genes.bed')

# Default blacklist regions for human mitochondrial genome
DEFAULT_MT_BLACKLIST = os.path.join(DATA_DIR, 'mt_blacklist_regions.bed')