#!/bin/env python


import collections
import gc
import gzip
import intervaltree
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pysam
import subprocess
import sys
import re
from scipy import stats
import shutil
import tempfile

from Bio import SeqIO
import Bio.bgzf

import argparse
PIPELINE_DIR = os.path.dirname(os.path.abspath(__file__)).replace('scripts', '')
sys.path.append(PIPELINE_DIR)  # pavlib
sys.path.append(os.path.join(PIPELINE_DIR, 'dep', 'svpop'))  # svpoplib
sys.path.append(os.path.join(PIPELINE_DIR, 'dep'))  # kanapy

import pavlib
import svpoplib
import kanapy



parser = argparse.ArgumentParser()

parser.add_argument("--bed", "-b", type=str, required=True)
parser.add_argument("--min_trim_tig_len", "-m", type=int, required=False, default=1000)
parser.add_argument("--tig_fai", "-t", type=str, required=True)
parser.add_argument("--output", "-o", type=str, required=True)



args = parser.parse_args()


df = pavlib.align.trim_alignments(
    pd.read_csv(args.bed, sep='\t'),  # Untrimmed alignment BED
    args.min_trim_tig_len,  # Minimum contig length
    args.tig_fai  # Path to alignment FASTA FAI
)

# Add batch ID for CIGAR calling (calls in batches)
df['CALL_BATCH'] = df['INDEX'].apply(lambda val: val % pavlib.cigarcall.CALL_CIGAR_BATCH_COUNT)

# Write
df.to_csv(args.output, sep='\t', index=False, compression='gzip')
