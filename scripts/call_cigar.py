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

parser.add_argument("--input_bed", "-i", type=str, required=True)
parser.add_argument("--input_tig_fa", "-f", type=str, required=True)
parser.add_argument("--batch", "-b", type=int, required=True)
parser.add_argument("--out_insdel", "-o", type=str, required=True)
parser.add_argument("--out_snv", "-s", type=str, required=True)
parser.add_argument("--ref_fa", "-r", type=str, required=True)
parser.add_argument("--hap", "-x", type=str, required=True)

args = parser.parse_args()



batch = int(args.batch)

os.makedirs(os.path.dirname(args.out_insdel), exist_ok=True)  # Random crashes with "FileNotFoundError", Snakemake not creating output directory?

# Read
df_align = pd.read_csv(args.input_bed, sep='\t', dtype={'#CHROM': str})

df_align = df_align.loc[df_align['CALL_BATCH'] == batch]

# Call
df_snv, df_insdel = pavlib.cigarcall.make_insdel_snv_calls(df_align, args.ref_fa, args.input_tig_fa, args.hap)

# Write
df_insdel.to_csv(args.out_insdel, sep='\t', index=False, compression='gzip')
df_snv.to_csv(args.out_snv, sep='\t', index=False, compression='gzip')