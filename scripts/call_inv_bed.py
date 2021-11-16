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

parser.add_argument("--bed", "-i", type=str, required=True)
parser.add_argument("--out_bed", "-o", type=str, required=True)
parser.add_argument("--dropped", "-d", dtype=str, required=True)
parser.add_argument("--min_svlen", "-m", type=int, required=True)

args = parser.parse_args()



        # Read inversions
df = pd.read_csv(args.bed, sep='\t')

# Filter
df_drop = df.loc[df['SVLEN'] < args.min_svlen]
df_drop.to_csv(args.dropped, sep='\t', index=False, compression='gzip')

df = df.loc[[index for index in df.index if index not in df_drop.index]]

df.drop_duplicates('ID', inplace=True)

# Write BED
df.to_csv(args.out_bed, sep='\t', index=False, compression='gzip')