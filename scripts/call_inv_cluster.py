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


import pavlib
import svpoplib
import kanapy


parser = argparse.ArgumentParser()

parser.add_argument("--bed", "-i", nargs='+', type=str, required=True)
parser.add_argument("--out_bed", "-o", type=str, required=True)
parser.add_argument("--cluster_win", "-w", type=int, required=False, default=200)
parser.add_argument("--cluster_win_min", "-x", type=int, required=False, default=500)
parser.add_argument("--cluster_min_snv", "-y", type=int, required=False, default=20)
parser.add_argument("--cluster_min_indel", "-z", type=int, required=False, default=10)
parser.add_argument("--vartype", "-v", type=str, required=True)

args = parser.parse_args()



# Params
cluster_win = args.cluster_win
cluster_win_min = args.cluster_win

if args.vartype == 'indel':
    cluster_min = args.cluster_min_indel
elif args.vartype == 'snv':
    cluster_min = args.cluster_min_snv
else:
    raise RuntimeError('Bad variant type {}: Expected "indel" or "snv"')

# Read
df = pd.concat(
    [pd.read_csv(input_file_name, sep='\t', usecols=('#CHROM', 'POS', 'END', 'SVTYPE', 'SVLEN')) for input_file_name in args.bed],
    axis=0
).sort_values(['#CHROM', 'POS'])

if args.vartype == 'indel':
    df = df.loc[df['SVLEN'] < 50]

# DEL to midpoint
df['POS'] = (df['END'] + df['POS']) // 2

# Find clusters
cluster_list = list()

chrom = None
cluster_pos = 0
cluster_end = 0
cluster_count = 0

for index, row in df.iterrows():

    if (row['POS'] < cluster_end + cluster_win) and (row['#CHROM'] == chrom):

        # Add to existing cluster
        cluster_count += 1
        cluster_end = row['POS']

    else:

        # Save last cluster
        if (cluster_count >= cluster_min) and (cluster_end - cluster_pos >= cluster_win_min):
            cluster_list.append(pd.Series(
                [chrom, cluster_pos, cluster_end, cluster_count],
                index=['#CHROM', 'POS', 'END', 'COUNT']
            ))

        # Start new cluster
        cluster_count = 1
        cluster_pos = row['POS']
        cluster_end = row['POS']
        chrom = row['#CHROM']

# Save final cluster
if (cluster_count >= cluster_min) and (cluster_end - cluster_pos >= cluster_win_min):
    cluster_list.append(pd.Series(
        [chrom, cluster_pos, cluster_end, cluster_count],
        index=['#CHROM', 'POS', 'END', 'COUNT']
    ))

# Merge records
if len(cluster_list) > 0:
    df_cluster = pd.concat(cluster_list, axis=1).T
else:
    df_cluster = pd.DataFrame([], columns=['#CHROM', 'POS', 'END', 'COUNT'])

# Write
os.makedirs(os.path.dirname(args.out_bed), exist_ok=True)

df_cluster.to_csv(args.out_bed, sep='\t', index=False, compression='gzip')

