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

parser.add_argument("--bed", "-e", type=str, required=True)
parser.add_argument("--tsv_group", "-t", type=str, required=True)
parser.add_argument("--fa", "-f", type=str, required=True)
parser.add_argument("--fai", "-i", type=str, required=True)
parser.add_argument("--bed_n", "-n", type=str, required=True)
parser.add_argument("--bed_ins", "-o", type=str, required=True)
parser.add_argument("--bed_del", "-d", type=str, required=True)
parser.add_argument("--bed_inv", "-v", type=str, required=True)
parser.add_argument("--k_size", "-k", type=int, required=True, default=31)
parser.add_argument("--inv_threads_lg", "-l", type=int, required=False, default=12)
parser.add_argument("--inv_region_limit", "-m", type=str, required=False, default=None)
parser.add_argument("--srs_list", "-s", nargs='*', type=str, required=True)
parser.add_argument("--hap", "-x", type=str, required=True)
parser.add_argument("--ref", "-r", type=str, required=True)
parser.add_argument("--log", "-z", type=str, required=True)
parser.add_argument("--batch", "-b", type=int, required=True)
parser.add_argument("--asm_name", "-a", type=str, required=True)


args = parser.parse_args()

if args.srs_list == ['None']:
	srs_list = None
else:
	srs_list = tuple([int(x) for x in args.srs_list])

if args.inv_region_limit == 'None':
    args.inv_region_limit = None
else:
    args.inv_region_limit = int(args.inv_region_limit)



srs_tree = pavlib.inv.get_srs_tree(srs_list)  # If none, tree contains a default for all region sizes

# Read
df = pd.read_csv(args.bed, sep='\t')
df_tig_fai = svpoplib.ref.get_df_fai(args.fai)

# Subset to alignment records in this batch
df_group = pd.read_csv(args.tsv_group, sep='\t')
df_group = df_group.loc[df_group['BATCH'] == int(args.batch)]

group_set = set(df_group[['CHROM', 'TIG']].apply(tuple, axis=1))

if df.shape[0] > 0:
    df = df.loc[df.apply(lambda row: (row['#CHROM'], row['QUERY_ID']) in group_set, axis=1)]

# Get trees of N bases
n_tree = collections.defaultdict(intervaltree.IntervalTree)

df_n = pd.read_csv(args.bed_n, sep='\t')

for index, row in df_n.iterrows():
    n_tree[row['#CHROM']][row['POS']:row['END']] = True

# Make density table output directory
density_out_dir = f'results/{args.asm_name}/inv_caller/density_table_lg'
os.makedirs(density_out_dir, exist_ok=True)

# Get large events
with open(args.log, 'w') as log_file:
    df_ins, df_del, df_inv = pavlib.lgsv.scan_for_events(
        df, df_tig_fai, args.hap, args.ref, args.fa,
        k_size=args.k_size,
        n_tree=n_tree,
        srs_tree=srs_tree,
        threads=args.inv_threads_lg,
        log=log_file,
        density_out_dir=density_out_dir,
        max_region_size=args.inv_region_limit
    )

# Write
df_ins.to_csv(args.bed_ins, sep='\t', index=False, compression='gzip')
df_del.to_csv(args.bed_del, sep='\t', index=False, compression='gzip')
df_inv.to_csv(args.bed_inv, sep='\t', index=False, compression='gzip')


