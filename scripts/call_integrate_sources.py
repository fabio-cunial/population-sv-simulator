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

parser.add_argument("--bed_cigar_insdel", "-c", type=str, required=True)
parser.add_argument("--bed_cigar_snv", "-b", type=str, required=True)
parser.add_argument("--bed_lg_ins", "-d", type=str, required=True)
parser.add_argument("--bed_lg_del", "-e", type=str, required=True)
parser.add_argument("--bed_lg_inv", "-f", type=str, required=True)
parser.add_argument("--bed_inv", "-g", type=str, required=True)
parser.add_argument("--bed_ins", "-i", type=str, required=True)
parser.add_argument("--bed_del", "-j", type=str, required=True)
parser.add_argument("--bed_inv_out", "-k", type=str, required=True)
parser.add_argument("--bed_snv", "-l", type=str, required=True)
parser.add_argument("--min_inv", "-m", type=str, required=False, default='300')
parser.add_argument("--max_inv", "-n", type=str, required=False, default='2000000')
parser.add_argument("--tig_filter_pattern", "-t", type=str, required=True, default='None')
parser.add_argument("--asm_name", "-a", type=str, required=True)
parser.add_argument("--hap", "-x", type=str, required=True)

args = parser.parse_args()




# Set parameters
if args.min_inv is not 'None' and args.min_inv != 'unlimited':
    min_inv = int(args.min_inv)
else:
    min_inv = None

if args.max_inv is not 'None' and args.max_inv != 'unlimited':
    max_inv = int(args.max_inv)
else:
    max_inv = None

tig_filter_tree = None

# Read tig filter (if present)


if args.tig_filter_pattern is not 'None':

    tig_filter_file = args.tig_filter_pattern.format(hap=args.hap, asm_name=args.asm_name)

    if os.path.isfile(tig_filter_file):
        tig_filter_tree = collections.defaultdict(intervaltree.IntervalTree)
        df_filter = pd.read_csv(tig_filter_file, sep='\t', header=None, comment='#', usecols=(0, 1, 2))
        df_filter.columns = ['#CHROM', 'POS', 'END']

        for index, row in df_filter.iterrows():
            tig_filter_tree[row['#CHROM']][row['POS']:row['END']] = True 


# Read INV calls
df_inv = pd.concat(
    [
        pd.read_csv(args.bed_inv, sep='\t', low_memory=False),
        pd.read_csv(args.bed_lg_inv, sep='\t', low_memory=False)
    ],
    axis=0
).sort_values(
    ['#CHROM', 'POS']
).reset_index(drop=True)

if min_inv is not None:
    df_inv = df_inv.loc[df_inv['SVLEN'] >= min_inv]

if max_inv is not None:
    df_inv = df_inv.loc[df_inv['SVLEN'] <= max_inv]

# Apply contig filter to INV
df_inv = pavlib.call.filter_by_tig_tree(df_inv, tig_filter_tree)

# Filter overlapping inversion calls
inv_tree = collections.defaultdict(intervaltree.IntervalTree)
inv_index_set = set()

for index, row in df_inv.sort_values(['SVLEN', 'POS']).iterrows():
    if len(inv_tree[row['#CHROM']][row['POS']:row['END']]) == 0:
        inv_index_set.add(index)
        inv_tree[row['#CHROM']][row['POS']:row['END']] = row['ID']

df_inv = df_inv.loc[inv_index_set].sort_values(['#CHROM', 'POS'])

# Initialize filter with inversions
filter_tree = collections.defaultdict(intervaltree.IntervalTree)

for index, row in df_inv.iterrows():
    filter_tree[row['#CHROM']][row['POS']:row['END']] = row['ID']

# Read large variants and filter by inversions
df_lg_ins = pd.read_csv(args.bed_lg_ins, sep='\t', low_memory=False)
df_lg_del = pd.read_csv(args.bed_lg_del, sep='\t', low_memory=False)

if df_lg_ins.shape[0] > 0:
    df_lg_ins = df_lg_ins.loc[
        df_lg_ins.apply(
            lambda row: len(filter_tree[row['#CHROM'][row['POS']:row['END']]]) == 0,
            axis=1
        )
    ]

    # Apply contig filter to large INS
    df_lg_ins = pavlib.call.filter_by_tig_tree(df_lg_ins, tig_filter_tree)

if df_lg_del.shape[0] > 0:
    df_lg_del = df_lg_del.loc[
        df_lg_del.apply(
            lambda row: len(filter_tree[row['#CHROM'][row['POS']:row['END']]]) == 0,
            axis=1
        )
    ]

    # Apply contig filter to large DEL
    df_lg_del = pavlib.call.filter_by_tig_tree(df_lg_del, tig_filter_tree)

    # Add large deletions to filter
    for index, row in df_lg_del.iterrows():
        filter_tree[row['#CHROM']][row['POS']:row['END']] = row['ID']

# Read CIGAR calls
df_cigar_insdel = pd.read_csv(args.bed_cigar_insdel, sep='\t', low_memory=False)
df_snv = pd.read_csv(args.bed_cigar_snv, sep='\t', low_memory=False)

# Check column conformance among INS/DEL callsets (required for merging)
if list(df_cigar_insdel.columns) != list(df_lg_ins.columns):
    raise RuntimeError('Columns from CIGAR and large SV INS callsets do not match')

if list(df_cigar_insdel.columns) != list(df_lg_del.columns):
    raise RuntimeError('Columns from CIGAR and large SV DEL callsets do not match')

# Filter CIGAR calls
if df_cigar_insdel.shape[0] > 0:
    df_cigar_insdel = df_cigar_insdel.loc[
        df_cigar_insdel.apply(
            lambda row: len(filter_tree[row['#CHROM'][row['POS']:row['END']]]) == 0,
            axis=1
        )
    ]

if df_snv.shape[0] > 0:
    df_snv = df_snv.loc[
        df_snv.apply(
            lambda row: len(filter_tree[row['#CHROM'][row['POS']:row['END']]]) == 0,
            axis=1
        )
    ]

# Apply contig filter to small variants
df_cigar_insdel = pavlib.call.filter_by_tig_tree(df_cigar_insdel, tig_filter_tree)
df_snv = pavlib.call.filter_by_tig_tree(df_snv, tig_filter_tree)

# Merge insertion/deletion variants
df_insdel = pd.concat(
    [
        df_lg_ins,
        df_lg_del,
        df_cigar_insdel
    ],
    axis=0
).sort_values(['#CHROM', 'POS'])

# Write
df_insdel.loc[
    df_insdel['SVTYPE'] == 'INS'
].to_csv(
    args.bed_ins, sep='\t', index=False, compression='gzip'
)

df_insdel.loc[
    df_insdel['SVTYPE'] == 'DEL'
].to_csv(
    args.bed_del, sep='\t', index=False, compression='gzip'
)

df_inv.to_csv(args.bed_inv_out, sep='\t', index=False, compression='gzip')
df_snv.to_csv(args.bed_snv, sep='\t', index=False, compression='gzip')
