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
parser.add_argument("--flank_cluster", "-c", type=int, required=False, default=2)
parser.add_argument("--flank_merge", "-d", type=int, required=False, default=2000)
parser.add_argument("--cluster_min_svlen", "-e", type=int, required=False, default=4)
parser.add_argument("--vartype", "-v", type=str, required=True)


args = parser.parse_args()


flank_cluster = args.flank_cluster
flank_merge = args.flank_merge
svlen_min = args.cluster_min_svlen if args.vartype == 'indel' else 50

# Input
df = pd.read_csv(args.bed, sep='\t', header=0, low_memory=False)
df = df.loc[df['SVLEN'] >= svlen_min]

if args.vartype == 'indel':
    df = df.loc[df['SVLEN'] < 50]

# Stop if variants are empty
if df.shape[0] == 0:
    pd.DataFrame(
        [],
        columns=['#CHROM', 'POS', 'END']
    ).to_csv(args.out_bed, sep='\t', index=False, compression='gzip')

else:

    # Subset by SVTYPE
    df_ins = df.loc[df['SVTYPE'] == 'INS']
    df_del = df.loc[df['SVTYPE'] == 'DEL']

    # Load deletion intervals into a tree
    deltree = collections.defaultdict(intervaltree.IntervalTree)

    for index, row in df_del.iterrows():
        deltree[row['#CHROM']][row['POS']:row['END']] = (row['ID'], row['SVLEN'], row['POS'], row['END'])

    # Flag matched INS/DELs
    match_list = list()

    for index, row in df_ins.iterrows():
        flank = row['SVLEN'] * flank_cluster

        match_set = deltree[row['#CHROM']][row['POS'] - flank : row['POS'] + flank]

        if match_set:
            match_list.append(pd.Series(
                [
                    row['#CHROM'],
                    np.min([record.data[2] for record in match_set]),
                    np.max([record.data[3] for record in match_set])
                ],
                index=['#CHROM', 'POS', 'END']
            ))

    if len(match_list) == 0:
        pd.DataFrame(
            [],
            columns=['#CHROM', 'POS', 'END']
        ).to_csv(args.out_bed, sep='\t', index=False, compression='gzip')
    else:

        # Merge overlapping intervals
        df_match = pd.concat(match_list, axis=1).T.sort_values(['#CHROM', 'POS'])

        match_merged_list = list()

        chrom = None
        pos = None
        end = None

        for index, row in df_match.iterrows():

            # Switch chromosomes
            if chrom != row['#CHROM']:

                # Record current record
                if chrom is not None:
                    match_merged_list.append(pd.Series(
                        [chrom, pos, end],
                        index=['#CHROM', 'POS', 'END']
                    ))

                chrom = row['#CHROM']
                pos = row['POS']
                end = row['END']

            # Same chromosome, check record
            if row['POS'] - flank_merge <= end:
                end = np.max([end, row['END']])

            else:
                match_merged_list.append(pd.Series(
                    [chrom, pos, end],
                    index=['#CHROM', 'POS', 'END']
                ))

                pos = row['POS']
                end = row['END']

        df_match = pd.concat(match_merged_list, axis=1).T.sort_values(['#CHROM', 'POS'])

        # Write
        df_match.to_csv(args.out_bed, sep='\t', index=False, compression='gzip')

