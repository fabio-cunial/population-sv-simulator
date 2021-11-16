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

parser.add_argument("--bed", "-i", nargs='+', type=str, required=True)
parser.add_argument("--fai", "-f", type=str, required=True)
parser.add_argument("--out_bed", "-o", type=str, required=True)
parser.add_argument("--asfile", "-a", type=str, required=True)
parser.add_argument("--asm_name", "-n", type=str, required=True)
parser.add_argument("--vartype", "-v", type=str, required=True)
parser.add_argument("--svtype", "-s", type=str, required=True)
parser.add_argument("--hap", "-x", type=str, required=True)


args = parser.parse_args()

field_table_file_name = os.path.join(PIPELINE_DIR, 'files/tracks/variant_track_fields.tsv')

# Read variants
df = pd.concat(
    [pd.read_csv(file_name, sep='\t') for file_name in args.bed],
    axis=0
)

df.sort_values(['#CHROM', 'POS'], inplace=True)

# Select table columns
del(df['QUERY_ID'])
del(df['QUERY_POS'])
del(df['QUERY_END'])
del(df['QUERY_STRAND'])

if args.vartype not in {'snv', 'indel'}:
    del(df['SEQ'])

# Read FAI and table columns
df_fai = svpoplib.ref.get_df_fai(args.fai)

# Filter columns that have track annotations
field_set = set(
    pd.read_csv(
        field_table_file_name,
        sep='\t', header=0
    )['FIELD']
)

df = df.loc[:, [col for col in df.columns if col in field_set]]

# Make BigBed
track_name = 'VariantTable'
track_description = f'{args.asm_name} - {args.vartype}-{args.svtype} - {args.hap}'

svpoplib.tracks.variant.make_bb_track(df, df_fai, args.out_bed, args.asfile, track_name, track_description, field_table_file_name)


