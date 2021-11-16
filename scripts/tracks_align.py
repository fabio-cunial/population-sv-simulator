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
parser.add_argument("--asfile", "-a", type=str, required=True)
parser.add_argument("--align_stage", "-l", type=str, required=True)

args = parser.parse_args()


if args.align_stage == 'pre-cut':
    track_desc_short = 'PreTrim'
    track_description = 'Pre-trimmed alignments'

elif args.align_stage == 'post-cut':
    track_desc_short = 'PostTrim'
    track_description = 'Post-trimmed alignments'

else:
    raise RuntimeError(f'Unknown align_stage wildcard: {args.align_stage}')

# Read field table
df_as = pd.read_csv(
    os.path.join(PIPELINE_DIR, 'files/tracks/alignment_track_fields.tsv'),
    sep='\t'
).set_index('FIELD')

# Read variants
df = pd.concat(
    [pd.read_csv(file_name, sep='\t') for file_name in args.bed],
    axis=0
).sort_values(['#CHROM', 'POS', 'END', 'QUERY_ID'])

# Add BED fields
df['POS_THICK'] = df['POS']
df['END_THICK'] = df['END']
df['ID'] = df['QUERY_ID']
df['SCORE'] = 1000
df['COL'] = df['HAP'].apply(lambda val: ALIGN_COLOR.get(val, '0,0,0'))
df['STRAND'] = df['REV'].apply(lambda val: '-' if val else '+')

del(df['CIGAR'])

# Sort columns
head_cols = ['#CHROM', 'POS', 'END', 'ID', 'SCORE', 'STRAND', 'POS_THICK', 'END_THICK', 'COL']
tail_cols = [col for col in df.columns if col not in head_cols]

df = df[head_cols + tail_cols]

# Check AS fields
missing_fields = [col for col in df.columns if col not in df_as.index]

if missing_fields:
    raise RuntimeError('Missing {} fields in AS definition: {}{}'.format(
        len(missing_fields), ', '.join(missing_fields[:3]), '...' if len(missing_fields) else ''
    ))

# Write AS file
with open(args.asfile, 'w') as out_file:

    # Heading
    out_file.write('table Align{}\n"{}"\n(\n'.format(track_desc_short, track_description))

    # Column definitions
    for col in df.columns:
        out_file.write('{TYPE} {NAME}; "{DESC}"\n'.format(**df_as.loc[col]))

    # Closing
    out_file.write(')\n')

# Write BED
df.to_csv(args.out_bed, sep='\t', index=False)
