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

parser.add_argument("--tig_fa", "-f", type=str, required=True)
parser.add_argument("--bed_inv", "-i", type=str, required=True)
parser.add_argument("--bed_density", "-d", type=str, required=True)
parser.add_argument("--fig_dot", "-p", type=str, required=True)
parser.add_argument("--fig_den", "-q", type=str, required=True)
parser.add_argument("--inv_id", "-d", type=str, required=True)
parser.add_argument("--asm_name", "-a", type=str, required=True)
parser.add_argument("--hap", "-x", type=str, required=True)

args = parser.parse_args()


df_inv = pd.read_csv(args.bed_inv, sep='\t', index_col='ID')

if args.inv_id not in df_inv.index:
    raise RuntimeError(
        f'Inversion not found in {args.inv_id} callset '
        f'(assembly={args.asm_name}, hap={args.hap})'
    )

inv_row = df_inv.loc[args.inv_id]

df_density = pd.read_csv(args.bed_density, sep='\t')

# Get inversion call
inv_call = pavlib.inv.get_inv_from_record(inv_row, df_density)

# Make plots
fig_dot = pavlib.plot.dotplot_inv_call(
    inv_call, REF_FA, seq_tig=inv_row['SEQ']
)

fig_density = pavlib.plot.kmer_density_plot(
    inv_call, hap=args.hap, flank_whiskers=True
)

# Write plots
fig_dot.savefig(args.fig_dot, bbox_inches='tight')
fig_density.savefig(args.fig_den, bbox_inches='tight')
