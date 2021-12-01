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

parser.add_argument("--bed_ins", "-i", type=str, required=True)
parser.add_argument("--bed_del", "-d", type=str, required=True)
parser.add_argument("--bed_inv", "-n", type=str, required=True)
parser.add_argument("--bed_snv", "-s", type=str, required=True)
parser.add_argument("--bed_snv_snv", "-o", type=str, required=True)
parser.add_argument("--bed_indel_ins", "-p", type=str, required=True)
parser.add_argument("--bed_indel_del", "-q", type=str, required=True)
parser.add_argument("--bed_sv_ins", "-r", type=str, required=True)
parser.add_argument("--bed_sv_del", "-z", type=str, required=True)
parser.add_argument("--bed_sv_inv", "-t", type=str, required=True)
parser.add_argument("--fa_indel_ins", "-u", type=str, required=True)
parser.add_argument("--fa_indel_del", "-v", type=str, required=True)
parser.add_argument("--fa_sv_ins", "-w", type=str, required=True)
parser.add_argument("--fa_sv_del", "-x", type=str, required=True)
parser.add_argument("--fa_sv_inv", "-y", type=str, required=True)

args = parser.parse_args()


df_ins = pd.read_csv(args.bed_ins, sep='\t', low_memory=False)
df_del = pd.read_csv(args.bed_del, sep='\t', low_memory=False)
df_inv = pd.read_csv(args.bed_inv, sep='\t', low_memory=False)

df_svtype_dict = {
    'ins': df_ins,
    'del': df_del,
    'inv': df_inv
}

for vartype, svtype in [('sv', 'ins'), ('sv', 'del'), ('sv', 'inv'), ('indel', 'ins'), ('indel', 'del')]:

    df = df_svtype_dict[svtype]

    # Subset
    if vartype == 'sv':
        df = df.loc[df['SVLEN'] >= 50]

    elif vartype == 'indel':
        df = df.loc[df['SVLEN'] < 50]

    else:
        raise RuntimeError('Program Bug: Unknown variant type: {}'.format(vartype))

    # Get output file names
    bed_file_name = vars(args)[f'bed_{vartype}_{svtype}']
    fa_file_name = vars(args)[f'fa_{vartype}_{svtype}']

    # Write FASTA
    with Bio.bgzf.BgzfWriter(fa_file_name, 'wb') as out_file:
        SeqIO.write(svpoplib.seq.bed_to_seqrecord_iter(df), out_file, 'fasta')

    del(df['SEQ'])

    # Arrange columns
    df = svpoplib.variant.order_variant_columns(
        df,
        head_cols=[
            '#CHROM', 'POS', 'END', 'ID',
            'SVTYPE', 'SVLEN', 'REF', 'ALT',
            'HAP', 'GT', 'CLUSTER_MATCH', 'CALL_SOURCE',
            'TIG_REGION', 'QUERY_STRAND', 'ALIGN_INDEX'
            'CI',
        ],
        tail_cols=[
            'CALL_SOURCE',
            'HAP_VARIANTS',
            'CI',
            'HAP_RO', 'HAP_OFFSET', 'HAP_SZRO', 'HAP_OFFSZ'
        ],
        allow_missing=True
    )

    # Write BED
    df.to_csv(bed_file_name, sep='\t', index=False, compression='gzip')

# Process SNVs
df = pd.read_csv(args.bed_snv, sep='\t', low_memory=False)
df.to_csv(args.bed_snv_snv, sep='\t', index=False, compression='gzip')




