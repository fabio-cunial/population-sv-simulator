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

parser.add_argument("--sam", "-s", type=str, required=True)
parser.add_argument("--align_head", "-a", type=str, required=True)
parser.add_argument("--bed_out", "-o", type=str, required=True)
parser.add_argument("--tig_fai", "-t", type=str, required=True)
parser.add_argument("--chrom_cluster", "-c", type=bool, required=False, default=False)
parser.add_argument("--hap", "-x", type=str, required=True)


args = parser.parse_args()


if os.stat(args.sam).st_size == 0:
    pd.DataFrame(
        [],
        columns=[
            '#CHROM', 'POS', 'END',
            'INDEX',
            'QUERY_ID', 'QUERY_POS', 'QUERY_END',
            'QUERY_TIG_POS', 'QUERY_TIG_END',
            'RG', 'AO',
            'MAPQ',
            'REV', 'FLAGS', 'HAP',
            'CIGAR'
        ]
    ).to_csv(
        args.bed_out, sep='\t', index=False, compression='gzip'
    )

    with open(args.align_head, 'w') as out_file:
        pass

else:

    # Read FAI
    df_tig_fai = svpoplib.ref.get_df_fai(args.tig_fai)
    df_tig_fai.index = df_tig_fai.index.astype(str)

    # Read alignments as a BED file.
    df = pavlib.align.get_align_bed(args.sam, df_tig_fai, args.hap, args.chrom_cluster)

    # Write SAM headers
    with gzip.open(args.sam, 'rt') as in_file:
        with gzip.open(args.align_head, 'wt') as out_file:

            line = next(in_file)

            while True:

                if not line.strip():
                    continue

                if not line.startswith('@'):
                    break

                out_file.write(line)

                try:
                    line = next(in_file)
                except StopIteration:
                    break

    # Write
    df.to_csv(args.bed_out, sep='\t', index=False, compression='gzip')
