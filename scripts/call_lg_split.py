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

parser.add_argument("--bed", "-i", type=str, required=True)
parser.add_argument("--tsv", "-o", type=str, required=True)
parser.add_argument("--batch_count", "-b", type=int, required=False, default=10)

args = parser.parse_args()


df = pd.read_csv(args.bed, sep='\t')

# Get ref/tig pairs with multiple mappings
tig_map_count = collections.Counter(df[['#CHROM', 'QUERY_ID']].apply(tuple, axis=1))

df_group_list = list()

index = 0

for chrom, tig in [(chrom, tig_id) for (chrom, tig_id), count in tig_map_count.items() if count > 1]:
    df_group_list.append(pd.Series(
        [chrom, tig, index % args.batch_count],
        index=['CHROM', 'TIG', 'BATCH']
    ))

    index += 1

# Merge (CHROM, TIG, BATCH)
if len(df_group_list) > 0:
    df_group = pd.concat(df_group_list, axis=1).T
else:
    df_group = pd.DataFrame([], columns=['CHROM', 'TIG', 'BATCH'])

# Write
df_group.to_csv(args.tsv, sep='\t', index=False, compression='gzip')

