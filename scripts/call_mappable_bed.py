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

parser.add_argument("--bed_align", "-a", type=str, required=True)
parser.add_argument("--bed_lg_del", "-d", type=str, required=True)
parser.add_argument("--bed_lg_ins", "-i", type=str, required=True)
parser.add_argument("--bed_lg_inv", "-n", type=str, required=True)
parser.add_argument("--flank", "-f", type=str, required=True)
parser.add_argument("--bed", "-o", type=str, required=True)

args = parser.parse_args()


try:
    flank = np.int32(args.flank)

except ValueError:
    raise RuntimeError(f'Flank parameter is not an integer: {args.flank}')

# Merge
df = pavlib.util.region_merge(
    [
        args.bed_align,
        args.bed_lg_del,
        args.bed_lg_ins,
        args.bed_lg_inv
    ],
    pad=flank
)

# Write
df.to_csv(args.bed, sep='\t', index=False, compression='gzip')