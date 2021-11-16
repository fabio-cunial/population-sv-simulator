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

parser.add_argument("--input_insdel", "-i", nargs='+', required=True)
parser.add_argument("--input_snv", "-n", nargs='+', required=True)
parser.add_argument("--output_insdel", "-o", type=str, required=True)
parser.add_argument("--output_snv", "-s", type=str, required=True)

args = parser.parse_args()


pd.concat(
    [pd.read_csv(file_name, sep='\t') for file_name in args.input_insdel],
    axis=0
).sort_values(
    ['#CHROM', 'POS']
).to_csv(
    args.output_insdel, sep='\t', index=False, compression='gzip'
)

pd.concat(
    [pd.read_csv(file_name, sep='\t') for file_name in args.input_snv],
    axis=0
).sort_values(
    ['#CHROM', 'POS']
).to_csv(
    args.output_snv, sep='\t', index=False, compression='gzip'
)
