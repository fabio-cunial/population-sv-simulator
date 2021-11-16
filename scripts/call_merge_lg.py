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
parser.add_argument("--out_bed", "-o", type=str, required=True)

args = parser.parse_args()


pd.concat(
    [pd.read_csv(file_name, sep='\t') for file_name in args.bed],
    axis=0
).to_csv(
    args.out_bed, sep='\t', index=False, compression='gzip'
)