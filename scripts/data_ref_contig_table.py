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

parser.add_argument("--ref", "-i", type=str, required=True)
parser.add_argument("--tsv", "-o", type=str, required=True)

args = parser.parse_args()


svpoplib.ref.get_ref_info(
    args.ref
).to_csv(
    args.tsv, sep='\t', index=True, compression='gzip'
)