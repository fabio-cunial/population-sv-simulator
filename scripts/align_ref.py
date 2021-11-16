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

parser.add_argument("--ref_fa", "-i", type=str, required=True)
parser.add_argument("--out_ref_fa", "-o", type=str, required=True)
parser.add_argument("--ref_fai", "-x", type=str, required=True)

args = parser.parse_args()



# Copy FASTA to FA/GZ
pavlib.seq.copy_fa_to_gz(args.ref_fa, args.out_ref_fa)

# Index
if os.stat(args.out_ref_fa).st_size > 0:
    os.system(f"""samtools faidx {args.out_ref_fa}""")

else:
    with open(args.ref_fai, 'w') as out_file:
        pass