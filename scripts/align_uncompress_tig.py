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

parser.add_argument("--input", "-i", type=str, required=True)
parser.add_argument("--output", "-o", type=str, required=False, default='/dev/stdout')

args = parser.parse_args()


if os.stat(args.input).st_size > 0:
    os.system(
        f"""zcat {args.input} > {args.output}"""
    )
else:
    with open(args.output, 'w') as out_file:
        pass
