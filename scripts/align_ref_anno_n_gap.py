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
import argparse

from Bio import SeqIO
import Bio.bgzf
import argparse


import pavlib
import svpoplib
import kanapy


parser = argparse.ArgumentParser()

parser.add_argument("--ref_fa", "-i", type=str, required=True)
parser.add_argument("--bed", "-o", type=str, required=True)

args = parser.parse_args()


with gzip.open(args.bed, 'wt') as out_file:
    out_file.write('#CHROM\tPOS\tEND\n')

    with gzip.open(args.ref_fa, 'rt') as in_file:
        for record in Bio.SeqIO.parse(in_file, 'fasta'):

            pos = None
            end = None

            enum_list = [i for i, val in enumerate(str(record.seq).upper()) if val == 'N']

            for index in enum_list:

                if pos is None:
                    pos = end = index

                elif index == end + 1:
                    end = index

                else:
                    out_file.write(f'{record.id}\t{pos}\t{end + 1}\n')
                    pos = end = index

            if pos is not None:
                out_file.write(f'{record.id}\t{pos}\t{end + 1}\n')