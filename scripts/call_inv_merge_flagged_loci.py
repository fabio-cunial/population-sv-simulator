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


def _call_inv_accept_flagged_region(row, allow_single_cluster=False, match_any=set()):
    """
    Annotate which flagged regions are "accepted" (will try to call INV).

    If `allow_single_cluster` is `False` and `match_any` is an empty set, then the only signatures accepted are
    matched SV events or matched indel events.

    :param row: Row from merged flagged regions.
    :param allow_single_cluster: Try to resolve inversions if True for loci with only signatures of clustered SNVs
        and/or clustered indels. Will try many regions and may increase false-positives.
    :param match_any: If defined, contains a set of signatures where at least one must match. Can be used to
        restrict accepted regions to SV-supported signatures only, but inversions with a small uniquely-inverted
        region will likely be missed.

    :return: `True` if the inversion caller should try the region.
    """

    if not allow_single_cluster and (row['TYPE'] == {'CLUSTER_SNV'} or row['TYPE'] == {'CLUSTER_INDEL'}):
        return False

    if match_any and not row['TYPE'] & match_any:
        return False

    return True


parser = argparse.ArgumentParser()

parser.add_argument("--bed_insdel_sv", "-a", type=str, required=True)
parser.add_argument("--bed_insdel_indel", "-b", type=str, required=True)
parser.add_argument("--bed_cluster_indel", "-c", type=str, required=True)
parser.add_argument("--bed_cluster_snv", "-d", type=str, required=True)
parser.add_argument("--bed", "-o", type=str, required=True)
parser.add_argument("--flank", "-f", type=int, required=False, default=500)
parser.add_argument("--batch_count", "-g", type=int, required=False, default=60)
parser.add_argument("--inv_sig_filter", "-i", type=str, required=False, default='svindel')

args = parser.parse_args()


flank = args.flank

# Get region filter parameters
allow_single_cluster = False
match_any = set()

if args.inv_sig_filter is not 'None':
    if args.inv_sig_filter == 'single_cluster':
        allow_single_cluster = True

    elif args.inv_sig_filter == 'svindel':
        match_any.add('MATCH_SV')
        match_any.add('MATCH_INDEL')

    elif args.inv_sig_filter == 'sv':
        match_any.add('MATCH_SV')

    else:
        raise RuntimeError(f'Unrecognized region filter: {args.inv_sig_filter} (must be "single_cluster", "svindel", or "sv")')

# Read
df_indsdel_sv = pd.read_csv(args.bed_insdel_sv, sep='\t')
df_indsdel_indel = pd.read_csv(args.bed_insdel_indel, sep='\t')
df_cluster_indel = pd.read_csv(args.bed_cluster_indel, sep='\t')
df_cluster_snv = pd.read_csv(args.bed_cluster_snv, sep='\t')


# Annotate counts
df_indsdel_sv['COUNT_INDEL'] = 0
df_indsdel_sv['COUNT_SNV'] = 0
df_indsdel_sv['TYPE'] = 'MATCH_SV'

df_indsdel_indel['COUNT_INDEL'] = 0
df_indsdel_indel['COUNT_SNV'] = 0
df_indsdel_indel['TYPE'] = 'MATCH_INDEL'

df_cluster_indel['COUNT_INDEL'] = df_cluster_indel['COUNT']
df_cluster_indel['COUNT_SNV'] = 0
df_cluster_indel['TYPE'] = 'CLUSTER_INDEL'
del(df_cluster_indel['COUNT'])

df_cluster_snv['COUNT_INDEL'] = 0
df_cluster_snv['COUNT_SNV'] = df_cluster_snv['COUNT']
df_cluster_snv['TYPE'] = 'CLUSTER_SNV'
del(df_cluster_snv['COUNT'])

# Merge
df = pd.concat([
    df_indsdel_sv,
    df_indsdel_indel,
    df_cluster_indel,
    df_cluster_snv
], axis=0).sort_values(['#CHROM', 'POS'])

# Merge flagged regions
region_list = list()

chrom = None
pos = 0
end = 0

indel_count = 0
snv_count = 0

type_set = set()

for index, row in df.iterrows():

    if (row['POS'] < end + flank) and (row['#CHROM'] == chrom):

        # Add to existing region
        type_set.add(row['TYPE'])
        end = row['END']

        indel_count += row['COUNT_INDEL']
        snv_count += row['COUNT_SNV']

    else:

        # Write region
        if type_set:
            region_list.append(pd.Series(
                [
                    chrom, pos, end,
                    '{}-{}-RGN-{}'.format(chrom, pos, end - pos),
                    'RGN', end - pos,
                    type_set,
                    #','.join(sorted(type_set)),
                    indel_count, snv_count
                ],
                index=['#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'TYPE', 'COUNT_INDEL', 'COUNT_SNV']
            ))

        # Start new region
        type_set = {row['TYPE']}
        pos = row['POS']
        end = row['END']
        chrom = row['#CHROM']

        indel_count = row['COUNT_INDEL']
        snv_count = row['COUNT_SNV']

# Final region
if type_set:
    region_list.append(pd.Series(
        [
            chrom, pos, end,
            '{}-{}-RGN-{}'.format(chrom, pos, end - pos),
            'RGN', end - pos,
            type_set,
            #','.join(sorted(type_set)),
            indel_count, snv_count
        ],
        index=['#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'TYPE', 'COUNT_INDEL', 'COUNT_SNV']
    ))

# Merge
if len(region_list) > 0:
    df_merged = pd.concat(region_list, axis=1).T.sort_values(['#CHROM', 'POS'])

    # Annotate accepted regions
    df_merged['TRY_INV'] = df_merged.apply(
        _call_inv_accept_flagged_region, allow_single_cluster=allow_single_cluster, match_any=match_any,
        axis=1
    )

    # Group into batches
    df_merged['BATCH'] = -1

    batch = 0

    for index, row in df_merged.iterrows():
        if row['TRY_INV']:
            df_merged.loc[index, 'BATCH'] = batch
            batch = (batch + 1) % args.batch_count

else:
    df_merged = pd.DataFrame([], columns=['#CHROM', 'POS', 'END', 'ID', 'SVTYPE', 'SVLEN', 'TYPE', 'COUNT_INDEL', 'COUNT_SNV', 'TRY_INV', 'BATCH'])

# Write
df_merged.to_csv(args.bed, sep='\t', index=False, compression='gzip')


