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

import pavlib
import svpoplib
import kanapy


def handle_arg_list(arglist):
	if arglist == 'None':
		return None
	else:
		return int(arglist)

parser = argparse.ArgumentParser()

parser.add_argument("--bed_flag", "-f", type=str, required=True)
parser.add_argument("--bed_aln", "-a", type=str, required=True)
parser.add_argument("--tig_fa", "-t", type=str, required=True)
parser.add_argument("--log", "-l", type=str, required=True)
parser.add_argument("--fai", "-i", type=str, required=True)
parser.add_argument("--bed", "-o", type=str, required=True)
parser.add_argument("--k_size", "-k", type=int, required=False, default=31)
parser.add_argument("--inv_threads", "-d", type=int, required=False, default=4)
parser.add_argument("--inv_region_limit", "-r", type=str, required=True)
parser.add_argument("--ref", "-R", type=str, required=True)
parser.add_argument("--inv_min_expand", "-e", type=str, required=True)
parser.add_argument("--batch", "-b", type=int, required=True)
parser.add_argument("--asm_name", "-n", type=str, required=True)
parser.add_argument("--hap", "-p", type=str, required=True)
parser.add_argument("--srs_list", "-s", nargs='*', type=str, required=True)

args = parser.parse_args()

# Get params
batch = int(args.batch)
k_size = args.k_size

density_out_dir = f'results/{args.asm_name}/inv_caller/density_table'

os.makedirs(density_out_dir, exist_ok=True)

if args.srs_list == ['None']:
    srs_list = None
else:
    srs_list = tuple([int(x) for x in args.srs_list])


inv_region_limit = handle_arg_list(args.inv_region_limit)
inv_min_expand = handle_arg_list(args.inv_min_expand)


# Get SRS (state-run-smooth)
srs_tree = pavlib.inv.get_srs_tree(srs_list)  # If none, tree contains a default for all region sizes

# Read and subset table to records in this batch
df_flag = pd.read_csv(args.bed_flag, sep='\t', header=0)
df_flag = df_flag.loc[df_flag['BATCH'] == batch]

if df_flag.shape[0] == 0:
    # No records in batch

    df_bed = pd.DataFrame(
        [],
        columns=[
            '#CHROM', 'POS', 'END',
            'ID', 'SVTYPE', 'SVLEN',
            'HAP',
            'TIG_REGION', 'QUERY_STRAND',
            'CI',
            'RGN_REF_INNER', 'RGN_TIG_INNER',
            'RGN_REF_DISC', 'RGN_TIG_DISC',
            'FLAG_ID', 'FLAG_TYPE',
            'ALIGN_INDEX', 'CLUSTER_MATCH',
            'CALL_SOURCE',
            'SEQ'
        ]
    )

else:

    # Init
    k_util = kanapy.util.kmer.KmerUtil(k_size)

    align_lift = pavlib.align.AlignLift(
        pd.read_csv(args.bed_aln, sep='\t'),
        svpoplib.ref.get_df_fai(args.fai)
    )

    # Read alignment BED
    df_aln = pd.read_csv(
        args.bed_aln,
        sep='\t',
        usecols=('INDEX', 'CLUSTER_MATCH'),
        index_col='INDEX',
        squeeze=False
    )

    # Call inversions
    call_list = list()

    with open(args.log, 'w') as log_file:
        for index, row in df_flag.iterrows():

            # Scan for inversions
            region_flag = pavlib.seq.Region(row['#CHROM'], row['POS'], row['END'])

            try:
                inv_call = pavlib.inv.scan_for_inv(
                    region_flag, args.ref, args.tig_fa, align_lift, k_util,
                    max_region_size=inv_region_limit,
                    threads=args.inv_threads, log=log_file, srs_tree=srs_tree,
                    min_exp_count=inv_min_expand
                )

            except RuntimeError as ex:
                log_file.write('RuntimeError in scan_for_inv(): {}\n'.format(ex))
                inv_call = None

            # Save inversion call
            if inv_call is not None:

                # Get seq
                seq = pavlib.seq.region_seq_fasta(
                    inv_call.region_tig_outer,
                    args.tig_fa,
                    rev_compl=inv_call.region_tig_outer.is_rev
                )

                # Get alignment record data
                aln_index_set = {
                    val for val_list in [
                        inv_call.region_ref_outer.pos_aln_index,
                        inv_call.region_ref_outer.end_aln_index,
                        inv_call.region_ref_inner.pos_aln_index,
                        inv_call.region_ref_inner.end_aln_index

                    ] for val in val_list
                }

                cluster_match_set = {df_aln.loc[index, 'CLUSTER_MATCH'] for index_tuple in aln_index_set for index in index_tuple}

                if False in cluster_match_set:
                    cluster_match = False
                elif np.any(pd.isnull(list(cluster_match_set))):
                    cluster_match = np.nan
                else:
                    if cluster_match_set != {True}:
                        raise RuntimeError(
                            'Found unexpected values in cluster_match: Expected True, False, and np.nan: {}'.format(
                                ', '.join([str(val) for val in cluster_match_set])
                            )
                        )

                    cluster_match = True

                # Save call
                call_list.append(pd.Series(
                    [
                        inv_call.region_ref_outer.chrom,
                        inv_call.region_ref_outer.pos,
                        inv_call.region_ref_outer.end,

                        inv_call.id,
                        'INV',
                        inv_call.svlen,

                        args.hap,

                        inv_call.region_tig_outer.to_base1_string(),
                        '-' if inv_call.region_tig_outer.is_rev else '+',

                        0,

                        inv_call.region_ref_inner.to_base1_string(),
                        inv_call.region_tig_inner.to_base1_string(),

                        inv_call.region_ref_discovery.to_base1_string(),
                        inv_call.region_tig_discovery.to_base1_string(),

                        inv_call.region_flag.region_id(),
                        row['TYPE'],

                        ','.join([str(val) for val in sorted(aln_index_set)]), cluster_match,

                        pavlib.inv.CALL_SOURCE,

                        seq
                    ],
                    index=[
                        '#CHROM', 'POS', 'END',
                        'ID', 'SVTYPE', 'SVLEN',
                        'HAP',
                        'TIG_REGION', 'QUERY_STRAND',
                        'CI',
                        'RGN_REF_INNER', 'RGN_TIG_INNER',
                        'RGN_REF_DISC', 'RGN_TIG_DISC',
                        'FLAG_ID', 'FLAG_TYPE',
                        'ALIGN_INDEX', 'CLUSTER_MATCH',
                        'CALL_SOURCE',
                        'SEQ'
                    ]
                ))

                # Save density table
                inv_call.df.to_csv(
                    os.path.join(density_out_dir, 'density_{}_{}.tsv.gz'.format(inv_call.id, args.hap)),
                    sep='\t', index=False, compression='gzip'
                )

                # Call garbage collector
                gc.collect()

    # Merge records
    if len(call_list) > 0:
        df_bed = pd.concat(call_list, axis=1).T

    else:
        # Create emtpy data frame
        df_bed = pd.DataFrame(
            [],
            columns=[
                '#CHROM', 'POS', 'END',
                'ID', 'SVTYPE', 'SVLEN',
                'HAP',
                'TIG_REGION', 'QUERY_STRAND',
                'CI',
                'RGN_REF_INNER', 'RGN_TIG_INNER',
                'RGN_REF_DISC', 'RGN_TIG_DISC',
                'FLAG_ID', 'FLAG_TYPE',
                'ALIGN_INDEX', 'CLUSTER_MATCH',
                'CALL_SOURCE',
                'SEQ'
            ]
        )

# Write
df_bed.to_csv(args.bed, sep='\t', index=False, compression='gzip')
