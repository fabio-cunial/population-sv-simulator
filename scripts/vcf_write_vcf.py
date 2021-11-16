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

parser.add_argument("--bed_snv_snv", "-a", type=str, required=True)
parser.add_argument("--bed_indel_ins", "-b", type=str, required=True)
parser.add_argument("--bed_indel_del", "-c", type=str, required=True)
parser.add_argument("--bed_sv_ins", "-d", type=str, required=True)
parser.add_argument("--bed_sv_del", "-e", type=str, required=True)
parser.add_argument("--bed_sv_inv", "-f", type=str, required=True)
parser.add_argument("--fa_indel_ins", "-g", type=str, required=True)
parser.add_argument("--fa_indel_del", "-i", type=str, required=True)
parser.add_argument("--fa_sv_ins", "-j", type=str, required=True)
parser.add_argument("--fa_sv_del", "-k", type=str, required=True)
parser.add_argument("--fa_sv_inv", "-l", type=str, required=True)
parser.add_argument("--ref_tsv", "-m", type=str, required=True)
parser.add_argument("--vcf", "-o", type=str, required=True)
parser.add_argument("--asm_name", "-n", type=str, required=True)
parser.add_argument("--reference", "-r", type=str, required=True)


args = parser.parse_args()


if args.asm_name in {'#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'}:
    raise RuntimeError(f'Assembly name conflicts with a VCF header column name: {args.asm_name}')


_VCF_SVTYPE = {
    'sv': ('ins', 'del', 'inv'),
    'indel': ('ins', 'del'),
    'snv': ('snv', )
}

# Check alt format
symbolic_alt = False
# if args.alt_fmt == 'alt':
#     symbolic_alt = False
# elif args.alt_fmt == 'sym':
#     symbolic_alt = True
# else:
#     raise RuntimeError(f'Unknown alt format wildcard (alt_fmt): {alt_fmt}')

# Process variant types
df_list = list()

for vartype in ('sv', 'indel', 'snv'):
    for svtype in _VCF_SVTYPE[vartype]:
        print(f'{vartype} - {svtype}')  # DBGTMP

        # Read variants
        bed_file_name = vars(args)[f'bed_{vartype}_{svtype}']
        df = pd.read_csv(bed_file_name, sep='\t')

        df['VARTYPE'] = vartype.upper()

        # Read sequence from FASTA
        if vartype in ('sv', 'indel'):

            df.set_index('ID', inplace=True)

            # Get Sequence from FASTA and assign to SEQ column (not for SNVs)
            fa_file_name = vars(args)[f'fa_{vartype}_{svtype}']
            with gzip.open(fa_file_name, 'rt') as fa_in:
                df_seq_dict = {
                    record.name: str(record.seq) for record in SeqIO.parse(fa_in, 'fasta')
                }

            df['SEQ'] = pd.Series(df_seq_dict)

            del(df_seq_dict)

            df.reset_index(inplace=True)

            # Check for missing sequences
            df_null_seq = df.loc[pd.isnull(df['SEQ'])]

            if df_null_seq.shape[0] > 0:
                id_list = ', '.join(df_null_seq.iloc[:3]['ID'])

                raise RuntimeError(
                    'Missing FASTA sequence for {} variants (vartype={}, svtype={}): {}{}'.format(
                        df_null_seq.shape[0],
                        vartype,
                        svtype,
                        ', '.join([str(val) for val in df_null_seq.iloc[:3]['ID']]),
                        '...' if df_null_seq.shape[0] > 3 else ''
                    )
                )

        # Reformat fields for INFO
        for col in ('CALL_SOURCE', 'TIG_REGION', 'QUERY_STRAND', 'RGN_REF_INNER', 'RGN_TIG_INNER'):
            if col in df.columns:
                df[col] = df[col].apply(lambda val: val.replace(';', ','))

        if svtype == 'del':
            df['SVLEN'] = - np.abs(df['SVLEN'])

        # INFO: Base
        df['INFO'] = df.apply(lambda row: 'ID={ID};SVTYPE={SVTYPE}'.format(**row), axis=1)

        # INFO: Add SV/INDEL annotations
        if vartype != 'snv':
            df['INFO'] = df.apply(lambda row: row['INFO'] + ';SVLEN={SVLEN}'.format(**row), axis=1)

        # INFO: Add contig placement info
        df['INFO'] = df.apply(lambda row: row['INFO'] + ';TIG_REGION={TIG_REGION};QUERY_STRAND={QUERY_STRAND}'.format(**row), axis=1)

        # INFO: Add INV
        if svtype == 'inv':
            df['INFO'] = df.apply(lambda row: row['INFO'] + ';INNER_REF={RGN_REF_INNER};INNER_TIG={RGN_TIG_INNER}'.format(**row), axis=1)

        # INFO: Add breakpoint homology
        if svtype in {'ins', 'del'}:
            df['INFO'] = df.apply(lambda row: row['INFO'] + ';HOM_REF={HOM_REF};HOM_TIG={HOM_TIG}'.format(**row), axis=1)

        # REF
        if 'REF' not in df.columns:
            df['REF'] = list(svpoplib.vcf.ref_base(df, args.reference))

        if svtype == 'inv' and not symbolic_alt:
            df_ref_base = df['REF']
            df['REF'] = df_ref_base + svpoplib.ref.get_ref_region(df, args.reference).apply(lambda val: val.upper())

        # ALT
        if vartype != 'snv':
            if symbolic_alt:
                df['ALT'] = df['SVTYPE'].apply(lambda val: f'<{val}>')

                df['INFO'] = df.apply(lambda row: row['INFO'] + ';SEQ={SEQ}'.format(**row))

            else:

                # Check for sequence types that cannot be placed in ALT (may need symbolic ALTs)
                if svtype != 'inv':
                    df['REF'] = df.apply(lambda row:
                        ((row['REF'] + row['SEQ']) if row['POS'] > 0 else (row['SEQ'] + row['REF'])) if row['SVTYPE'] == 'DEL' else row['REF'],
                        axis=1
                    )

                    df['ALT'] = df.apply(lambda row:
                        ((row['REF'] + row['SEQ']) if row['POS'] > 0 else (row['SEQ'] + row['REF'])) if row['SVTYPE'] == 'INS' else row['REF'][0],
                        axis=1
                    )

                else:
                    df['ALT'] = df_ref_base + df['SEQ']

                del df['SEQ']

                df['ALT'] = df['ALT'].apply(lambda val: val.upper())

        else:
            # Fix position for SNVs (0-based BED to 1-based VCF)
            df['POS'] += 1

            df['ALT'] = df['ALT'].apply(lambda val: val.upper())

        # Save columns needed for VCF
        df = df[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'INFO', 'GT']]

        df_list.append(df)


# Merge
df = pd.concat(df_list, axis=0)
df.sort_values(['#CHROM', 'POS'], inplace=True)

# INFO headers
info_header_list = list()

info_header_list.append(('ID', '1', 'String', 'Variant ID'))
info_header_list.append(('SVTYPE', '1', 'String', 'Variant type'))
info_header_list.append(('SVLEN', '.', 'String', 'Variant length'))
info_header_list.append(('TIG_REGION', '.', 'String', 'Contig region where variant was found (one per alt with h1 before h2 for homozygous calls)'))
info_header_list.append(('QUERY_STRAND', '.', 'String', 'Strand of variant in the contig relative to the reference (order follows TIG_REGION)'))
info_header_list.append(('INNER_REF', '.', 'String', 'Inversion inner breakpoint in reference coordinates (order follows TIG_REGION)'))
info_header_list.append(('INNER_TIG', '.', 'String', 'Inversion inner breakpoint in contig coordinates (order follows TIG_REGION)'))
info_header_list.append(('HOM_REF', '.', 'String', 'Perfect breakpoint homology (SV sequence vs reference). Format \'X,Y\' where X homology upstream, and Y is homology downstream. Homology vs reference is often better for DEL.'))
info_header_list.append(('HOM_TIG', '.', 'String', 'Perfect breakpoint homology (SV sequence vs contig). Format \'X,Y\' where X homology upstream, and Y is homology downstream. Homology vs contig is often better for INS.'))

if symbolic_alt:
    info_header_list.append(('SEQ', '.', 'String', 'SV or indel sequence'))

# ALT headers
alt_header_list = list()

if symbolic_alt:
    alt_header_list.append(('INS', 'Sequence insertion'))
    alt_header_list.append(('DEL', 'Sequence deletion'))
    alt_header_list.append(('INV', 'Inversion'))

# QUAL, FILTER, FORMAT
if 'QUAL' not in df.columns:
    df['QUAL'] = '.'

if 'FILTER' in df.columns:
    raise RuntimeError('FILTER is defined in dataframe, but FILTER headers are not yet implemented')

filter_header_list = list()
df['FILTER'] = '.'

df['FORMAT'] = 'GT'

format_header_list = [
    ('GT', '1', 'String', 'Genotype')
]

# VCF order
df = df[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'GT']]

df.columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', args.asm_name]

# Write
df_ref = pd.read_csv(args.ref_tsv, sep='\t')

with Bio.bgzf.open(args.vcf, 'wt') as out_file:
    for line in svpoplib.vcf.header_list(
        df_ref,
        info_header_list,
        format_header_list,
        alt_header_list,
        filter_header_list,
        variant_source='PAV {}'.format(pavlib.constants.get_version_string()),
        ref_file_name=os.path.basename(args.reference)
    ):
        out_file.write(line)

    out_file.write('\t'.join(df.columns))
    out_file.write('\n')

    for index, row in df.iterrows():
        out_file.write('\t'.join(row.astype(str)))
        out_file.write('\n')

# Write tabix index if possible
try:
    os.system(f"""tabix {args.vcf} && touch -r {args.vcf} {args.vcf}.tbi""")
except:
    pass


