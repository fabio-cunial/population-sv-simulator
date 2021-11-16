"""
Data files including reference and data tables for the reference.
"""


#
# Reference contig data table
#

# data_ref_contig_table
#
# Contig table.
rule data_ref_contig_table:
    input:
        ref = config['reference']
    output:
        tsv='data/ref/contig_info.tsv.gz'
    shell:
        '''
        python3 {PIPELINE_DIR}/scripts/data_ref_contig_table.py -i {input.ref} -o {output.tsv}
        '''



#
# Reference annotation
#

# align_ref_anno_n_gap
#
# Find locations of N-gaps.
rule data_align_ref_anno_n_gap:
    input:
        ref_fa='data/ref/ref.fa.gz'
    output:
        bed='data/ref/n_gap.bed.gz'
    shell:
        '''
        python3 {PIPELINE_DIR}/scripts/align_ref_anno_n_gap.py -i {input.ref_fa} -o {output.bed}
        '''

#
# Reference
#

# align_ref_lra_index
#
# Index reference for LRA.
rule data_align_ref_lra_index:
    input:
        fa='data/ref/ref.fa.gz',
    output:
        gli='data/ref/ref.fa.gz.gli',
        mmi='data/ref/ref.fa.gz.mmi'
    shell:
        """lra index -CONTIG {input.fa}"""

# align_ref
#
# Setup reference.
rule data_align_ref:
    input:
        ref_fa=config['reference']
    output:
        ref_fa='data/ref/ref.fa.gz',
        ref_fai='data/ref/ref.fa.gz.fai'
    shell:
        '''
        python3 {PIPELINE_DIR}/scripts/align_ref.py -i {input.ref_fa} -o {output.ref_fa} -x {output.ref_fai}
        '''