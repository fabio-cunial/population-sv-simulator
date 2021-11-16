"""
Call inversions from aligned contigs.

Inversion calling has two key steps:
1) Flag: Find signatures of inversions from alignments and/or variant calls.
2) Call: Use flagged loci to call inversions. Calling starts with the flagged locus, then it expands until it
   finds the full inversion including unique reference sequence on each flank.
"""

###################
### Definitions ###
###################


def _input_call_inv_cluster(wildcards):
    """
    Get input for flagging by clustered variant calls. Return both ins & del indels or snvs.

    :param wildcards: Wildcards.

    :return: A list of input files.
    """

    if wildcards.vartype == 'indel':
        return [
            'temp/{asm_name}/cigar/pre_inv/svindel_insdel_{hap}.bed.gz'.format(**wildcards),
        ]

    elif wildcards.vartype == 'snv':
        return [
            'temp/{asm_name}/cigar/pre_inv/snv_snv_{hap}.bed.gz'.format(**wildcards),
        ]


    raise RuntimeError('Unknown variant type {} for input into call_inv_cluster: Expected "indel" or "snv"'.format(
        wildcards.vartype
    ))

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

BATCH_COUNT_DEFAULT = 60


#############
### Rules ###
#############


#
# Call inversions
#

# call_inv_batch_merge
#
# Merge batches.
rule call_inv_batch_merge:
    input:
        bed=expand('temp/{{asm_name}}/inv_caller/batch/{{hap}}/inv_call_{batch}.bed.gz', batch=range(int(config.get('inv_sig_batch_count', BATCH_COUNT_DEFAULT))))
    output:
        bed=temp('temp/{asm_name}/inv_caller/sv_inv_{hap}.bed.gz')
    shell:
        '''
        python3 {PIPELINE_DIR}/scripts/call_inv_batch_merge.py -i {input.bed} -o {output.bed}
        '''

# call_inv_batch
#
# Call inversions in batches of flagged regions.
rule call_inv_batch:
    input:
        bed_flag='results/{asm_name}/inv_caller/flagged_regions_{hap}.bed.gz',
        bed_aln='results/{asm_name}/align/aligned_tig_{hap}.bed.gz',
        tig_fa='temp/{asm_name}/align/contigs_{hap}.fa.gz',
        fai='temp/{asm_name}/align/contigs_{hap}.fa.gz.fai'
    output:
        bed=temp('temp/{asm_name}/inv_caller/batch/{hap}/inv_call_{batch}.bed.gz')
    log:
        log='log/{asm_name}/inv_caller/log/{hap}/inv_call_{batch}.log'
    params:
        k_size=config.get('inv_k_size', 31),
        inv_threads=config.get('inv_threads', 4),
        inv_region_limit=config.get('inv_region_limit', None),
        inv_min_expand=config.get('inv_min_expand', None),
        srs_list = config.get('srs_list', None)
    shell:
        """
        python3 {PIPELINE_DIR}/scripts/call_inv_batch.py -f {input.bed_flag} -a {input.bed_aln} -t {input.tig_fa} -i {input.fai} -o {output.bed} -k {params.k_size} -d {params.inv_threads} --inv_region_limit {params.inv_region_limit} --inv_min_expand {params.inv_min_expand} --srs_list {params.srs_list} -b {wildcards.batch} -n {wildcards.asm_name} -p {wildcards.hap} -l {log.log} -R {REF_FA}
        """

#
# Flag regions to scan for inversions
#

# call_inv_merge_flagged_loci
#
# Merge flagged overlapping regions (SV ins/del match, indel ins/del match, indel cluster, SNV cluster) into single
# flagged record regions and annotate each. Column "TRY_INV" will be used by the inversion caller to determine which
# regions to try calling.
rule call_inv_merge_flagged_loci:
    input:
        bed_insdel_sv='temp/{asm_name}/inv_caller/flag/insdel_sv_{hap}.bed.gz',
        bed_insdel_indel='temp/{asm_name}/inv_caller/flag/insdel_indel_{hap}.bed.gz',
        bed_cluster_indel='temp/{asm_name}/inv_caller/flag/cluster_indel_{hap}.bed.gz',
        bed_cluster_snv='temp/{asm_name}/inv_caller/flag/cluster_snv_{hap}.bed.gz'
    output:
        bed='results/{asm_name}/inv_caller/flagged_regions_{hap}.bed.gz'
    params:
        flank=config.get('inv_sig_merge_flank', 500) , # Merge windows within this many bp
        batch_count=int(config.get('inv_sig_batch_count', BATCH_COUNT_DEFAULT)),  # Batch signature regions into this many batches for the caller. Marked here so that this file can be cross-referenced with the inversion caller log
        inv_sig_filter=config.get('inv_sig_filter', 'svindel')   # Filter flagged regions
    shell:
        '''
        python3 {PIPELINE_DIR}/scripts/call_inv_merge_flagged_loci.py -a {input.bed_insdel_sv} -b {input.bed_insdel_indel} -c {input.bed_cluster_indel} -d {input.bed_cluster_snv} -o {output.bed} -f {params.flank} -g {params.batch_count} -i {params.inv_sig_filter}
        '''


# call_inv_flag_insdel_cluster
#
# Flag inversion regions by matched INS/DEL calls, which often occur inside inversions. Aligners will generate a
# deletion over the inversion with an insertion of the same size flanking it where the inversion is the inverted
# sequence.
rule call_inv_flag_insdel_cluster:
    input:
        bed='temp/{asm_name}/cigar/pre_inv/svindel_insdel_{hap}.bed.gz'
    output:
        bed=temp('temp/{asm_name}/inv_caller/flag/insdel_{vartype}_{hap}.bed.gz')
    params:
        flank_cluster=int(config.get('inv_sig_insdel_cluster_flank', 2)), # For each INS, multiply SVLEN by this to get the distance to the nearest DEL it may intersect
        flank_merge=int(config.get('inv_sig_insdel_merge_flank', 2000)),  # Merge clusters within this distance
        cluster_min_svlen=int(config.get('inv_sig_cluster_svlen_min', 4))    # Discard indels less than this size
    wildcard_constraints:
        vartype='sv|indel'
    shell:
        '''
        python3 {PIPELINE_DIR}/scripts/call_inv_flag_insdel_cluster.py -i {input.bed} -o {output.bed} -v {wildcards.vartype} -c {params.flank_cluster} -d {params.flank_merge} -e {params.cluster_min_svlen}
        '''


# call_inv_cluster
#
# Detect clusters of indels and SNVs that common occur when contigs are aligned through inversions in direct orientation
rule call_inv_cluster:
    input:
        bed=_input_call_inv_cluster
    output:
        bed=temp('temp/{asm_name}/inv_caller/flag/cluster_{vartype}_{hap}.bed.gz')
    params:
        cluster_win=config.get('inv_sig_cluster_win', 200),            # Cluster variants within this many bases
        cluster_win_min=config.get('inv_sig_cluster_win_min', 500),    # Window must reach this size
        cluster_min_snv=config.get('inv_sig_cluster_snv_min', 20),     # Minimum number if SNVs in window (if vartype == snv)
        cluster_min_indel=config.get('inv_sig_cluster_indel_min', 10)  # Minimum number of indels in window (if vartype == indel)
    wildcard_constraints:
        vartype='indel|snv'
    shell:
        '''
        python3 {PIPELINE_DIR}/scripts/call_inv_cluster.py -i {input.bed} -o {output.bed} -v {wildcards.vartype} -w {params.cluster_win} -x {params.cluster_win_min} -y {params.cluster_min_snv} -z {params.cluster_min_indel}
        '''