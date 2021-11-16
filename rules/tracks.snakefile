"""
Prepare UCSC tracks for data.
"""

#
# Definitions
#

VARTYPE_TO_SVTYPE_TUPLE = {
    'snv': ('snv'),
    'sv': ('ins', 'del'),
    'indel': ('ins', 'del')
}

ALIGN_COLOR = {
    'h1': '160,00,144',
    'h2': '64,00,160'
}

def _track_get_input_bed(wildcards):
    """
    Get one or more input files for tracks. If "svtype" is "all", collect all relevant input files.

    :param wildcards: Wildcards.

    :return: List of input file(s).
    """

    if wildcards.svtype == 'insdel':

        if wildcards.vartype not in {'sv', 'indel'}:
            raise RuntimeError('Cannot collect "all" for variant type: {} (must be a valid variant type with more than one svtype)'.format(wildcards.vartype))

        return [
            'results/{{asm_name}}/bed/pre_merge/{{vartype}}_{svtype}_{{hap}}.bed.gz'.format(svtype=svtype)
            for svtype in VARTYPE_TO_SVTYPE_TUPLE[wildcards.vartype]
        ]

    return ['results/{asm_name}/bed/pre_merge/{vartype}_{svtype}_{hap}.bed.gz']

def _get_align_bed(wildcards):

    if wildcards.align_stage == 'pre-cut':
        return [
            'results/{asm_name}/align/pre-cut/aligned_tig_h1.bed.gz'.format(asm_name=wildcards.asm_name),
            'results/{asm_name}/align/pre-cut/aligned_tig_h2.bed.gz'.format(asm_name=wildcards.asm_name)
        ]

    if wildcards.align_stage == 'post-cut':
        return [
            'results/{asm_name}/align/aligned_tig_h1.bed.gz'.format(asm_name=wildcards.asm_name),
            'results/{asm_name}/align/aligned_tig_h2.bed.gz'.format(asm_name=wildcards.asm_name)
        ]

    raise RuntimeError('Unknown align_stage wildcard: '.format(wildcards.align_stage))

#
# Variant calls
#

# tracks_hap_call_bb
#
# BigBed for one variant set.
rule tracks_hap_call_bb:
    input:
        bed='temp/{asm_name}/tracks/bed_pre_merge/{vartype}_{svtype}_{hap}.bed',
        asfile='temp/{asm_name}/tracks/bed_pre_merge/{vartype}_{svtype}_{hap}.as',
        fai=REF_FAI
    output:
        bb='tracks/{asm_name}/variant/pre_merge/{vartype}_{svtype}_{hap}.bb'
    shell:
        """bedToBigBed -tab -as={input.asfile} -type=bed9+ {input.bed} {input.fai} {output.bb}"""

# tracks_hap_call
#
# Tracks for one variant set.
rule tracks_hap_call:
    input:
        bed=_track_get_input_bed,
        fai=REF_FAI
    output:
        bed=temp('temp/{asm_name}/tracks/bed_pre_merge/{vartype}_{svtype}_{hap}.bed'),
        asfile=temp('temp/{asm_name}/tracks/bed_pre_merge/{vartype}_{svtype}_{hap}.as')
    shell:
        '''
        python3 {PIPELINE_DIR}/scripts/tracks_align.py -i {input.bed} -f {input.fai} -o {output.bed} -a {output.asfile} -n {wildcards.asm_name} -v {wildcards.vartype} -s {wildcards.svtype} -x {wildcards.hap}
        '''


#
# Alignments
#

# tracks_align_bb
#
# Alignment track BED to BigBed.
rule tracks_align_bb:
    input:
        bed='temp/{asm_name}/tracks/align/{align_stage}.bed',
        asfile='temp/{asm_name}/tracks/align/{align_stage}.as',
        fai=REF_FAI
    output:
        bb='tracks/{asm_name}/align/{align_stage}.bb'
    shell:
        """bedToBigBed -tab -as={input.asfile} -type=bed9+ {input.bed} {input.fai} {output.bb}"""

# tracks_align
#
# Alignment tracks.
rule tracks_align:
    input:
        bed=_get_align_bed
    output:
        bed=temp('temp/{asm_name}/tracks/align/{align_stage}.bed'),
        asfile=temp('temp/{asm_name}/tracks/align/{align_stage}.as')
    wildcard_constraints:
        align_stage='(pre|post)-cut'
    shell:
        '''
        python3 {PIPELINE_DIR}/scripts/tracks_align.py -i {input.bed} -o {output.bed} -a {output.asfile} -l {wildcards.align_stage}
        '''