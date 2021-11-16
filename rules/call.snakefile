"""
Call variants from aligned contigs.
"""


#
# Definitions
#


#
# Finalize variant calls
#

# call_final_bed
#
# Separate variants into final BED and FA files and write to results.
rule call_final_bed:
    input:
        bed_ins='temp/{asm_name}/bed/merged/svindel_ins.bed.gz',
        bed_del='temp/{asm_name}/bed/merged/svindel_del.bed.gz',
        bed_inv='temp/{asm_name}/bed/merged/sv_inv.bed.gz',
        bed_snv='temp/{asm_name}/bed/merged/snv_snv.bed.gz'
    output:
        bed_snv_snv='results/{asm_name}/bed/snv_snv.bed.gz',
        bed_indel_ins='results/{asm_name}/bed/indel_ins.bed.gz',
        bed_indel_del='results/{asm_name}/bed/indel_del.bed.gz',
        bed_sv_ins='results/{asm_name}/bed/sv_ins.bed.gz',
        bed_sv_del='results/{asm_name}/bed/sv_del.bed.gz',
        bed_sv_inv='results/{asm_name}/bed/sv_inv.bed.gz',
        fa_indel_ins='results/{asm_name}/bed/fa/indel_ins.fa.gz',
        fa_indel_del='results/{asm_name}/bed/fa/indel_del.fa.gz',
        fa_sv_ins='results/{asm_name}/bed/fa/sv_ins.fa.gz',
        fa_sv_del='results/{asm_name}/bed/fa/sv_del.fa.gz',
        fa_sv_inv='results/{asm_name}/bed/fa/sv_inv.fa.gz'
    shell:
        '''
        python3 {PIPELINE_DIR}/scripts/call_final_bed.py -i {input.bed_ins} -d {input.bed_del} -n {input.bed_inv} -s {input.bed_snv} -o {output.bed_snv_snv} -p {output.bed_indel_ins} -q {output.bed_indel_del} -r {output.bed_sv_ins} -z {output.bed_sv_del} -t {output.bed_sv_inv} -u {output.fa_indel_ins} -v {output.fa_indel_del} -w {output.fa_sv_ins} -x {output.fa_sv_del} -y {output.fa_sv_inv}
        '''

#
# Merge haplotypes
#

# pg_variant_bed
#
# Make variant BED and FASTA. Write all variant calls regardless of consensus loci.
#
# If merging is not done by chromosome (default), then this rule reads from each haplotype and merges in one step
# (wildcards.merge_level = merged). If merging is done per chromosome (config['merge_by_chrom'] is True), then this
# rule calls itself recursively first by chromosome (wildcards.merge_level = bychrom) then by concatenating the merged
# chromosomes (wildcards.merge_level = merged). The code will know which step its on based on the wildcards and config.
rule call_merge_haplotypes:
    input:
        bed_var_h1='temp/{asm_name}/bed/integrated/h1/{vartype_svtype}.bed.gz',
        bed_var_h2='temp/{asm_name}/bed/integrated/h2/{vartype_svtype}.bed.gz',
        callable_h1='results/{asm_name}/callable/callable_regions_h1_500.bed.gz',
        callable_h2='results/{asm_name}/callable/callable_regions_h2_500.bed.gz',
        bed_chrom=lambda wildcards: [
            'temp/{asm_name}/bed/bychrom/{vartype_svtype}/{chrom}.bed.gz'.format(
                asm_name=wildcards.asm_name, vartype_svtype=wildcards.vartype_svtype, chrom=chrom
            ) for chrom in sorted(svpoplib.ref.get_df_fai(config['reference'] + '.fai').index)
        ] if pavlib.util.as_bool(config.get('merge_by_chrom', True)) else []
    output:
        bed=temp('temp/{asm_name}/bed/merged/{vartype_svtype}.bed.gz')
    params:
        ro_min=float(config.get('ro_min', 0.5)),
        offset_max=int(config.get('offset_max', 200)),
        merge_threads=int(config.get('merge_threads', 12))
    shell:
        '''
        python3 {PIPELINE_DIR}/scripts/call_merge_haplotypes.py -v {input.bed_var_h1} -w {input.bed_var_h2} -c {input.callable_h1} -d {input.callable_h2} -o {output.bed} -r {params.ro_min} -f {params.offset_max} -m {params.merge_threads} -s {wildcards.vartype_svtype} --chrom {input.bed_chrom}
        '''

# call_merge_haplotypes_chrom
#
# Merge by chromosome. This rule is used if "merge_by_chrom" is True.
rule call_merge_haplotypes_chrom:
    input:
        bed_var_h1='temp/{asm_name}/bed/integrated/h1/{vartype_svtype}.bed.gz',
        bed_var_h2='temp/{asm_name}/bed/integrated/h2/{vartype_svtype}.bed.gz',
        callable_h1='results/{asm_name}/callable/callable_regions_h1_500.bed.gz',
        callable_h2='results/{asm_name}/callable/callable_regions_h2_500.bed.gz',
    output:
        bed='temp/{asm_name}/bed/bychrom/{vartype_svtype}/{chrom}.bed.gz'
    params:
        ro_min=float(config.get('ro_min', 0.5)),
        offset_max=int(config.get('offset_max', 200)),
        merge_threads=int(config.get('merge_threads', 12))
    shell:
        '''
        python3 {PIPELINE_DIR}/scripts/call_merge_haplotypes_chrom.py -v {input.bed_var_h1} -w {input.bed_var_h2} -c {input.callable_h1} -d {input.callable_h2} -o {output.bed} -r {params.ro_min} -f {params.offset_max} -m {params.merge_threads} -s {wildcards.vartype_svtype} -x {wildcards.chrom}
        '''


#
# Integrate variant calls
#

# call_mappable_bed
#
# Make a table of mappable regions by merging aligned loci with loci covered by alignment-truncating events.
# "flank" parameter is an integer describing how far away records may be to merge (similar to the "bedtools merge"
# "slop" parameter). The flank is not added to the regions that are output.
rule call_mappable_bed:
    input:
        bed_align='results/{asm_name}/align/aligned_tig_{hap}.bed.gz',
        bed_lg_del='temp/{asm_name}/lg_sv/sv_del_{hap}.bed.gz',
        bed_lg_ins='temp/{asm_name}/lg_sv/sv_ins_{hap}.bed.gz',
        bed_lg_inv='temp/{asm_name}/lg_sv/sv_inv_{hap}.bed.gz'
    output:
        bed='results/{asm_name}/callable/callable_regions_{hap}_{flank}.bed.gz'
    shell:
        '''
        python3 {PIPELINE_DIR}/scripts/call_mappable_bed.py -a {input.bed_align} -d {input.bed_lg_del} -i {input.bed_lg_ins} -n {input.bed_lg_inv} -f {wildcards.flank} -o {output.bed}
        '''

# call_correct_inter_inv
#
# Filter variants from inside inversions
rule call_integrate_sources:
    input:
        bed_cigar_insdel='temp/{asm_name}/cigar/pre_inv/svindel_insdel_{hap}.bed.gz',
        bed_cigar_snv='temp/{asm_name}/cigar/pre_inv/snv_snv_{hap}.bed.gz',
        bed_lg_ins='temp/{asm_name}/lg_sv/sv_ins_{hap}.bed.gz',
        bed_lg_del='temp/{asm_name}/lg_sv/sv_del_{hap}.bed.gz',
        bed_lg_inv='temp/{asm_name}/lg_sv/sv_inv_{hap}.bed.gz',
        bed_inv='temp/{asm_name}/inv_caller/sv_inv_{hap}.bed.gz'
    output:
        bed_ins=temp('temp/{asm_name}/bed/integrated/{hap}/svindel_ins.bed.gz'),
        bed_del=temp('temp/{asm_name}/bed/integrated/{hap}/svindel_del.bed.gz'),
        bed_inv=temp('temp/{asm_name}/bed/integrated/{hap}/sv_inv.bed.gz'),
        bed_snv=temp('temp/{asm_name}/bed/integrated/{hap}/snv_snv.bed.gz')
    params:
        min_inv=config.get('min_inv', 300),
        max_inv=config.get('max_inv', 2000000),
        tig_filt=config.get('tig_filter_pattern')
    shell:
        '''
        python3 {PIPELINE_DIR}/scripts/call_integrate_sources.py -c {input.bed_cigar_insdel} -b {input.bed_cigar_snv} -d {input.bed_lg_ins} -e {input.bed_lg_del} -f {input.bed_lg_inv} -g {input.bed_inv} -i {output.bed_ins} -j {output.bed_del} -k {output.bed_inv} -l {output.bed_snv} -m {params.min_inv} -n {params.max_inv} -a {wildcards.asm_name} -x {wildcards.hap} -t {params.tig_filt}
        '''


# call_inv_bed
#
# Make inversion call BED.
rule call_inv_bed:
    input:
        bed='temp/{asm_name}/inv_caller/sv_inv_{hap}.bed.gz'
    output:
        bed=temp('temp/{asm_name}/bed/pre_merge/sv_inv_{hap}.bed.gz'),
        bed_dropped='results/{asm_name}/bed/dropped/shortinv_sv_inv_{hap}.bed.gz'
    params:
        min_svlen=config.get('inv_min_svlen', 300)
    shell:
        '''
        python3 {PIPELINE_DIR}/scripts/call_inv_bed.py -i {input.bed} -o {output.bed} -d {output.bed_dropped} -m {params.min_svlen}
        '''

#
# Call from CIGAR
#

# call_cigar_merge
#
# Merge discovery sets from each batch.
rule call_cigar_merge:
    input:
        bed_insdel=expand('temp/{{asm_name}}/cigar/batched/insdel_{{hap}}_{batch}.bed.gz', batch=range(pavlib.cigarcall.CALL_CIGAR_BATCH_COUNT)),
        bed_snv=expand('temp/{{asm_name}}/cigar/batched/snv.bed_{{hap}}_{batch}.gz', batch=range(pavlib.cigarcall.CALL_CIGAR_BATCH_COUNT))
    output:
        bed_insdel=temp('temp/{asm_name}/cigar/pre_inv/svindel_insdel_{hap}.bed.gz'),
        bed_snv=temp('temp/{asm_name}/cigar/pre_inv/snv_snv_{hap}.bed.gz')
    shell:
        """
        python3 {PIPELINE_DIR}/scripts/call_cigar_merge.py -i {input.bed_insdel} -n {input.bed_snv} -o {output.bed_insdel} -s {output.bed_snv}
        """


# call_cigar
#
# Call variants by alignment CIGAR parsing.
rule call_cigar:
    input:
        bed='results/{asm_name}/align/aligned_tig_{hap}.bed.gz',
        tig_fa_name='temp/{asm_name}/align/contigs_{hap}.fa.gz'
    output:
        bed_insdel=temp('temp/{asm_name}/cigar/batched/insdel_{hap}_{batch}.bed.gz'),
        bed_snv=temp('temp/{asm_name}/cigar/batched/snv.bed_{hap}_{batch}.gz')
    shell:
        """
        python3 {PIPELINE_DIR}/scripts/call_cigar.py -i {input.bed} -f {input.tig_fa_name} -b {wildcards.batch} -o {output.bed_insdel} -s {output.bed_snv} -r {REF_FA} -x {wildcards.hap}
        """

