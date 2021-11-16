"""
Call alignment-truncating events (large SVs).
"""

# call_merge_lg
#
# Merge variant calls from large SVs.
rule call_merge_lg:
    input:
        bed=expand(
            'temp/{{asm_name}}/lg_sv/batch/sv_{{svtype}}_{{hap}}_{batch}.bed.gz',
            batch=range(config.get('lg_batch_count', 10))
        )
    output:
        bed=temp('temp/{asm_name}/lg_sv/sv_{svtype}_{hap}.bed.gz')
    shell:
        '''
        python3 {PIPELINE_DIR}/scripts/call_merge_lg.py -i {input.bed} -o {output.bed}
        '''

# call_lg_discover
#
# Call alignment-truncating SVs.
rule call_lg_discover:
    input:
        bed='results/{asm_name}/align/aligned_tig_{hap}.bed.gz',
        tsv_group='temp/{asm_name}/lg_sv/batch_{hap}.tsv.gz',
        fa='temp/{asm_name}/align/contigs_{hap}.fa.gz',
        fai='temp/{asm_name}/align/contigs_{hap}.fa.gz.fai',
        bed_n='data/ref/n_gap.bed.gz'
    output:
        bed_ins=temp('temp/{asm_name}/lg_sv/batch/sv_ins_{hap}_{batch}.bed.gz'),
        bed_del=temp('temp/{asm_name}/lg_sv/batch/sv_del_{hap}_{batch}.bed.gz'),
        bed_inv=temp('temp/{asm_name}/lg_sv/batch/sv_inv_{hap}_{batch}.bed.gz')
    log:
        log='log/{asm_name}/lg_sv/log/lg_sv_{hap}_{batch}.log'
    params:
        k_size=int(config.get('inv_k_size', 31)),
        inv_threads_lg=int(config.get('inv_threads_lg', 12)),
        inv_region_limit=config.get('inv_region_limit', None),
        srs_list = config.get('srs_list', None)
    shell:
        '''
        python3 {PIPELINE_DIR}/scripts/call_lg_discover.py -e {input.bed} -t {input.tsv_group} -f {input.fa} -i {input.fai} -n {input.bed_n} -o {output.bed_ins} -d {output.bed_del} -v {output.bed_inv} -k {params.k_size} -l {params.inv_threads_lg} --inv_region_limit {params.inv_region_limit} --srs_list {params.srs_list}  -x {wildcards.hap} -b {wildcards.batch} -a {wildcards.asm_name} -z {log.log} -r {REF_FA}
        '''


# call_lg_split
#
# Split chromosome/tig records into batches for chromosome/tig pairs with multiple alignments.
rule call_lg_split:
    input:
        bed='results/{asm_name}/align/aligned_tig_{hap}.bed.gz'
    output:
        tsv=temp('temp/{asm_name}/lg_sv/batch_{hap}.tsv.gz')
    params:
        batch_count=config.get('lg_batch_count', 10)
    shell:
        '''
        python3 {PIPELINE_DIR}/scripts/call_lg_split.py -i {input.bed} -o {output.tsv} -b {params.batch_count}
        '''
