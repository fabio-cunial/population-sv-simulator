"""
Process alignments and alignment tiling paths.
"""

#
# Definitions
#

def _align_map_cpu(wildcards, config):

    if 'map_threads' in config:
        try:
            return int(config['map_threads'])
        except ValueError as ex:
            raise ValueError('Config parameter "map_threads" is not an integer: {map_threads}'.format(**config))

    return 12


#
# Alignment generation and processing
#

# align_cut_tig_overlap
#
# Cut contig alignment overlaps
rule align_cut_tig_overlap:
    input:
        bed='results/{asm_name}/align/pre-cut/aligned_tig_{hap}.bed.gz',
        tig_fai='temp/{asm_name}/align/contigs_{hap}.fa.gz.fai'
    output:
        bed='results/{asm_name}/align/aligned_tig_{hap}.bed.gz'
    params:
        min_trim_tig_len=np.int32(config.get('min_trim_tig_len', 1000))  # Minimum aligned tig length
    shell:
        '''        
        python3 {PIPELINE_DIR}/scripts/align_cut_tig_overlap.py -b {input.bed} -t {input.tig_fai} -m {params.min_trim_tig_len} -o {output.bed}
        '''


# align_get_read_bed
#
# Get alignment BED for one part (one aligned cell or split BAM) in one assembly.
rule align_get_read_bed:
    input:
        sam='temp/{asm_name}/align/pre-cut/aligned_tig_{hap}.sam.gz',
        tig_fai='temp/{asm_name}/align/contigs_{hap}.fa.gz.fai'
    output:
        bed='results/{asm_name}/align/pre-cut/aligned_tig_{hap}.bed.gz',
        align_head='results/{asm_name}/align/pre-cut/aligned_tig_{hap}.headers.gz'
    params:
        chrom_cluster=pavlib.util.as_bool(config.get('chrom_cluster', False))  # Assembly was clustered by chromosome and first part of chromosome name before "_" is the cluster name.
    wildcard_constraints:
        hap='h(0|1|2)'
    shell:
        """
        python3 {PIPELINE_DIR}/scripts/align_get_read_bed.py -s {input.sam} -t {input.tig_fai} -o {output.bed} -a {output.align_head} -c {params.chrom_cluster} -x {wildcards.hap}
        """



rule align_process_alignments:
    input:
        sam = expand("temp/{{asm_name}}/iter_{iter}/{{hap}}.sam", iter=['1', '2'])
    output:
        sam=temp('temp/{asm_name}/align/pre-cut/aligned_tig_{hap}.sam.gz')
    params:
        cpu=lambda wildcards: _align_map_cpu(wildcards, config)
    shell:
        '''
        cat <( samtools view -h {input.sam[0]} ) <( samtools view {input.sam[1]} ) | awk -vOFS="\\t" '($1 !~ /^@/) {{$10 = "*"; $11 = "*"}} {{print}}' | gzip > {output.sam}
        '''


rule align_map_inv:
    input:
        ref_fa='data/ref/ref.fa.gz',
        fa='temp/{asm_name}/align/contigs_{hap}.fa.gz' if config.get('aligner', 'minimap2') != 'lra' else 'temp/{asm_name}/align/contigs_{hap}.fa',
        aln="temp/{asm_name}/iter_1/{hap}.sam",
    output:
        aln=temp("temp/{asm_name}/iter_2/{hap}.sam"),
    params:
        cpu=lambda wildcards: _align_map_cpu(wildcards, config)
    shell:
        """
        minimap2 -K 8G -t {threads} \
            -ax asm20 -Y \
            --secondary=no --eqx -s 25000 \
            <(seqtk seq \
                -M <(samtools view -h {input.aln} | paftools.js sam2paf - | cut -f 6,8,9 | bedtools sort -i -) \
                -n "N" {input.ref_fa} \
            ) \
            <(seqtk seq \
                -M <(samtools view -h {input.aln} | paftools.js sam2paf - | cut -f 1,3,4 | bedtools sort -i -) \
                -n "N" {input.fa} \
            ) \
            | samtools view -F 4 -h - \
            > {output.aln}
        """

rule align_alignment:
    input:
        ref_fa='data/ref/ref.fa.gz',
        fa='temp/{asm_name}/align/contigs_{hap}.fa.gz' if config.get('aligner', 'minimap2') != 'lra' else 'temp/{asm_name}/align/contigs_{hap}.fa',
        gli='data/ref/ref.fa.gz.gli' if config.get('aligner', 'minimap2') == 'lra' else [],
        mmi='data/ref/ref.fa.gz.mmi' if config.get('aligner', 'minimap2') == 'lra' else []
    output:
        aln=temp("temp/{asm_name}/iter_1/{hap}.sam"),
    params:
        cpu=lambda wildcards: _align_map_cpu(wildcards, config)
    shell:
        """
        minimap2 -K 8G -t {params.cpu} \
            -ax asm20 \
            --secondary=no --eqx -Y -s 25000 \
            {input.ref_fa} {input.fa} \
            | samtools view -F 4 -h - \
            > {output.aln}
        """



# align_uncompress_tig
#
# Uncompress contig for aligners that cannot read gzipped FASTAs.
rule align_uncompress_tig:
    input:
        fa='temp/{asm_name}/align/contigs_{hap}.fa.gz'
    output:
        fa='temp/{asm_name}/align/contigs_{hap}.fa'
    shell:
        '''
        python3 {PIPELINE_DIR}/scripts/align_uncompress_tig.py -i {input.fa} -o {output.fa}
        '''



# align_get_tig_fa
#
# Get FASTA files.
rule align_get_tig_fa:
    input:
        fa=align_input_fasta
    output:
        fa=temp('temp/{asm_name}/align/contigs_{hap}.fa.gz'),
        fai=temp('temp/{asm_name}/align/contigs_{hap}.fa.gz.fai')
    shell:
        '''
        python3 {PIPELINE_DIR}/scripts/align_get_tig_fa.py -i {input.fa} -o {output.fa} -x {output.fai}
        '''

#
# Utilities
#
# Not needed for PAV, but useful rules for troubleshooting.

# align_get_cram_postcut
#
# Reconstruct CRAM from alignment BED files after trimming redundantly mapped bases (post-cut).
rule align_get_cram_postcut:
    input:
        bed='results/{asm_name}/align/aligned_tig_{hap}.bed.gz',
        fa=align_input_fasta,
        align_head='results/{asm_name}/align/pre-cut/aligned_tig_{hap}.headers.gz',
        ref_fa='data/ref/ref.fa.gz'

    output:
        cram='results/{asm_name}/align/aligned_tig_{hap}.cram'
    shell:
        """python3 {PIPELINE_DIR}/scripts/reconstruct_sam.py """
            """--bed {input.bed} --fasta {input.fa} --headers {input.align_head} | """
        """samtools view -T {input.ref_fa} -O CRAM -o {output.cram} && """
        """samtools index {output.cram}"""

# align_get_cram_precut
#
# Reconstruct CRAM from alignment BED files before trimming redundantly mapped bases (post-cut).
rule align_get_cram_precut:
    input:
        bed='results/{asm_name}/align/pre-cut/aligned_tig_{hap}.bed.gz',
        fa=align_input_fasta,
        align_head='results/{asm_name}/align/pre-cut/aligned_tig_{hap}.headers.gz',
        ref_fa='data/ref/ref.fa.gz'
    output:
        cram='results/{asm_name}/align/pre-cut/aligned_tig_{hap}.cram'
    shell:
        """python3 {PIPELINE_DIR}/scripts/reconstruct_sam.py """
            """--bed {input.bed} --fasta {input.fa} --headers {input.align_head} | """
        """samtools view -T {input.ref_fa} -O CRAM -o {output.cram} && """
        """samtools index {output.cram}"""
