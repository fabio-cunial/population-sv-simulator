"""
Rules for writing VCF output.
"""

# vcf_write_vcf
#
# Make VCF headers.
rule vcf_write_vcf:
    input:
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
        fa_sv_inv='results/{asm_name}/bed/fa/sv_inv.fa.gz',
        ref_tsv='data/ref/contig_info.tsv.gz',
        ref_fa=config['reference']
    output:
        vcf='pav_{asm_name}.vcf.gz'
    wildcard_constraints:
        alt_fmt='alt|sym'
    shell:
        '''
        python3 {PIPELINE_DIR}/scripts/vcf_write_vcf.py -a {input.bed_snv_snv} -b {input.bed_indel_ins} -c {input.bed_indel_del} -d {input.bed_sv_ins} -e {input.bed_sv_del} -f {input.bed_sv_inv} -g {input.fa_indel_ins} -i {input.fa_indel_del} -j {input.fa_sv_ins} -k {input.fa_sv_del} -l {input.fa_sv_inv} -m {input.ref_tsv} -o {output.vcf} -n {wildcards.asm_name} -r {input.ref_fa}
        '''