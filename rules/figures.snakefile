"""
Generate figures from variant discovery.

The rules in this snakefile may not generate variant calls. The pipeline should be run first, then variants can be
plotted with rules in this file. These may be useful for reporting results or troubleshooting.
"""

#
# Inversions
#

# figures_inv_dot_density
#
# Make dot and density plot for an inversion call. May plot with inverted repeat whiskers (whiskers=whisk) or without
# (whiskers=nowhisk). These whiskers show k-mers unique to each inverted repeat flanking the inversion and whether they
# are in reference orientation (pointing up) or inverted orientation (pointing down).
rule figures_inv_dot_density:
    input:
        tig_fa='temp/{asm_name}/align/contigs_{hap}.fa.gz',
        bed_inv='results/{asm_name}/{flag_source}/sv_inv_{hap}.bed.gz',
        bed_density='results/{asm_name}/{flag_source}/density_table/density_{inv_id}_{hap}.tsv.gz'
    output:
        fig_dot='figures/inv/{asm_name}/{hap}/{flag_source}/{inv_id}_dot.{ext}',
        fig_den='figures/inv/{asm_name}/{hap}/{flag_source}/{inv_id}_density.{ext}'
    wildcard_constraints:
        ext='pdf|png'
    shell:
        '''
        python3 {PIPELINE_DIR}/scripts/figures_inv_dot_density.py -f {input.tig_fa} -i {input.bed_inv} -d {input.bed_density} -p {output.fig_dot} -q {output.fig_den} -d {wildcards.inv_id} -a {wildcards.asm_name} -x {wildcards.hap}
        '''

