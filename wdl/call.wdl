task call_cigar_h1 {
    File trimBed
    File asmGz
    String batch
    String threads
    String sample
    String mem_gb
    String hap
  command {
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/cigar/batched/insdel_${hap}_${batch}.bed.gz temp/${sample}/cigar/batched/snv.bed_${hap}_${batch}.gz

  }
  output {
    File snvBed = "temp/${sample}/cigar/batched/snv.bed_${hap}_${batch}.gz"
    File insdelBed = "temp/${sample}/cigar/batched/insdel_${hap}_${batch}.bed.gz"
  }
}

task call_cigar_h2 {
    String sample
    File trimBed
    File asmGz
    String batch
    String threads
    String mem_gb
    String hap
  command {
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/cigar/batched/insdel_${hap}_${batch}.bed.gz temp/${sample}/cigar/batched/snv.bed_${hap}_${batch}.gz
  }
  output {
    File snvBed = "temp/${sample}/cigar/batched/snv.bed_${hap}_${batch}.gz"
    File insdelBed = "temp/${sample}/cigar/batched/insdel_${hap}_${batch}.bed.gz"
  }
}


task call_cigar_merge_h1 {
    String sample
    String hap
    File insdelBatch
    File snvBatch
    String threads
    String mem_gb
  command {
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/cigar/pre_inv/svindel_insdel_${hap}.bed.gz temp/${sample}/cigar/pre_inv/snv_snv_${hap}.bed.gz
  }
  output {
    File insdelBedMerge = "temp/${sample}/cigar/pre_inv/svindel_insdel_${hap}.bed.gz"
    File snvBedMerge = "temp/${sample}/cigar/pre_inv/snv_snv_${hap}.bed.gz"
  }
}

task call_cigar_merge_h2 {
    String sample
    String hap
    File insdelBatch
    File snvBatch
    String threads
    String mem_gb
  command {
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/cigar/pre_inv/svindel_insdel_${hap}.bed.gz temp/${sample}/cigar/pre_inv/snv_snv_${hap}.bed.gz
  }
  output {
    File insdelBedMerge = "temp/${sample}/cigar/pre_inv/svindel_insdel_${hap}.bed.gz"
    File snvBedMerge = "temp/${sample}/cigar/pre_inv/snv_snv_${hap}.bed.gz"
  }
}


task call_mappable_bed_h1 {
    String sample
    File delBed
    File insBed
    File invBed
    File trimBed
    String hap
    String threads
    String mem_gb
  command {
    snakemake -s pav/Snakefile --cores ${threads} results/${sample}/callable/callable_regions_${hap}_500.bed.gz
  }
  output {
    File bed = "results/${sample}/callable/callable_regions_${hap}_500.bed.gz"
  }
}

task call_mappable_bed_h2 {
    String sample
    File delBed
    File insBed
    File invBed
    File trimBed
    String hap
    String threads
    String mem_gb
  command {
    snakemake -s pav/Snakefile --cores ${threads} results/${sample}/callable/callable_regions_${hap}_500.bed.gz
  }
  output {
    File bed = "results/${sample}/callable/callable_regions_${hap}_500.bed.gz"
  }
}

task call_integrate_sources_h1 {
    String sample
    String hap
    File preInvSvBed
    File preInvSnvBed
    File delBedIn
    File insBedIn
    File invBedIn
    File invBatch
    String threads
    String mem_gb
  command {
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/bed/integrated/${hap}/svindel_ins.bed.gz temp/${sample}/bed/integrated/${hap}/svindel_del.bed.gz temp/${sample}/bed/integrated/${hap}/snv_snv.bed.gz temp/${sample}/bed/integrated/${hap}/svindel_inv.bed.gz
  }
  output {
    File insBed = "temp/${sample}/bed/integrated/${hap}/svindel_ins.bed.gz"
    File delBed = "temp/${sample}/bed/integrated/${hap}/svindel_del.bed.gz"
    File snvBed = "temp/${sample}/bed/integrated/${hap}/snv_snv.bed.gz"
    File invBed = "temp/${sample}/bed/integrated/${hap}/svindel_inv.bed.gz"
  }
}

task call_integrate_sources_h2 {
    String sample
    String hap
    File preInvSvBed
    File preInvSnvBed
    File delBedIn
    File insBedIn
    File invBedIn
    File invBatch
    String threads
    String mem_gb
  command {
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/bed/integrated/${hap}/svindel_ins.bed.gz temp/${sample}/bed/integrated/${hap}/svindel_del.bed.gz temp/${sample}/bed/integrated/${hap}/snv_snv.bed.gz temp/${sample}/bed/integrated/${hap}/svindel_inv.bed.gz
  }
  output {
    File insBed = "temp/${sample}/bed/integrated/${hap}/svindel_ins.bed.gz"
    File delBed = "temp/${sample}/bed/integrated/${hap}/svindel_del.bed.gz"
    File snvBed = "temp/${sample}/bed/integrated/${hap}/snv_snv.bed.gz"
    File invBed = "temp/${sample}/bed/integrated/${hap}/svindel_inv.bed.gz"
  }
}

task call_merge_haplotypes_chrom_svindel_del {
    String svtype
    File delBed_h1
    String sample
    File delBed_h2
    File callable_h1
    File callable_h2
    String chrom
    String threads
    String mem_gb
  command {
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/bed/bychrom/${svtype}/${chrom}.bed.gz
  }
  output {
    File bed = "temp/${sample}/bed/bychrom/${svtype}/${chrom}.bed.gz"
  }
}

task call_merge_haplotypes_chrom_svindel_ins {
    String svtype
    File insBed_h1
    File insBed_h2
    File callable_h1
    String sample
    File callable_h2
    String chrom
    String threads
    String mem_gb
  command {
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/bed/bychrom/${svtype}/${chrom}.bed.gz
  }
  output {
    File bed = "temp/${sample}/bed/bychrom/${svtype}/${chrom}.bed.gz"
  }
}

task call_merge_haplotypes_chrom_svinv {
    String svtype
    File invBed_h1
    File invBed_h2
    File callable_h1
    File callable_h2
    String chrom
    String sample
    String threads
    String mem_gb
  command {
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/bed/bychrom/${svtype}/${chrom}.bed.gz
  }
  output {
    File bed = "temp/${sample}/bed/bychrom/${svtype}/${chrom}.bed.gz"
  }
}

task call_merge_haplotypes_chrom_snv {
    String svtype
    File snvBed_h1
    File snvBed_h2
    File callable_h1
    String sample
    File callable_h2
    String chrom
    String threads
    String mem_gb
  command {
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/bed/bychrom/${svtype}/${chrom}.bed.gz
  }
  output {
    File bed = "temp/${sample}/bed/bychrom/${svtype}/${chrom}.bed.gz"
  }
}

task call_merge_haplotypes_inv {
    String svtype
    File inbed
    String sample
    String threads
    String mem_gb
  command {
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/bed/merged/${svtype}.bed.gz
  }
  output {
    File bed = "temp/${sample}/bed/merged/${svtype}.bed.gz"
  }
}

task call_merge_haplotypes_svindel_del {
    String svtype
    File inbed
    String sample
    String threads
    String mem_gb
  command {
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/bed/merged/${svtype}.bed.gz
  }
  output {
    File bed = "temp/${sample}/bed/merged/${svtype}.bed.gz"
  }
}

task call_merge_haplotypes_svindel_ins {
    String svtype
    File inbed
    String sample
    String threads
    String mem_gb
  command {
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/bed/merged/${svtype}.bed.gz
  }
  output {
    File bed = "temp/${sample}/bed/merged/${svtype}.bed.gz"
  }
}

task call_merge_haplotypes_snv {
    String svtype
    String sample
    File inbed
    String threads
    String mem_gb
  command {
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/bed/merged/${svtype}.bed.gz
  }
  output {
    File bed = "temp/${sample}/bed/merged/${svtype}.bed.gz"
  }
}

task call_final_bed {
    File invBed
    File insBed
    File delBed
    File snvBed
    String threads
    String mem_gb
    String sample
  command {
    snakemake -s pav/Snakefile --cores ${threads} results/${sample}/bed/snv_snv.bed.gz results/${sample}/bed/indel_ins.bed.gz results/${sample}/bed/indel_del.bed.gz results/${sample}/bed/sv_ins.bed.gz results/${sample}/bed/sv_del.bed.gz results/${sample}/bed/sv_inv.bed.gz results/${sample}/bed/fa/indel_ins.fa.gz results/${sample}/bed/fa/indel_del.fa.gz results/${sample}/bed/fa/sv_ins.fa.gz results/${sample}/bed/fa/sv_del.fa.gz results/${sample}/bed/fa/sv_inv.fa.gz
  }
  output {
    File snvBedOut = "results/${sample}/bed/snv_snv.bed.gz"
    File indelInsBed = "results/${sample}/bed/indel_ins.bed.gz"
    File indelDelBed = "results/${sample}/bed/indel_del.bed.gz"
    File svInsBed = "results/${sample}/bed/sv_ins.bed.gz"
    File svDelBed = "results/${sample}/bed/sv_del.bed.gz"
    File invBedOut = "results/${sample}/bed/sv_inv.bed.gz"
    File indelInsFasta = "results/${sample}/bed/fa/indel_ins.fa.gz"
    File indelDelFasta = "results/${sample}/bed/fa/indel_del.fa.gz"
    File svInsFasta = "results/${sample}/bed/fa/sv_ins.fa.gz"
    File svDelFasta = "results/${sample}/bed/fa/sv_del.fa.gz"
    File invFasta = "results/${sample}/bed/fa/sv_inv.fa.gz"
  }
}
