task align_ref {
    File pav_conf
    File pav_sw
    File pav_asm
    String sample
    String threads
    String mem_gb
  command {
    cp ${pav_conf} ./
    tar zxvf ${pav_sw}
    tar zxvf ${pav_asm}
    snakemake -s pav/Snakefile --cores ${threads} data/ref/ref.fa.gz data/ref/ref.fa.gz.fai
    tar zcvf align_ref_${sample}.tgz data/ref/ref.fa.gz data/ref/ref.fa.gz.fai
  }
  output {
    File refGz = "align_ref_${sample}.tgz"
  }
}

task align_get_tig_fa_h1 {
    File pav_conf
    File pav_sw
    File pav_asm
    String sample
    String hap
    String threads
    String mem_gb
  command {
    cp ${pav_conf} ./
    tar zxvf ${pav_sw}
    tar zxvf ${pav_asm}
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/align/contigs_${hap}.fa.gz temp/${sample}/align/contigs_${hap}.fa.gz.fai
    tar zcvf align_get_tig_fa_${hap}_${sample}.tgz temp/${sample}/align/contigs_${hap}.fa.gz temp/${sample}/align/contigs_${hap}.fa.gz.fai
  }
  output {
    File asmGz = "align_get_tig_fa_${hap}_${sample}.tgz"
  }
}

task align_get_tig_fa_h2 {
    File pav_conf
    File pav_sw
    File pav_asm
    String hap
    String sample
    String threads
    String mem_gb
  command {
    cp ${pav_conf} ./
    tar zxvf ${pav_sw}
    tar zxvf ${pav_asm}
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/align/contigs_${hap}.fa.gz temp/${sample}/align/contigs_${hap}.fa.gz.fai
    tar zcvf align_get_tig_fa_${hap}_${sample}.tgz temp/${sample}/align/contigs_${hap}.fa.gz temp/${sample}/align/contigs_${hap}.fa.gz.fai
  }
  output {
    File asmGz = "align_get_tig_fa_${hap}_${sample}.tgz"
  }
}

task align_ref_anno_n_gap {
    File pav_conf
    File pav_sw
    File pav_asm
    String sample
    String threads
    String mem_gb
  command {
    cp ${pav_conf} ./
    tar zxvf ${pav_sw}
    tar zxvf ${pav_asm}
    snakemake -s pav/Snakefile --cores ${threads} data/ref/n_gap.bed.gz
    tar zcvf align_ref_anno_n_gap_${sample}.tgz data/ref/n_gap.bed.gz
  }
  output {
    File gaps = "align_ref_anno_n_gap_${sample}.tgz"
  }
}

task align_map_h1 {
  File pav_conf
  File pav_sw
  File pav_asm
  String hap
  String sample
  File asmGz
  File refGz
  String threads
  String mem_gb
  command {
    cp ${pav_conf} ./
    tar zxvf ${pav_sw}
    tar zxvf ${pav_asm}
    tar zxvf ${asmGz}
    tar zxvf ${refGz}
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/align/pre-cut/aligned_tig_${hap}.sam.gz
    tar czvf align_map_${hap}_${sample}.tgz temp/${sample}/align/pre-cut/aligned_tig_${hap}.sam.gz
  }
  output {
    File  samGz = "align_map_${hap}_${sample}.tgz"
  }
}

task align_map_h2 {
  File pav_conf
  File pav_sw
  File pav_asm
  String hap
  String sample
  File refGz
  File asmGz
  String threads
  String mem_gb
  command {
    cp ${pav_conf} ./
    tar zxvf ${pav_sw}
    tar zxvf ${pav_asm}
    tar zxvf ${refGz}
    tar zxvf ${asmGz}
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/align/pre-cut/aligned_tig_${hap}.sam.gz
    tar czvf align_map_${hap}_${sample}.tgz temp/${sample}/align/pre-cut/aligned_tig_${hap}.sam.gz
  }
  output {
    File samGz = "align_map_${hap}_${sample}.tgz"
  }
}

task align_get_read_bed_h1 {
    File pav_conf
    File pav_sw
    File pav_asm
    File refGz
    String sample
    File samGz
    String hap
    String threads
    String mem_gb
  command {
    cp ${pav_conf} ./
    tar zxvf ${pav_sw}
    tar zxvf ${pav_asm}
    snakemake -s pav/Snakefile --cores ${threads} results/${sample}/align/pre-cut/aligned_tig_${hap}.bed.gz
    tar czvf align_get_read_bed_${hap}_${sample}.tgz results/${sample}/align/pre-cut/aligned_tig_${hap}.bed.gz
  }
  output {
    File bedGz = "align_get_read_bed_${hap}_${sample}.tgz"
  }
}

task align_get_read_bed_h2 {
    File pav_conf
    File pav_sw
    File pav_asm
    File refGz
    String sample
    File samGz
    String hap
    String threads
    String mem_gb
  command {
    cp ${pav_conf} ./
    tar zxvf ${pav_sw}
    tar zxvf ${pav_asm}
    snakemake -s pav/Snakefile --cores ${threads} results/${sample}/align/pre-cut/aligned_tig_${hap}.bed.gz
    tar czvf align_get_read_bed_${hap}_${sample}.tgz results/${sample}/align/pre-cut/aligned_tig_${hap}.bed.gz
  }
  output {
    File bedGz = "align_get_read_bed_${hap}_${sample}.tgz"
  }
}

task align_cut_tig_overlap_h1 {
    File pav_conf
    File pav_sw
    File pav_asm
    File asmGz
    String hap
    File bedGz
    String threads
    String mem_gb
    String sample
  command {
    cp ${pav_conf} ./
    tar zxvf ${pav_sw}
    tar zxvf ${pav_asm}
    tar zxvf ${asmGz}
    tar zxvf ${bedGz}
    snakemake -s pav/Snakefile --cores ${threads} results/${sample}/align/aligned_tig_${hap}.bed.gz
    tar czvf align_cut_tig_overlap_${hap}_${sample}.tgz results/${sample}/align/aligned_tig_${hap}.bed.gz
  }
  output {
    File trimBed = "align_cut_tig_overlap_${hap}_${sample}.tgz"
  }
}

task align_cut_tig_overlap_h2 {
    File pav_conf
    File pav_sw
    File pav_asm
    String hap
    File asmGz
    File bedGz
    String threads
    String mem_gb
    String sample
  command {
    cp ${pav_conf} ./
    tar zxvf ${pav_sw}
    tar zxvf ${pav_asm}
    tar zxvf ${asmGz}
    tar zxvf ${bedGz}
    snakemake -s pav/Snakefile --cores ${threads} results/${sample}/align/aligned_tig_${hap}.bed.gz
    tar czvf align_cut_tig_overlap_${hap}_${sample}.tgz results/${sample}/align/aligned_tig_${hap}.bed.gz
  }
  output {
    File trimBed = "align_cut_tig_overlap_${hap}_${sample}.tgz"
  }
}
