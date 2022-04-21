task align_ref {
    File pav_conf
    File pav_sw
    File pav_asm
    File ref
    String threads
    String mem_gb
  command {
    cp ~{pav_conf} ./
    tar zxvf ~{pav_sw}
    tar zxvf ~{pav_asm}
    snakemake -s pav/Snakefile --cores ${threads} data/ref/ref.fa.gz data/ref/ref.fa.gz.fai
    tar zxvf align_ref.tgz data/ref/ref.fa.gz data/ref/ref.fa.gz.fai
  }
  output {
    File refGz = "align_ref.tgz"
  }
}

task align_get_tig_fa_h1 {
    File asm
    String sample
    String hap
    String threads
    String mem_gb
  command {
    cp ~{pav_conf} ./
    tar zxvf ~{pav_sw}
    tar zxvf ~{pav_asm}
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/align/contigs_${hap}.fa.gz temp/${sample}/align/contigs_${hap}.fa.gz.fai
    tar zxvf align_get_tig_fa_h1.tgz temp/${sample}/align/contigs_${hap}.fa.gz temp/${sample}/align/contigs_${hap}.fa.gz.fai
  }
  output {
    File asmGz = "align_get_tig_fa_h1.tgz"
  }
}

task align_get_tig_fa_h2 {
  File asm
  String hap
  String sample
  String threads
  String mem_gb
  command {
    cp ~{pav_conf} ./
    tar zxvf ~{pav_sw}
    tar zxvf ~{pav_asm}
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/align/contigs_${hap}.fa.gz temp/${sample}/align/contigs_${hap}.fa.gz.fai
    tar zxvf align_get_tig_fa_h2.tgz temp/${sample}/align/contigs_${hap}.fa.gz temp/${sample}/align/contigs_${hap}.fa.gz.fai
  }
  output {
    File asmGz = "align_get_tig_fa_h2.tgz"
  }
}

task align_ref_anno_n_gap {
  File refGz
  String sample
  File refFai
  String threads
  String mem_gb
  command {
    ln -s /net/eichler/vol27/projects/hprc/nobackups/variants/pav-wdl/pav
    snakemake -s pav/Snakefile --cores ${threads} data/ref/n_gap.bed.gz
  }
  output {
    File gaps = "data/ref/n_gap.bed.gz"
  }
}

task align_map_h1 {
  String hap
  String sample
  File refGz
  File asmGz
  String threads
  String mem_gb
  command {
    ln -s /net/eichler/vol27/projects/hprc/nobackups/variants/pav-wdl/pav
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/align/pre-cut/aligned_tig_${hap}.sam.gz
  }
  output {
    File  samGz = "temp/${sample}/align/pre-cut/aligned_tig_${hap}.sam.gz"
  }
}

task align_map_h2 {
    String hap
    String sample
    File refGz
    File asmGz
    String threads
    String mem_gb
  command {
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/align/pre-cut/aligned_tig_${hap}.sam.gz
  }
  output {
    File  samGz = "temp/${sample}/align/pre-cut/aligned_tig_${hap}.sam.gz"
  }
}

task align_get_read_bed_h1 {
    File refGz
    String sample
    File samGz
    String hap
    String threads
    String mem_gb
  command {
    snakemake -s pav/Snakefile --cores ${threads} results/${sample}/align/pre-cut/aligned_tig_${hap}.bed.gz
  }
  output {
    File bedGz = "results/${sample}/align/pre-cut/aligned_tig_${hap}.bed.gz"
    File headerGz = "results/${sample}/align/pre-cut/aligned_tig_${hap}.headers.gz"
  }
}

task align_get_read_bed_h2 {
    File refGz
    String sample
    File samGz
    String hap
    String threads
    String mem_gb
  command {
    snakemake -s pav/Snakefile --cores ${threads} results/${sample}/align/pre-cut/aligned_tig_${hap}.bed.gz
  }
  output {
    File bedGz = "results/${sample}/align/pre-cut/aligned_tig_${hap}.bed.gz"
    File headerGz = "results/${sample}/align/pre-cut/aligned_tig_${hap}.headers.gz"
  }
}

task align_cut_tig_overlap_h1 {
    String hap
    File asmFai
    File bedGz
    String threads
    String mem_gb
    String sample
  command {
    snakemake -s pav/Snakefile --cores ${threads} results/${sample}/align/aligned_tig_${hap}.bed.gz
  }
  output {
    File trimBed = "results/${sample}/align/aligned_tig_${hap}.bed.gz"
  }
}

task align_cut_tig_overlap_h2 {
    String hap
    File asmFai
    File bedGz
    String threads
    String mem_gb
    String sample
  command {
    snakemake -s pav/Snakefile --cores ${threads} results/${sample}/align/aligned_tig_${hap}.bed.gz
  }
  output {
    File trimBed = "results/${sample}/align/aligned_tig_${hap}.bed.gz"
  }
}
