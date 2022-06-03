version 1.0


task tar_asm {
  input {
    File ref
    File hapOne
    File hapTwo
    String sample
    String threads
    String mem_gb
  } command <<<
    mkdir -p asm/${sample}
    cp ${ref} asm/ref.fa
    samtools faidx asm/ref.fa
    cp ${hapOne} asm/${sample}/h1.fa.gz
    cp ${hapTwo} asm/${sample}/h2.fa.gz
    tar zcvf asm.tgz asm/
  >>>
  output {
    File asm_tar = "asm.tgz"
  }
}

task call_final_bed {
  input {
    File pav_conf
    File pav_sw
    File pav_asm
    File invBed
    File insBed
    File delBed
    File snvBed
    String threads
    String mem_gb
    String sample
  } command <<<
    cp ${pav_conf} ./
    tar zxvf ${pav_sw}
    tar zxvf ${pav_asm}
    tar zxvf ${invBed}
    tar zxvf ${snvBed}
    tar zxvf ${insBed}
    tar zxvf ${delBed}
    snakemake -s pav/Snakefile --cores ${threads} results/${sample}/bed/snv_snv.bed.gz results/${sample}/bed/indel_ins.bed.gz results/${sample}/bed/indel_del.bed.gz results/${sample}/bed/sv_ins.bed.gz results/${sample}/bed/sv_del.bed.gz results/${sample}/bed/sv_inv.bed.gz results/${sample}/bed/fa/indel_ins.fa.gz results/${sample}/bed/fa/indel_del.fa.gz results/${sample}/bed/fa/sv_ins.fa.gz results/${sample}/bed/fa/sv_del.fa.gz results/${sample}/bed/fa/sv_inv.fa.gz
  >>>
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
