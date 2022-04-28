task call_lg_split_h1 {
    File pav_conf
    File pav_sw
    File pav_asm
    String sample
    File trimBed
    String hap
    String threads
    String mem_gb
  command {
    cp ${pav_conf} ./
    tar zxvf ${pav_sw}
    tar zxvf ${pav_asm}
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/lg_sv/batch_${hap}.tsv.gz
  }
  output {
    File batch = "temp/${sample}/lg_sv/batch_${hap}.tsv.gz"
  }
}

task call_lg_split_h2 {
    String sample
    File pav_conf
    File pav_sw
    File pav_asm
    File trimBed
    String hap
    String threads
    String mem_gb
  command {
    cp ${pav_conf} ./
    tar zxvf ${pav_sw}
    tar zxvf ${pav_asm}
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/lg_sv/batch_${hap}.tsv.gz
  }
  output {
    File batch = "temp/${sample}/lg_sv/batch_${hap}.tsv.gz"
  }
}

task call_lg_discover_h1 {
    String sample
    File trimBed
    File pav_conf
    File pav_sw
    File pav_asm
    String batch
    File asmGz
    File gaps
    String hap
    String threads
    String mem_gb
  command {
    cp ${pav_conf} ./
    tar zxvf ${pav_sw}
    tar zxvf ${pav_asm}
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/lg_sv/batch/sv_ins_${hap}_${batch}.bed.gz temp/${sample}/lg_sv/batch/sv_del_${hap}_${batch}.bed.gz temp/${sample}/lg_sv/batch/sv_inv_${hap}_${batch}.bed.gz
  }
  output {
    File insBed = "temp/${sample}/lg_sv/batch/sv_ins_${hap}_${batch}.bed.gz"
    File delBed = "temp/${sample}/lg_sv/batch/sv_del_${hap}_${batch}.bed.gz"
    File invBed = "temp/${sample}/lg_sv/batch/sv_inv_${hap}_${batch}.bed.gz"
  }
}

task call_lg_discover_h2 {
    File pav_conf
    File pav_sw
    File pav_asm
    String sample
    File trimBed
    String batch
    File asmGz
    File gaps
    String hap
    String threads
    String mem_gb
  command {
    cp ${pav_conf} ./
    tar zxvf ${pav_sw}
    tar zxvf ${pav_asm}
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/lg_sv/batch/sv_ins_${hap}_${batch}.bed.gz temp/${sample}/lg_sv/batch/sv_del_${hap}_${batch}.bed.gz temp/${sample}/lg_sv/batch/sv_inv_${hap}_${batch}.bed.gz
  }
  output {
    File insBed = "temp/${sample}/lg_sv/batch/sv_ins_${hap}_${batch}.bed.gz"
    File delBed = "temp/${sample}/lg_sv/batch/sv_del_${hap}_${batch}.bed.gz"
    File invBed = "temp/${sample}/lg_sv/batch/sv_inv_${hap}_${batch}.bed.gz"
  }
}


task call_merge_lg_del_h1 {
    File pav_conf
    File pav_sw
    File pav_asm
    String sample
    File inbed
    String hap
    String threads
    String mem_gb
    String svtype
  command {
    cp ${pav_conf} ./
    tar zxvf ${pav_sw}
    tar zxvf ${pav_asm}
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/lg_sv/sv_${svtype}_${hap}.bed.gz
  }
  output {
    File mergeBed = "temp/${sample}/lg_sv/sv_${svtype}_${hap}.bed.gz"
  }
}

task call_merge_lg_ins_h1 {
    File pav_conf
    File pav_sw
    File pav_asm
    String sample
    File inbed
    String hap
    String threads
    String mem_gb
    String svtype
  command {
    cp ${pav_conf} ./
    tar zxvf ${pav_sw}
    tar zxvf ${pav_asm}
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/lg_sv/sv_${svtype}_${hap}.bed.gz
  }
  output {
    File mergeBed = "temp/${sample}/lg_sv/sv_${svtype}_${hap}.bed.gz"
  }
}

task call_merge_lg_inv_h1 {
    File pav_conf
    File pav_sw
    File pav_asm
    String sample
    File inbed
    String hap
    String threads
    String mem_gb
    String svtype
  command {
    cp ${pav_conf} ./
    tar zxvf ${pav_sw}
    tar zxvf ${pav_asm}
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/lg_sv/sv_${svtype}_${hap}.bed.gz
  }
  output {
    File mergeBed = "temp/${sample}/lg_sv/sv_${svtype}_${hap}.bed.gz"
  }
}

task call_merge_lg_del_h2 {
    File pav_conf
    File pav_sw
    File pav_asm
    String sample
    File inbed
    String hap
    String threads
    String mem_gb
    String svtype
  command {
    cp ${pav_conf} ./
    tar zxvf ${pav_sw}
    tar zxvf ${pav_asm}
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/lg_sv/sv_${svtype}_${hap}.bed.gz
  }
  output {
    File mergeBed = "temp/${sample}/lg_sv/sv_${svtype}_${hap}.bed.gz"
  }
}

task call_merge_lg_ins_h2 {
    File pav_conf
    File pav_sw
    File pav_asm
    String sample
    File inbed
    String hap
    String threads
    String mem_gb
    String svtype
  command {
    cp ${pav_conf} ./
    tar zxvf ${pav_sw}
    tar zxvf ${pav_asm}
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/lg_sv/sv_${svtype}_${hap}.bed.gz
  }
  output {
    File mergeBed = "temp/${sample}/lg_sv/sv_${svtype}_${hap}.bed.gz"
  }
}

task call_merge_lg_inv_h2 {
    File pav_conf
    File pav_sw
    File pav_asm
    String sample
    File inbed
    String hap
    String threads
    String mem_gb
    String svtype
  command {
    cp ${pav_conf} ./
    tar zxvf ${pav_sw}
    tar zxvf ${pav_asm}
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/lg_sv/sv_${svtype}_${hap}.bed.gz
  }
  output {
    File mergeBed = "temp/${sample}/lg_sv/sv_${svtype}_${hap}.bed.gz"
  }
}
