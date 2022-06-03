version 1.0

task call_lg_split_hap {
  input {
    File pav_conf
    File pav_sw
    File pav_asm
    String sample
    File trimBed
    String hap
    String threads
    String mem_gb
  }
  command <<<
    cp ~{pav_conf} ./
    tar zxvf ~{pav_sw}
    tar zxvf ~{pav_asm}
    tar zxvf ~{trimBed}
    snakemake -s pav/Snakefile --cores ~{threads} temp/~{sample}/lg_sv/batch_~{hap}.tsv.gz
    tar zcvf call_lg_split_~{hap}_~{sample}.tgz temp/~{sample}/lg_sv/batch_~{hap}.tsv.gz
  >>>
  output {
    File batch = "call_lg_split_~{hap}_~{sample}.tgz"
  }
  ############################
  runtime {
      cpu:            threads
      memory:         mem_gb + " GiB"
      disks:          "local-disk " + 1000 + " HDD"
      bootDiskSizeGb: 50
      preemptible:    3
      maxRetries:     1
      docker:         "us.gcr.io/broad-dsp-lrma/lr-pav:1.2.1"
  }
}

task call_lg_discover_hap {
  input {
    String sample
    File trimBed
    File pav_conf
    File pav_sw
    File pav_asm
    String batch
    File batchFile
    File refGz
    File asmGz
    File gaps
    String hap
    String threads
    String mem_gb
  }
  command <<<
    cp ~{pav_conf} ./
    tar zxvf ~{pav_sw}
    tar zxvf ~{pav_asm}
    tar zxvf ~{refGz}
    tar zxvf ~{gaps}
    tar zxvf ~{asmGz}
    tar zxvf ~{trimBed}
    tar zxvf ~{batchFile}
    snakemake -s pav/Snakefile --cores ~{threads} temp/~{sample}/lg_sv/batch/sv_ins_~{hap}_~{batch}.bed.gz temp/~{sample}/lg_sv/batch/sv_del_~{hap}_~{batch}.bed.gz temp/~{sample}/lg_sv/batch/sv_inv_~{hap}_~{batch}.bed.gz
    tar zcvf call_lg_discover_~{sample}_~{hap}_~{batch}.tgz temp/~{sample}/lg_sv/batch/sv_ins_~{hap}_~{batch}.bed.gz temp/~{sample}/lg_sv/batch/sv_del_~{hap}_~{batch}.bed.gz temp/~{sample}/lg_sv/batch/sv_inv_~{hap}_~{batch}.bed.gz
  >>>
  output {
    File allBed = "call_lg_discover_~{sample}_~{hap}_~{batch}.tgz"
  }
  ############################
  runtime {
      cpu:            threads
      memory:         mem_gb + " GiB"
      disks:          "local-disk " + 1000 + " HDD"
      bootDiskSizeGb: 50
      preemptible:    3
      maxRetries:     1
      docker:         "us.gcr.io/broad-dsp-lrma/lr-pav:1.2.1"
  }
}

task call_merge_lg_del_h1 {
  input {
    File pav_conf
    File pav_sw
    File pav_asm
    String sample
    Array[File] inbed
    String hap
    String threads
    String mem_gb
    String svtype
  }
  command <<<
    cp ~{pav_conf} ./
    tar zxvf ~{pav_sw}
    tar zxvf ~{pav_asm}
    echo ~{sep=" " inbed} | tr " " "\n" | xargs -I '@' tar zxvf @
    snakemake -s pav/Snakefile --cores ~{threads} temp/~{sample}/lg_sv/sv_~{svtype}_~{hap}.bed.gz
    tar zcvf call_merge_lg_del_~{hap}_~{svtype}_~{sample}.tgz temp/~{sample}/lg_sv/sv_~{svtype}_~{hap}.bed.gz
  >>>
  output {
    File mergeBed = "call_merge_lg_del_~{hap}_~{svtype}_~{sample}.tgz"
  }
  ############################
  runtime {
      cpu:            threads
      memory:         mem_gb + " GiB"
      disks:          "local-disk " + 1000 + " HDD"
      bootDiskSizeGb: 50
      preemptible:    3
      maxRetries:     1
      docker:         "us.gcr.io/broad-dsp-lrma/lr-pav:1.2.1"
  }
}

task call_merge_lg_ins_h1 {
  input {
    File pav_conf
    File pav_sw
    File pav_asm
    String sample
    Array[File] inbed
    String hap
    String threads
    String mem_gb
    String svtype
  }
  command <<<
    cp ~{pav_conf} ./
    tar zxvf ~{pav_sw}
    tar zxvf ~{pav_asm}
    echo ~{sep=" " inbed} | tr " " "\n" | xargs -I '@' tar zxvf @
    snakemake -s pav/Snakefile --cores ~{threads} temp/~{sample}/lg_sv/sv_~{svtype}_~{hap}.bed.gz
    tar zcvf call_merge_lg_ins_~{hap}_~{svtype}_~{sample}.tgz temp/~{sample}/lg_sv/sv_~{svtype}_~{hap}.bed.gz
  >>>
  output {
    File mergeBed = "call_merge_lg_ins_~{hap}_~{svtype}_~{sample}.tgz"
  }
  ############################
  runtime {
      cpu:            threads
      memory:         mem_gb + " GiB"
      disks:          "local-disk " + 1000 + " HDD"
      bootDiskSizeGb: 50
      preemptible:    3
      maxRetries:     1
      docker:         "us.gcr.io/broad-dsp-lrma/lr-pav:1.2.1"
  }
}

task call_merge_lg_inv_h1 {
  input {
    File pav_conf
    File pav_sw
    File pav_asm
    String sample
    Array[File] inbed
    String hap
    String threads
    String mem_gb
    String svtype
  }
  command <<<
    cp ~{pav_conf} ./
    tar zxvf ~{pav_sw}
    tar zxvf ~{pav_asm}
    echo ~{sep=" " inbed} | tr " " "\n" | xargs -I '@' tar zxvf @
    snakemake -s pav/Snakefile --cores ~{threads} temp/~{sample}/lg_sv/sv_~{svtype}_~{hap}.bed.gz
    tar zcvf call_merge_lg_inv_~{hap}_~{svtype}_~{sample}.tgz temp/~{sample}/lg_sv/sv_~{svtype}_~{hap}.bed.gz
  >>>
  output {
    File mergeBed = "call_merge_lg_inv_~{hap}_~{svtype}_~{sample}.tgz"
  }
  ############################
  runtime {
      cpu:            threads
      memory:         mem_gb + " GiB"
      disks:          "local-disk " + 1000 + " HDD"
      bootDiskSizeGb: 50
      preemptible:    3
      maxRetries:     1
      docker:         "us.gcr.io/broad-dsp-lrma/lr-pav:1.2.1"
  }
}

task call_merge_lg_del_h2 {
  input {
    File pav_conf
    File pav_sw
    File pav_asm
    String sample
    Array[File] inbed
    String hap
    String threads
    String mem_gb
    String svtype
  }
  command <<<
    cp ~{pav_conf} ./
    tar zxvf ~{pav_sw}
    tar zxvf ~{pav_asm}
    echo ~{sep=" " inbed} | tr " " "\n" | xargs -I '@' tar zxvf @
    snakemake -s pav/Snakefile --cores ~{threads} temp/~{sample}/lg_sv/sv_~{svtype}_~{hap}.bed.gz
    tar zcvf call_merge_lg_del_~{hap}_~{svtype}_~{sample}.tgz temp/~{sample}/lg_sv/sv_~{svtype}_~{hap}.bed.gz
  >>>
  output {
    File mergeBed = "call_merge_lg_del_~{hap}_~{svtype}_~{sample}.tgz"
  }
  ############################
  runtime {
      cpu:            threads
      memory:         mem_gb + " GiB"
      disks:          "local-disk " + 1000 + " HDD"
      bootDiskSizeGb: 50
      preemptible:    3
      maxRetries:     1
      docker:         "us.gcr.io/broad-dsp-lrma/lr-pav:1.2.1"
  }
}

task call_merge_lg_ins_h2 {
  input {
    File pav_conf
    File pav_sw
    File pav_asm
    String sample
   Array[File] inbed
    String hap
    String threads
    String mem_gb
    String svtype
  }
  command <<<
    cp ~{pav_conf} ./
    tar zxvf ~{pav_sw}
    tar zxvf ~{pav_asm}
    echo ~{sep=" " inbed} | tr " " "\n" | xargs -I '@' tar zxvf @
    snakemake -s pav/Snakefile --cores ~{threads} temp/~{sample}/lg_sv/sv_~{svtype}_~{hap}.bed.gz
    tar zcvf call_merge_lg_ins_~{hap}_~{svtype}_~{sample}.tgz temp/~{sample}/lg_sv/sv_~{svtype}_~{hap}.bed.gz
  >>>
  output {
    File mergeBed = "call_merge_lg_ins_~{hap}_~{svtype}_~{sample}.tgz"
  }
  ############################
  runtime {
      cpu:            threads
      memory:         mem_gb + " GiB"
      disks:          "local-disk " + 1000 + " HDD"
      bootDiskSizeGb: 50
      preemptible:    3
      maxRetries:     1
      docker:         "us.gcr.io/broad-dsp-lrma/lr-pav:1.2.1"
  }
}

task call_merge_lg_inv_h2 {
  input {
    File pav_conf
    File pav_sw
    File pav_asm
    String sample
    Array[File] inbed
    String hap
    String threads
    String mem_gb
    String svtype
  }
  command <<<
    cp ~{pav_conf} ./
    tar zxvf ~{pav_sw}
    tar zxvf ~{pav_asm}
    echo ~{sep=" " inbed} | tr " " "\n" | xargs -I '@' tar zxvf @
    snakemake -s pav/Snakefile --cores ~{threads} temp/~{sample}/lg_sv/sv_~{svtype}_~{hap}.bed.gz
    tar zcvf call_merge_lg_inv_~{hap}_~{svtype}_~{sample}.tgz temp/~{sample}/lg_sv/sv_~{svtype}_~{hap}.bed.gz
  >>>
  output {
    File mergeBed = "call_merge_lg_inv_~{hap}_~{svtype}_~{sample}.tgz"
  }
  ############################
  runtime {
      cpu:            threads
      memory:         mem_gb + " GiB"
      disks:          "local-disk " + 1000 + " HDD"
      bootDiskSizeGb: 50
      preemptible:    3
      maxRetries:     1
      docker:         "us.gcr.io/broad-dsp-lrma/lr-pav:1.2.1"
  }
}
