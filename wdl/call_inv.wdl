task call_inv_flag_insdel_cluster_indel_h1 {
    String sample
    File pav_conf
    File pav_sw
    File pav_asm
    File inbed
    String hap
    String threads
    String mem_gb
    String vartype
  command {
    cp ${pav_conf} ./
    tar zxvf ${pav_sw}
    tar zxvf ${pav_asm}
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz
    tar zcvf call_inv_flag_insdel_cluster_indel_${sample}_${vartype}_${hap}.tgz temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz
  }
  output {
    File bed = "call_inv_flag_insdel_cluster_indel_${sample}_${vartype}_${hap}.tgz"
  }
}

task call_inv_flag_insdel_cluster_sv_h1 {
    String sample
    File pav_conf
    File pav_sw
    File pav_asm
    File inbed
    String hap
    String threads
    String mem_gb
    String vartype
  command {
    cp ${pav_conf} ./
    tar zxvf ${pav_sw}
    tar zxvf ${pav_asm}
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz
    tar zcvf call_inv_flag_insdel_cluster_indel_${sample}_${vartype}_${hap}.tgz temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz
  }
  output {
    File bed = "call_inv_flag_insdel_cluster_indel_${sample}_${vartype}_${hap}.tgz"
  }
}

task call_inv_flag_insdel_cluster_indel_h2 {
    String sample
    File inbed
    File pav_conf
    File pav_sw
    File pav_asm
    String hap
    String threads
    String mem_gb
    String vartype
  command {
    cp ${pav_conf} ./
    tar zxvf ${pav_sw}
    tar zxvf ${pav_asm}
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz
    tar zcvf call_inv_flag_insdel_cluster_indel_${sample}_${vartype}_${hap}tgz temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz
  }
  output {
    File bed = "call_inv_flag_insdel_cluster_indel_${sample}_${vartype}_${hap}.tgz"
  }
}

task call_inv_flag_insdel_cluster_sv_h2 {
    String sample
    File inbed
    String hap
    File pav_conf
    File pav_sw
    File pav_asm
    String threads
    String mem_gb
    String vartype
  command {
    cp ${pav_conf} ./
    tar zxvf ${pav_sw}
    tar zxvf ${pav_asm}
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz
    tar zcvf call_inv_flag_insdel_cluster_indel_${sample}_${vartype}_${hap}.tgz temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz
  }
  output {
    File bed = "call_inv_flag_insdel_cluster_indel_${sample}_${vartype}_${hap}.tgz"
  }
}

task call_inv_cluster_sv_h1 {
    String sample
    File pav_conf
    File pav_sw
    File pav_asm
    File inbed
    String threads
    String mem_gb
    String hap
    String vartype
  command {
    cp ${pav_conf} ./
    tar zxvf ${pav_sw}
    tar zxvf ${pav_asm}
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/inv_caller/flag/sv_${vartype}_${hap}.bed.gz
    tar zcvf call_inv_cluster_sv_${sample}_${vartype}_${hap}.tgz temp/${sample}/inv_caller/flag/sv_${vartype}_${hap}.bed.gz
  }
  output {
    File bed = "call_inv_cluster_sv_${sample}_${vartype}_${hap}.tgz"
  }
}

task call_inv_cluster_indel_h1 {
    String sample
    File pav_conf
    File pav_sw
    File pav_asm
    File inbed
    String threads
    String mem_gb
    String hap
    String vartype
  command {
    cp ${pav_conf} ./
    tar zxvf ${pav_sw}
    tar zxvf ${pav_asm}
    tar zxvf ${inbed}
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz
    tar zcvf call_inv_cluster_indel_${hap}_${sample}_${vartype}.tgz temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz
  }
  output {
    File bed = "call_inv_cluster_indel_${hap}_${sample}_${vartype}.tgz"
  }
}

task call_inv_cluster_snv_h1 {
    String sample
    File pav_conf
    File pav_sw
    File pav_asm
    File inbed
    String threads
    String mem_gb
    String hap
    String vartype
  command {
    cp ${pav_conf} ./
    tar zxvf ${pav_sw}
    tar zxvf ${pav_asm}
    tar zxvf ${inbed}
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz
    tar zcvf call_inv_cluster_snv_${hap}_${sample}_${vartype}.tgz temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz
  }
  output {
    File bed = "call_inv_cluster_snv_${hap}_${sample}_${vartype}.tgz"
  }
}

task call_inv_cluster_sv_h2 {
    String sample
    File pav_conf
    File pav_sw
    File pav_asm
    File inbed
    String threads
    String mem_gb
    String hap
    String vartype
  command {
    cp ${pav_conf} ./
    tar zxvf ${pav_sw}
    tar zxvf ${pav_asm}
    tar zxvf ${inbed}
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz
    tar zcvf call_inv_cluster_sv_${sample}_${vartype}_${hap}.tgz temp/${sample}/inv_caller/flag/sv_${vartype}_${hap}.bed.gz
  }
  output {
    File bed = "call_inv_cluster_sv_${sample}_${vartype}_${hap}.tgz"
  }
}

task call_inv_cluster_indel_h2 {
    String sample
    File pav_conf
    File pav_sw
    File pav_asm
    File inbed
    String threads
    String mem_gb
    String hap
    String vartype
  command {
    cp ${pav_conf} ./
    tar zxvf ${pav_sw}
    tar zxvf ${pav_asm}
    tar zxvf ${inbed}
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz
    tar zcvf call_inv_cluster_indel_${hap}_${sample}_${vartype}.tgz temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz
  }
  output {
    File bed = "call_inv_cluster_indel_${hap}_${sample}_${vartype}.tgz"
  }
}

task call_inv_cluster_snv_h2 {
    String sample
    File inbed
    File pav_conf
    File pav_sw
    File pav_asm
    String threads
    String mem_gb
    String hap
    String vartype
  command {
    cp ${pav_conf} ./
    tar zxvf ${pav_sw}
    tar zxvf ${pav_asm}
    tar zxvf ${inbed}
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz
    tar zcvf call_inv_cluster_snv_${hap}_${sample}_${vartype}.tgz temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz
  }
  output {
    File bed = "call_inv_cluster_snv_${hap}_${sample}_${vartype}.tgz"
  }
}



task call_inv_merge_flagged_loci_h1 {
    String sample
    String hap
    File pav_conf
    File pav_sw
    File pav_asm
    File indelFlag
    File svFlag
    File snvCluster
    File svCluster
    File indelCluster
    String threads
    String mem_gb
  command {
    cp ${pav_conf} ./
    tar zxvf ${pav_sw}
    tar zxvf ${pav_asm}
    tar zxvf ${indelFlag}
    tar zxvf ${svFlag}
    tar zxvf ${snvCluster}
    tar zxvf ${svCluster}
    tar zxvf ${indelCluster}
    snakemake -s pav/Snakefile --cores results/${sample}/inv_caller/flagged_regions_${hap}.bed.gz
    tar zcvf call_inv_merge_flagged_loci_${hap}_${sample}.tgz results/${sample}/inv_caller/flagged_regions_${hap}.bed.gz
  }
  output {
    File bed = "call_inv_merge_flagged_loci_${hap}_${sample}.tgz"
  }
}

task call_inv_merge_flagged_loci_h2 {
    String sample
    String hap
    File indelFlag
    File pav_conf
    File pav_sw
    File pav_asm
    File svFlag
    File snvCluster
    File svCluster
    File indelCluster
    String threads
    String mem_gb
   command {
    cp ${pav_conf} ./
    tar zxvf ${pav_sw}
    tar zxvf ${pav_asm}
    tar zxvf ${indelFlag}
    tar zxvf ${svFlag}
    tar zxvf ${snvCluster}
    tar zxvf ${svCluster}
    tar zxvf ${indelCluster}
    snakemake -s pav/Snakefile --cores results/${sample}/inv_caller/flagged_regions_${hap}.bed.gz
    tar zcvf call_inv_merge_flagged_loci_${hap}_${sample}.tgz results/${sample}/inv_caller/flagged_regions_${hap}.bed.gz
  }
  output {
    File bed = "call_inv_merge_flagged_loci_${hap}_${sample}.tgz"
  }
}

task call_inv_batch_h1 {
    String sample
    File pav_conf
    File pav_sw
    File pav_asm
    String hap
    File trimBed
    File flag
    File asmGz
    String batch
    String threads
    String mem_gb
  command {
    cp ${pav_conf} ./
    tar zxvf ${pav_sw}
    tar zxvf ${pav_asm}
    tar zxvf ${trimBed}
    tar zxvf ${flag}
    tar zxvf ${asmGz}
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/inv_caller/batch/${hap}/inv_call_${batch}.bed.gz
    tar zcvf call_inv_batch_${hap}_${batch}_${sample}.tgz temp/${sample}/inv_caller/batch/${hap}/inv_call_${batch}.bed.gz
  }
  output {
    File bed = "call_inv_batch_${hap}_${batch}_${sample}.tgz"
  }
}

task call_inv_batch_h2 {
    String sample
    String hap
    File trimBed
    File flag
    File pav_conf
    File pav_sw
    File pav_asm
    File asmGz
    String batch
    String threads
    String mem_gb
  command {
    cp ${pav_conf} ./
    tar zxvf ${pav_sw}
    tar zxvf ${pav_asm}
    tar zxvf ${trimBed}
    tar zxvf ${flag}
    tar zxvf ${asmGz}
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/inv_caller/batch/${hap}/inv_call_${batch}.bed.gz
    tar zcvf call_inv_batch_${hap}_${batch}_${sample}.tgz temp/${sample}/inv_caller/batch/${hap}/inv_call_${batch}.bed.gz
  }
  output {
    File bed = "call_inv_batch_${hap}_${batch}_${sample}.tgz"
  }
}

task call_inv_batch_merge_h1 {
    String sample
    File pav_conf
    File pav_sw
    File pav_asm
    String hap
    Array[File] invBed
    String threads
    String mem_gb
  command {
    cp ${pav_conf} ./
    tar zxvf ${pav_sw}
    tar zxvf ${pav_asm}
    echo ${invBed} | xargs -I '@' tar zxvf @
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/inv_caller/sv_inv_${hap}.bed.gz
    tar zcvf call_inv_batch_merge_${hap}_${sample}.tgz temp/${sample}/inv_caller/sv_inv_${hap}.bed.gz
  }
  output {
    File bed = "call_inv_batch_merge_${hap}_${sample}.tgz "
  }
}


task call_inv_batch_merge_h2 {
    String sample
    String hap
    File pav_conf
    File pav_sw
    File pav_asm
    Array[File] invBed
    String threads
    String mem_gb
  command {
    cp ${pav_conf} ./
    tar zxvf ${pav_sw}
    tar zxvf ${pav_asm}
    echo ${invBed} | xargs -I '@' tar zxvf @
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/inv_caller/sv_inv_${hap}.bed.gz
    tar zcvf call_inv_batch_merge_${hap}_${sample}.tgz temp/${sample}/inv_caller/sv_inv_${hap}.bed.gz
  }
  output {
    File bed = "call_inv_batch_merge_${hap}_${sample}.tgz "
  }
}