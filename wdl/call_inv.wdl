task call_inv_flag_insdel_cluster_indel_h1 {
    String sample
    File inbed
    String hap
    String threads
    String mem_gb
    String vartype
  command {
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz
  }
  output {
    File bed = "temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz"
  }
}

task call_inv_flag_insdel_cluster_sv_h1 {
    String sample
    File inbed
    String hap
    String threads
    String mem_gb
    String vartype
  command {
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz
  }
  output {
    File bed = "temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz"
  }
}

task call_inv_flag_insdel_cluster_indel_h2 {
    String sample
    File inbed
    String hap
    String threads
    String mem_gb
    String vartype
  command {
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz
  }
  output {
    File bed = "temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz"
  }
}

task call_inv_flag_insdel_cluster_sv_h2 {
    String sample
    File inbed
    String hap
    String threads
    String mem_gb
    String vartype
  command {
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz
  }
  output {
    File bed = "temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz"
  }
}

task call_inv_cluster_sv_h1 {
    String sample
    File inbed
    String threads
    String mem_gb
    String hap
    String vartype
  command {
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz
  }
  output {
    File bed = "temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz"
  }
}

task call_inv_cluster_indel_h1 {
    String sample
    File inbed
    String threads
    String mem_gb
    String hap
    String vartype
  command {
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz
  }
  output {
    File bed = "temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz"
  }
}

task call_inv_cluster_snv_h1 {
    String sample
    File inbed
    String threads
    String mem_gb
    String hap
    String vartype
  command {
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz
  }
  output {
    File bed = "temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz"
  }
}

task call_inv_cluster_sv_h2 {
    String sample
    File inbed
    String threads
    String mem_gb
    String hap
    String vartype
  command {
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz
  }
  output {
    File bed = "temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz"
  }
}

task call_inv_cluster_indel_h2 {
    String sample
    File inbed
    String threads
    String mem_gb
    String hap
    String vartype
  command {
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz
  }
  output {
    File bed = "temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz"
  }
}

task call_inv_cluster_snv_h2 {
    String sample
    File inbed
    String threads
    String mem_gb
    String hap
    String vartype
  command {
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz
  }
  output {
    File bed = "temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz"
  }
}



task call_inv_merge_flagged_loci_h1 {
    String sample
    String hap
    File indelFlag
    File svFlag
    File snvCluster
    File svCluster
    File indelCluster
    String threads
    String mem_gb
   command {
  snakemake -s pav/Snakefile --cores results/${sample}/inv_caller/flagged_regions_${hap}.bed.gz
   }
   output {
  File bed = "results/${sample}/inv_caller/flagged_regions_${hap}.bed.gz"
   }
}

task call_inv_merge_flagged_loci_h2 {
    String sample
    String hap
    File indelFlag
    File svFlag
    File snvCluster
    File svCluster
    File indelCluster
    String threads
    String mem_gb
   command {
  snakemake -s pav/Snakefile --cores results/${sample}/inv_caller/flagged_regions_${hap}.bed.gz
   }
   output {
  File bed = "results/${sample}/inv_caller/flagged_regions_${hap}.bed.gz"
   }
}

task call_inv_batch_h1 {
    String sample
    String hap
    File trimBed
    File flag
    File asmFai
    File asmGz
    String batch
    String threads
    String mem_gb
  command {
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/inv_caller/batch/${hap}/inv_call_${batch}.bed.gz
  }
  output {
    File bed = "temp/${sample}/inv_caller/batch/${hap}/inv_call_${batch}.bed.gz"
  }
}

task call_inv_batch_h2 {
    String sample
    String hap
    File trimBed
    File flag
    File asmFai
    File asmGz
    String batch
    String threads
    String mem_gb
  command {
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/inv_caller/batch/${hap}/inv_call_${batch}.bed.gz
  }
  output {
    File bed = "temp/${sample}/inv_caller/batch/${hap}/inv_call_${batch}.bed.gz"
  }
}

task call_inv_batch_merge_h1 {
    String sample
    String hap
    File invBed
    String threads
    String mem_gb
  command {
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/inv_caller/sv_inv_${hap}.bed.gz
  }
  output {
    File bed = "temp/${sample}/inv_caller/sv_inv_${hap}.bed.gz"
  }
}

task call_inv_batch_merge_h2 {
    String sample
    String hap
    File invBed
    String threads
    String mem_gb
  command {
    snakemake -s pav/Snakefile --cores ${threads} temp/${sample}/inv_caller/sv_inv_${hap}.bed.gz
  }
  output {
    File bed = "temp/${sample}/inv_caller/sv_inv_${hap}.bed.gz"
  }
}