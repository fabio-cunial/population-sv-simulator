import "wdl/setup.wdl" as setup
import "wdl/align.wdl" as align
import "wdl/call_pav.wdl" as call_pav
import "wdl/call_inv.wdl" as call_inv
import "wdl/call_lg.wdl" as call_lg


workflow pav {
  File ref
  File hapOne
  File hapTwo
  String sample
  File refFai
  File config
  File pav_tar
  Array[Array[String]] chroms = read_tsv(refFai)

  call setup.tar_asm {
    input:
      ref = ref,
      hapOne = hapOne,
      hapTwo = hapTwo,
      threads = "1",
      mem_gb = "2"
  }
  call align.align_ref {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      threads = "1",
      mem_gb = "8",
  }
  call align.align_get_tig_fa_h1 {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      asm = hapOne,
      hap = "h1",
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  call align.align_get_tig_fa_h2 {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      asm = hapTwo,
      hap = "h2",
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  call align.align_ref_anno_n_gap {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      refGz = align.align_ref.refGz,
      refFai = align.align_ref.refFai,
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  call align.align_map_h1 {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      refGz = align.align_ref.refGz,
      hap = "h1",
      asmGz = align.align_get_tig_fa_h2.asmGz,
      threads = "8",
      mem_gb = "12",
      sample = sample
  }
  call align.align_map_h2 {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      refGz = align.align_ref.refGz,
      hap = "h2",
      asmGz = align.align_get_tig_fa_h2.asmGz,
      threads = "8",
      mem_gb = "12",
      sample = sample
  }
  call align.align.align_get_read_bed_h1 {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      refGz = align.align_ref.refGz,
      hap = "h1",
      samGz = align.align_map_h1.samGz,
      threads = "1",
      mem_gb = "32",
      sample = sample
  }
  call align.align_get_read_bed_h2 {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      refGz = align.align_ref.refGz,
      hap = "h2",
      samGz = align.align_map_h2.samGz,
      threads = "1",
      mem_gb = "32",
      sample = sample
  }
  call align.align_cut_tig_overlap_h1 {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      hap = "h1",
      bedGz = align.align_get_read_bed_h1.bedGz,
      asmFai = align.align_get_tig_fa_h1.asmFai,
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  call align.align_cut_tig_overlap_h2 {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      hap = "h2",
      bedGz = align.align_get_read_bed_h2.bedGz,
      asmFai = align.align_get_tig_fa_h2.asmFai,
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  scatter(i in range(10)) {
     call call_pav.call_cigar_h1 {
      input:
        pav_conf = config,
        pav_sw = pav_tar,
        pav_asm = setup.tar_asm.asm_tar,
        hap = "h1",
        trimBed = align.align_cut_tig_overlap_h1.trimBed,
        asmGz = align.align_get_tig_fa_h1.asmGz,
        batch = i,
        threads = "1",
        mem_gb = "24",
        sample = sample
     }
  }
  scatter(i in range(10)) {
     call call_pav.call_cigar_h2 {
      input:
        pav_conf = config,
        pav_sw = pav_tar,
        pav_asm = setup.tar_asm.asm_tar,
        hap = "h2",
        trimBed = align.align_cut_tig_overlap_h2.trimBed,
        asmGz = align.align_get_tig_fa_h2.asmGz,
        batch = i,
        threads = "1",
        mem_gb = "24",
        sample = sample
     }
  }
  call call_pav.call_lg_split_h1 {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      hap = "h1",
      trimBed = align.align_cut_tig_overlap_h2.trimBed,
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  call call_pav.call_lg_split_h2 {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      hap = "h2",
      trimBed = align.align_cut_tig_overlap_h2.trimBed,
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  scatter(i in range(10)) {
     call call_pav.call_lg_discover_h1 {
      input:
        pav_conf = config,
        pav_sw = pav_tar,
        pav_asm = setup.tar_asm.asm_tar,
        hap = "h1",
        trimBed = align.align_cut_tig_overlap_h1.trimBed,
        gaps = align.align_ref_anno_n_gap.gaps,
        asmFai = align.align_get_tig_fa_h1.asmFai,
        asmGz = align.align_get_tig_fa_h1.asmGz,
        batch = call_pav.call_lg_split_h1.batch,
        threads = "8",
        mem_gb = "8",
        sample = sample
     }
  }
  scatter(i in range(10)) {
     call call_pav.call_lg_discover_h2 {
      input:
        pav_conf = config,
        pav_sw = pav_tar,
        pav_asm = setup.tar_asm.asm_tar,
        hap = "h2",
        trimBed = align.align_cut_tig_overlap_h2.trimBed,
        gaps = align.align_ref_anno_n_gap.gaps,
        asmFai = align.align_get_tig_fa_h2.asmFai,
        asmGz = align.align_get_tig_fa_h2.asmGz,
        batch = call_pav.call_lg_split_h2.batch,
        threads = "8",
        mem_gb = "8",
        sample = sample
     }
  }
  call call_pav.call_cigar_merge_h1 {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      hap = "h1",
      insdelBatch = call_pav.call_cigar_h1.insdelBed,
      snvBatch = call_pav.call_cigar_h1.snvBed,
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  call call_pav.call_cigar_merge_h2 {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      hap = "h2",
      insdelBatch = call_pav.call_cigar_h2.insdelBed,
      snvBatch = call_pav.call_cigar_h2.snvBed,
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  call call_pav.call_merge_lg_del_h1 {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      hap = "h1",
      svtype = "del",
      inbed = call_lg.call_lg_discover_h1.delBed,
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  call call_pav.call_merge_lg_ins_h1 {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      hap = "h1",
      svtype = "ins",
      inbed = call_lg.call_lg_discover_h1.insBed,
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  call call_pav.call_merge_lg_inv_h1 {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      hap = "h1",
      svtype = "inv",
      inbed = call_lg.call_lg_discover_h1.invBed,
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  call call_pav.call_merge_lg_del_h2 {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      hap = "h2",
      svtype = "del",
      inbed = call_lg.call_lg_discover_h2.delBed,
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  call call_pav.call_merge_lg_ins_h2 {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      hap = "h2",
      svtype = "ins",
      inbed = call_lg.call_lg_discover_h2.insBed,
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  call call_pav.call_merge_lg_inv_h2 {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      hap = "h2",
      svtype = "inv",
      inbed = call_lg.call_lg_discover_h2.invBed,
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  call call_inv.call_inv_cluster_indel_h1 {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      hap = "h1",
      vartype="indel",
      inbed = call_pav.call_cigar_merge_h1.insdelBedMerge,
      threads = "4",
      mem_gb = "16",
      sample = sample
  }
  call call_inv.call_inv_cluster_snv_h1 {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      hap = "h1",
      vartype="snv",
      inbed = call_pav.call_cigar_merge_h1.insdelBedMerge,
      threads = "4",
      mem_gb = "16",
      sample = sample
  }
  call call_inv.call_inv_cluster_sv_h1 {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      hap = "h1",
      vartype="sv",
      inbed = call_pav.call_cigar_merge_h1.insdelBedMerge,
      threads = "4",
      mem_gb = "16",
      sample = sample
  }
  call call_inv.call_inv_cluster_indel_h2 {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      hap = "h2",
      vartype="indel",
      inbed = call_pav.call_cigar_merge_h2.insdelBedMerge,
      threads = "4",
      mem_gb = "16",
      sample = sample
  }
  call call_inv.call_inv_cluster_snv_h2 {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      hap = "h2",
      vartype="snv",
      inbed = call_pav.call_cigar_merge_h2.insdelBedMerge,
      threads = "4",
      mem_gb = "16",
      sample = sample
  }
  call call_inv.call_inv_cluster_sv_h2 {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      hap = "h2",
      vartype="sv",
      inbed = call_pav.call_cigar_merge_h2.insdelBedMerge,
      threads = "4",
      mem_gb = "16",
      sample = sample
  }
  call call_inv.call_inv_flag_insdel_cluster_indel_h1 {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      hap = "h1",
      vartype="indel",
      inbed = call_pav.call_cigar_merge_h1.insdelBedMerge,
      threads = "4",
      mem_gb = "16",
      sample = sample
  }
  call call_inv.call_inv_flag_insdel_cluster_indel_h2 {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      hap = "h2",
      vartype="indel",
      inbed = call_pav.call_cigar_merge_h2.insdelBedMerge,
      threads = "4",
      mem_gb = "16",
      sample = sample
  }
  call call_inv.call_inv_flag_insdel_cluster_sv_h1 {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      hap = "h1",
      vartype="sv",
      inbed = call_pav.call_cigar_merge_h1.insdelBedMerge,
      threads = "4",
      mem_gb = "16",
      sample = sample
  }
  call call_inv.call_inv_flag_insdel_cluster_sv_h2 {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      hap = "h2",
      vartype="sv",
      inbed = call_pav.call_cigar_merge_h2.insdelBedMerge,
      threads = "4",
      mem_gb = "16",
      sample = sample
  }
  call call_pav.call_mappable_bed_h1 {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      hap = "h1",
      trimBed = align.align_cut_tig_overlap_h1.trimBed,
      delBed = call_lg.call_merge_lg_del_h1.mergeBed,
      insBed = call_lg.call_merge_lg_ins_h1.mergeBed,
      invBed = call_lg.call_merge_lg_inv_h1.mergeBed,
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  call call_pav.call_mappable_bed_h2 {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      hap = "h2",
      trimBed = align.align_cut_tig_overlap_h2.trimBed,
      delBed = call_lg.call_merge_lg_del_h2.mergeBed,
      insBed = call_lg.call_merge_lg_ins_h2.mergeBed,
      invBed = call_lg.call_merge_lg_inv_h2.mergeBed,
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  call call_inv.call_inv_merge_flagged_loci_h1 {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      hap = "h1",
      indelFlag = call_inv.call_inv_flag_insdel_cluster_indel_h1.bed,
      svFlag = call_inv.call_inv_flag_insdel_cluster_sv_h1.bed,
      snvCluster = call_inv.call_inv_cluster_snv_h1.bed,
      svCluster = call_inv.call_inv_cluster_sv_h1.bed,
      indelCluster = call_inv.call_inv_cluster_indel_h1.bed,
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  call call_inv.call_inv_merge_flagged_loci_h2 {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      hap = "h2",
      indelFlag = call_inv.call_inv_flag_insdel_cluster_indel_h2.bed,
      svFlag = call_inv.call_inv_flag_insdel_cluster_sv_h2.bed,
      snvCluster = call_inv.call_inv_cluster_snv_h2.bed,
      svCluster = call_inv.call_inv_cluster_sv_h2.bed,
      indelCluster = call_inv.call_inv_cluster_indel_h2.bed,
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  scatter(i in range(60)) {
     call call_inv.call_inv_batch_h1 {
      input:
        pav_conf = config,
        pav_sw = pav_tar,
        pav_asm = setup.tar_asm.asm_tar,
        hap = "h1",
        trimBed = align.align_cut_tig_overlap_h1.trimBed,
        flag = call_inv.call_inv_merge_flagged_loci_h1.bed,
        asmFai = align.align_get_tig_fa_h1.asmFai,
        asmGz = align.align_get_tig_fa_h1.asmGz,
        batch = i,
        threads = "8",
        mem_gb = "8",
        sample = sample
     }
  }
  scatter(i in range(60)) {
     call call_inv.call_inv_batch_h2 {
      input:
        pav_conf = config,
        pav_sw = pav_tar,
        pav_asm = setup.tar_asm.asm_tar,
        hap = "h2",
        trimBed = align.align_cut_tig_overlap_h2.trimBed,
        flag = call_inv.call_inv_merge_flagged_loci_h2.bed,
        asmFai = align.align_get_tig_fa_h2.asmFai,
        asmGz = align.align_get_tig_fa_h2.asmGz,
        batch = i,
        threads = "8",
        mem_gb = "8",
        sample = sample
     }
  }
  call call_inv.call_inv_batch_merge_h1 {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      hap = "h1",
      invBed = call_inv.call_inv_batch_h1.bed,
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  call call_inv.call_inv_batch_merge_h2 {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      hap = "h2",
      invBed = call_inv.call_inv_batch_h2.bed,
      threads = "1",
      mem_gb = "8",
      sample = sample
  }
  call call_pav.call_integrate_sources_h1 {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      hap = "h1",
      preInvSvBed = call_pav.call_cigar_merge_h1.insdelBedMerge,
      preInvSnvBed = call_pav.call_cigar_merge_h1.snvBedMerge,
      delBedIn = call_lg.call_merge_lg_del_h1.mergeBed,
      insBedIn = call_lg.call_merge_lg_ins_h1.mergeBed,
      invBedIn = call_lg.call_merge_lg_inv_h1.mergeBed,
      invBatch = call_inv.call_inv_batch_merge_h1.bed,
      threads = "2",
      mem_gb = "32",
      sample = sample
  }
  call call_pav.call_integrate_sources_h2 {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      hap = "h2",
      preInvSvBed = call_pav.call_cigar_merge_h2.insdelBedMerge,
      preInvSnvBed = call_pav.call_cigar_merge_h2.snvBedMerge,
      delBedIn = call_pav.call_merge_lg_del_h2.mergeBed,
      insBedIn = call_pav.call_merge_lg_ins_h2.mergeBed,
      invBedIn = call_lg.call_merge_lg_inv_h2.mergeBed,
      invBatch = call_inv.call_inv_batch_merge_h2.bed,
      threads = "2",
      mem_gb = "32",
      sample = sample
  }
  scatter(chrom in chroms) {
     call call_pav.call_merge_haplotypes_chrom_svindel_ins {
      input:
        pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
        svtype = "svindel_ins",
        insBed_h1 = call_pav.call_integrate_sources_h1.insBed,
        insBed_h2 = call_pav.call_integrate_sources_h2.insBed,
        callable_h1 = call_pav.call_mappable_bed_h2.bed,
        callable_h2 = call_pav.call_mappable_bed_h1.bed,
        chrom = chrom,
        threads = "8",
        mem_gb = "12",
        sample = sample
     }
  }
  scatter(chrom in chroms) {
     call call_pav.call_merge_haplotypes_chrom_svindel_del {
      input:
        pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
        svtype = "svindel_del",
        delBed_h1 = call_pav.call_integrate_sources_h1.delBed,
        delBed_h2 = call_pav.call_integrate_sources_h2.delBed,
        callable_h1 = call_pav.call_mappable_bed_h2.bed,
        callable_h2 = call_pav.call_mappable_bed_h1.bed,
        chrom = chrom,
        threads = "8",
        mem_gb = "12",
        sample = sample
     }
  }
  scatter(chrom in chroms) {
     call call_pav.call_merge_haplotypes_chrom_svinv {
      input:
        pav_conf = config,
        pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
        svtype = "sv_inv",
        invBed_h1 = call_pav.call_integrate_sources_h1.invBed,
        invBed_h2 = call_pav.call_integrate_sources_h2.invBed,
        callable_h1 = call_pav.call_mappable_bed_h2.bed,
        callable_h2 = call_pav.call_mappable_bed_h1.bed,
        chrom = chrom,
        threads = "8",
        mem_gb = "12",
        sample = sample
     }
  }
  scatter(chrom in chroms) {
     call call_pav.call_merge_haplotypes_chrom_snv {
      input:
        pav_conf = config,
        pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
        svtype = "snv_snv",
        snvBed_h1 = call_pav.call_integrate_sources_h1.snvBed,
        snvBed_h2 = call_pav.call_integrate_sources_h2.snvBed,
        callable_h1 = call_pav.call_mappable_bed_h2.bed,
        callable_h2 = call_pav.call_mappable_bed_h1.bed,
        chrom = chrom,
        threads = "8",
        mem_gb = "12",
        sample = sample
     }
  }
  call call_pav.call_merge_haplotypes_snv {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      svtype = "snv_snv",
      inbed = call_pav.call_merge_haplotypes_chrom_snv.bed,
      threads = "12",
      mem_gb = "24",
      sample = sample
  }
  call call_pav.call_merge_haplotypes_inv {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      svtype = "sv_inv",
      inbed = call_pav.call_merge_haplotypes_chrom_svinv.bed,
      threads = "12",
      mem_gb = "24",
      sample = sample
  }
  call call_pav.call_merge_haplotypes_svindel_ins {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      svtype = "sv_indel_ins",
      inbed = call_pav.call_merge_haplotypes_chrom_svindel_ins.bed,
      threads = "12",
      mem_gb = "24",
      sample = sample
  }
  call call_pav.call_merge_haplotypes_svindel_del {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      svtype = "sv_indel_del",
      inbed = call_pav.call_merge_haplotypes_chrom_svindel_del.bed,
      threads = "12",
      mem_gb = "24",
      sample = sample
  }
  call call_pav.call_final_bed {
    input:
      pav_conf = config,
      pav_sw = pav_tar,
      pav_asm = setup.tar_asm.asm_tar,
      invBed = call_pav.call_merge_haplotypes_inv.bed,
      insBed = call_pav.call_merge_haplotypes_svindel_ins.bed,
      delBed = call_pav.call_merge_haplotypes_svindel_del.bed,
      snvBed = call_pav.call_merge_haplotypes_snv.bed,
      threads = "4",
      mem_gb = "16",
      sample = sample
  }
  output {
    File snvBed = call_pav.call_final_bed.snvBedOut
    File invBed = call_pav.call_final_bed.invBedOut
    File svInsBed = call_pav.call_final_bed.svInsBed
    File svDelBed = call_pav.call_final_bed.svDelBed
    File indelInsBed = call_pav.call_final_bed.indelInsBed
    File indelDelBed = call_pav.call_final_bed.indelDelBed
    File invFasta = call_pav.call_final_bed.invFasta
    File svInsFasta = call_pav.call_final_bed.svInsFasta
    File svDelFasta = call_pav.call_final_bed.svDelFasta
    File indelInsFasta = call_pav.call_final_bed.indelInsFasta
    File indelDelFasta = call_pav.call_final_bed.indelDelFasta
  }
}
