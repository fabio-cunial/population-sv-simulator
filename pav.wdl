workflow pav {
   input {
     File ref
     File hapOne
     File hapTwo
     String sample
   }
   call align_ref {
   }
   call align_get_tig_fa {
      	input:
        	asm = hapOne,
        	hap = "h1"
   }
   call align_get_tig_fa {
      input:
         asm = hapTwo,
         hap = "h2"
   }
   call align_ref_anno_n_gap {
   		input:
   			refGz = align_ref.refGz,
   			refFai = align_ref.refFai
   }
   call align_map {
   		input:
   			refGz = align_ref.refGz,
   			hap = "h1",
   			asmGz = align_get_tig_fa.asmGz,
   			threads = "8"
   }
   call align_map {
   		input:
   			refGz = align_ref.refGz,
   			hap = "h2",
   			asmGz = align_get_tig_fa.asmGz,
   			threads = "8"
   }
   call align_get_read_bed {
   		input:
   			refFai = align_ref.refFai,
   			hap = "h1",
   			samGz = align_get_tig_fa.samGz,
   			threads = "1",
   }
   call align_get_read_bed {
   		input:
   			refFai = align_ref.refFai,
   			hap = "h2",
   			samGz = align_get_tig_fa.samGz,
   			threads = "1",
   }
   call align_cut_tig_overlap {
   		input:
   			hap = "h1",
   			bedGz = align_get_read_bed.bedGz,
   			asmFai = align_get_tig_fa.asmFai
   }
   call align_cut_tig_overlap {
   		input:
   			hap = "h2",
   			bedGz = align_get_read_bed.bedGz,
   			asmFai = align_get_tig_fa.asmFai
   }
   scatter(i in range(10)) {
	   call call_cigar {
	   		input:
		   		hap = "h1"
		   		trimBed = align_cut_tig_overlap.trimBed,
		   		asmGz = align_get_tig_fa.asmGz,
		   		batch = i,
		   		threads = "8"
	   }
	}
   scatter(i in range(10)) {
	   call call_cigar {
	   		input:
		   		hap = "h2",
		   		trimBed = align_cut_tig_overlap.trimBed,
		   		asmGz = align_get_tig_fa.asmGz,
		   		batch = i,
		   		threads = "8"
	   }
	}
   call call_lg_split {
   		input:
	   		hap = "h1",
	   		trimBed = align_cut_tig_overlap.trimBed,  			
   }
   call call_lg_split {
   		input:
	   		hap = "h2",
	   		trimBed = align_cut_tig_overlap.trimBed,  			
   }
   scatter(i in range(10)) {
	   call call_lg_discover {
	   		input:
		   		hap = "h2",
		   		trimBed = align_cut_tig_overlap.trimBed,
		   		gaps = align_ref_anno_n_gap.gaps,
		   		asmFai = align_get_tig_fa.asmFai,
		   		asmGz = align_get_tig_fa.asmGz,
		   		batch = call_lg_split.batch,
		   		threads = "8"
	   }
	}
   scatter(i in range(10)) {
	   call call_lg_discover {
	   		input:
		   		hap = "h2",
		   		trimBed = align_cut_tig_overlap.trimBed,
		   		gaps = align_ref_anno_n_gap.gaps,
		   		asmFai = align_get_tig_fa.asmFai,
		   		asmGz = align_get_tig_fa.asmGz,
		   		batch = call_lg_split.batch,
		   		threads = "8"
	   }
	}
   call call_cigar_merge {
   		input:
			threads = "8",
			hap = "h1",
			insdelBatch = call_cigar.insdelBed,
			snvBatch = call_cigar.snvBed,
   }
   call call_cigar_merge {
		input:
			threads = "8",
			hap = "h1",
			insdelBatch = call_cigar.insdelBed,
			snvBed = call_cigar.snvBed,
   }
   call call_merge_lg {
   		input:
   			hap = "h1",
   			svtype = "del",
   			bed = call_lg_discover.delBed
   }
   call call_merge_lg {
   		input:
   			hap = "h1",
   			svtype = "ins",
   			bed = call_lg_discover.insBed
   }
 	call call_merge_lg {
		input:
			hap = "h1",
			svtype = "inv",
			bed = call_lg_discover.invBed
   }
   call call_merge_lg {
   		input:
   			hap = "h2",
   			svtype = "del",
   			bed = call_lg_discover.delBed
   }
   call call_merge_lg {
   		input:
   			hap = "h2",
   			svtype = "ins",
   			bed = call_lg_discover.insBed
   }
   call call_merge_lg {
		input:
			hap = "h2",
			svtype = "inv",
			bed = call_lg_discover.invBed
   }
   call call_inv_cluster {
   		input:
   			hap = "h1",
   			vartype="indel",
   			bed = call_cigar_merge.insdelBedMerge,
   			threads = "1"   	
   }
   call call_inv_cluster {
   		input:
   			hap = "h1",
   			vartype="snv",
   			bed = call_cigar_merge.insdelBedMerge,
   			threads = "1"   	
   }
   call call_inv_cluster {
   		input:
   			hap = "h1",
   			vartype="sv",
   			bed = call_cigar_merge.insdelBedMerge
   			threads = "1"   	
   }
   call call_inv_cluster {
   		input:
   			hap = "h2",
   			vartype="indel",
   			bed = call_cigar_merge.insdelBedMerge,
   			threads = "1"   	
   }
   call call_inv_cluster {
   		input:
   			hap = "h2",
   			vartype="snv",
   			bed = call_cigar_merge.insdelBedMerge,
   			threads = "1"   	
   }
   call call_inv_cluster {
   		input:
   			hap = "h2",
   			vartype="sv",
   			bed = call_cigar_merge.insdelBedMerge
   			threads = "1"   	
   }
   call call_inv_flag_insdel_cluster {
   		input:
   			hap = "h1",
   			vartype="indel",
   			bed = call_cigar_merge.insdelBedMerge
   }
   call call_inv_flag_insdel_cluster {
   		input:
   			hap = "h2",
   			vartype="indel",
   			bed = call_cigar_merge.insdelBedMerge
   }
  call call_inv_flag_insdel_cluster {
   		input:
   			hap = "h1",
   			vartype="sv",
   			bed = call_cigar_merge.insdelBedMerge
   }
   call call_inv_flag_insdel_cluster {
   		input:
   			hap = "h2",
   			vartype="sv",
   			bed = call_cigar_merge.insdelBedMerge
   }
   call call_mappable_bed {

   }
   call call_inv_merge_flagged_loci {

   }
   call call_inv_batch {

   }
   call call_inv_batch_merge {

   }
   call call_integrate_sources {

   }
   call call_merge_haplotypes_chrom {

   }
   call call_merge_haplotypes {

   }
   call call_final_bed {

   }



  output {
     File outfile = ReadItBackToMe.repeated_name
  }
}

task align_ref {
	input {
		File ref
	}
	command {
		snakemake --threads ${threads} data/ref/ref.fa.gz data/ref/ref.fa.gz.fai
	}
	output {
		File refGz = "data/ref/ref.fa.gz"
		File refFai = "data/ref/ref.fa.gz.fai"
	}
	runtime {
		docker : ""
	}
}

task align_get_tig_fa {
	input {
		File asm
		String sample
		String hap
		String threads
	}
	command {
		snakemake --threads ${threads} temp/${sample}/align/contigs_${hap}.fa.gz temp/${sample}/align/contigs_${hap}.fa.gz.fai
	}
	output {
		File asmGz = "temp/${sample}/align/contigs_${hap}.fa.gz"
		File asmFai = "temp/${sample}/align/contigs_${hap}.fa.gz.fai"
	}
	runtime {
		docker : ""
	}
}

task align_ref_anno_n_gap {
	input {
		File refGZ
		File refFai
		String threads
	}
	command {
		snakemake --threads ${threads} data/ref/n_gap.bed.gz
	}
	output {
		File gaps = "data/ref/n_gap.bed.gz"
	}
}

task align_map {
	input {
		String hap
		String sample
		File refGz
		File asmGz
		String threads
	}
	command {
		snakemake --threads ${threads} temp/${sample}/align/pre-cut/aligned_tig_${hap}.sam.gz
	}
	output {
		File  samGz = "temp/${sample}/align/pre-cut/aligned_tig_${hap}.sam.gz"
	}
}

task align_get_read_bed {
	input {
		File refGz
		File samGz
		String hap
		String threads
	}
	command {
		snakemake --threads ${threads} results/${sample}/align/pre-cut/aligned_tig_${hap}.bed.gz
	}
	output {
		File bedGz = results/${sample}/align/pre-cut/aligned_tig_${hap}.bed.gz
		File headerGz = results/${sample}/align/pre-cut/aligned_tig_${hap}.headers.gz
	}
}

task align_cut_tig_overlap {
	input {
		File asmFai
		File bedGz
		String threads
	}
	command {
		snakemake --threads ${threads} results/${sample}/align/aligned_tig_${hap}.bed.gz
	}
	output {
		File trimBed = "results/${sample}/align/aligned_tig_${hap}.bed.gz"
	}
}

task call_cigar {
	input {
		File trimBed,
		File asmGz,
		String batch,
		String threads,
		String hap
	}
	command {
		snakemake --threads ${threads} temp/${sample}/cigar/batched/insdel_${hap}_${batch}.bed.gz temp/${sample}/cigar/batched/snv.bed_${hap}_${batch}.gz
	}
	output {
		File snvBed = "temp/${sample}/cigar/batched/snv.bed_${hap}_${batch}.gz"
		File insdelBed = "temp/${sample}/cigar/batched/insdel_${hap}_${batch}.bed.gz"
	}
}

task call_lg_split {
	input {
		File trimBed
		String hap
	}
	command {
		snakemake --threads ${threads} temp/${sample}/lg_sv/batch_${hap}.tsv.gz
	}
	output {
		File batch = "temp/${sample}/lg_sv/batch_${hap}.tsv.gz"
	}
}

task call_lg_discover {
	input {
		File trimBed
		File batch
		File asmGz
		File asmFai
		File gaps
		String hap
		String threads
		String batch
	}
	command {
		snakemake --threads ${threads} temp/${sample}/lg_sv/batch/sv_ins_${hap}_${batch}.bed.gz temp/${sample}/lg_sv/batch/sv_del_${hap}_${batch}.bed.gz temp/${sample}/lg_sv/batch/sv_inv_${hap}_${batch}.bed.gz
	}
	output {
		File insBed = "temp/${sample}/lg_sv/batch/sv_ins_${hap}_${batch}.bed.gz"
		File delBed = "temp/${sample}/lg_sv/batch/sv_del_${hap}_${batch}.bed.gz"
		File invBed = "temp/${sample}/lg_sv/batch/sv_inv_${hap}_${batch}.bed.gz"
	}
}

task call_cigar_merge {
	input {
		String hap
		File insdelBatch
		File snvBatch
		String threads
	}
	command {
		snakemake --threads ${threads} temp/${sample}/cigar/pre_inv/svindel_insdel_${hap}.bed.gz temp/${sample}/cigar/pre_inv/snv_snv_${hap}.bed.gz
	}
	output {
		File insdelBedMerge = "temp/${sample}/cigar/pre_inv/svindel_insdel_${hap}.bed.gz"
		File snvBedMerge = "temp/${sample}/cigar/pre_inv/snv_snv_${hap}.bed.gz"
	}
}

task call_merge_lg {
	input {
		File bed
		String hap
		String threads
		String svtype
	}
	command {
		snakemake --threads ${threads} temp/${sample}/lg_sv/sv_${svtype}_${hap}.bed.gz
	}
	output {
		File mergeBed = "temp/${sample}/lg_sv/sv_${svtype}_${hap}.bed.gz"
	}
}

task call_inv_flag_insdel_cluster {
	input {
		File bed
		String hap
		String threads
		String vartype
	}
	command {
		snakemake --threads ${threads} temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz
	}
	output {
		File bed = "temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz"
	}
}

task call_inv_cluster {
	input {
		File bed
		String threads
		String hap
		String vartype
	}
	command {
		snakemake --threads ${threads} temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz
	}
	output {
		File bed = "temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz"
	}
}