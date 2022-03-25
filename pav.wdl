workflow pav {
	input {
		File ref
		File hapOne
		File hapTwo
		String sample
		Array[String] chroms = read_tsv(refFai)
	}
	call align_ref {
		input:
			threads = "1",
			mem_gb = "8"
	}
	call align_get_tig_fa_h1 {
   	input:
	     	asm = hapOne,
	     	hap = "h1",
			threads = "1",
			mem_gb = "8"
	}
	call align_get_tig_fa_h2 {
	   input:
	      asm = hapTwo,
	      hap = "h2",
			threads = "1",
			mem_gb = "8"
	}
	call align_ref_anno_n_gap {
			input:
				refGz = align_ref.refGz,
				refFai = align_ref.refFai,
				threads = "1",
				mem_gb = "8"
	}
	call align_map_h1 {
			input:
				refGz = align_ref.refGz,
				hap = "h1",
				asmGz = align_get_tig_fa_h2.asmGz,
				threads = "8",
				mem_gb = "12"
	}
	call align_map_h2 {
			input:
				refGz = align_ref.refGz,
				hap = "h2",
				asmGz = align_get_tig_fa_h2.asmGz,
				threads = "8",
				mem_gb = "12"
	}
	call align_get_read_bed_h1 {
		input:
			refFai = align_ref.refFai,
			hap = "h1",
			samGz = align_get_tig_fa_h1.samGz,
			threads = "1",
			mem_gb = "32"
	}
	call align_get_read_bed_h2 {
		input:
			refFai = align_ref.refFai,
			hap = "h2",
			samGz = align_get_tig_fa_h2.samGz,
			threads = "1",
			mem_gb = "32"
	}
	call align_cut_tig_overlap_h1 {
		input:
			hap = "h1",
			bedGz = align_get_read_bed_h1.bedGz,
			asmFai = align_get_tig_fa_h1.asmFai,
			threads = "1",
			mem_gb = "8"
	}
	call align_cut_tig_overlap_h2 {
		input:
			hap = "h2",
			bedGz = align_get_read_bed_h2.bedGz,
			asmFai = align_get_tig_fa_h2.asmFai,
			threads = "1",
			mem_gb = "8"
	}
	scatter(i in range(10)) {
	   call call_cigar_h1 {
	   		input:
		   		hap = "h1"
		   		trimBed = align_cut_tig_overlap_h1.trimBed,
		   		asmGz = align_get_tig_fa_h1.asmGz,
		   		batch = i,
		   		threads = "1"
		   		mem_gb = "24"
	   }
	}
	scatter(i in range(10)) {
	   call call_cigar_h2 {
	   		input:
		   		hap = "h2",
		   		trimBed = align_cut_tig_overlap_h2.trimBed,
		   		asmGz = align_get_tig_fa_h2.asmGz,
		   		batch = i,
		   		threads = "1"
		   		mem_gb = "24"
	   }
	}
	call call_lg_split_h1 {
		input:
   		hap = "h1",
   		trimBed = align_cut_tig_overlap_h2.trimBed,
			threads = "1",
			mem_gb = "8"
	}
	call call_lg_split_h2 {
		input:
   		hap = "h2",
   		trimBed = align_cut_tig_overlap_h2.trimBed,
			threads = "1",
			mem_gb = "8"
	}
	scatter(i in range(10)) {
	   call call_lg_discover_h1 {
	   		input:
		   		hap = "h1",
		   		trimBed = align_cut_tig_overlap_h1.trimBed,
		   		gaps = align_ref_anno_n_gap.gaps,
		   		asmFai = align_get_tig_fa_h1.asmFai,
		   		asmGz = align_get_tig_fa_h1.asmGz,
		   		batch = call_lg_split_h1.batch,
		   		threads = "8",
					mem_gb = "8"
	   }
	}
	scatter(i in range(10)) {
	   call call_lg_discover_h2 {
	   		input:
		   		hap = "h2",
		   		trimBed = align_cut_tig_overlap_h2.trimBed,
		   		gaps = align_ref_anno_n_gap.gaps,
		   		asmFai = align_get_tig_fa_h2.asmFai,
		   		asmGz = align_get_tig_fa_h2.asmGz,
		   		batch = call_lg_split_h2.batch,
		   		threads = "8",
					mem_gb = "8"
	   }
	}
	call call_cigar_merge_h1 {
		input:
			hap = "h1",
			insdelBatch = call_cigar_h1.insdelBed,
			snvBatch = call_cigar_h1.snvBed,
			threads = "1",
			mem_gb = "8"
	}
	call call_cigar_merge_h2 {
		input:
			hap = "h2",
			insdelBatch = call_cigar_h2.insdelBed,
			snvBed = call_cigar_h2.snvBed,
			threads = "1",
			mem_gb = "8"
	}
	call call_merge_lg_del_h1 {
		input:
			hap = "h1",
			svtype = "del",
			bed = call_lg_discover_h1.delBed,
			threads = "1",
			mem_gb = "8"
	}
	call call_merge_lg_ins_h1 {
		input:
			hap = "h1",
			svtype = "ins",
			bed = call_lg_discover_h1.insBed,
			threads = "1",
			mem_gb = "8"
	}
	call call_merge_lg_inv_h1 {
		input:
			hap = "h1",
			svtype = "inv",
			bed = call_lg_discover_h1.invBed,
			threads = "1",
			mem_gb = "8"
	}
	call call_merge_lg_del_h2 {
		input:
			hap = "h2",
			svtype = "del",
			bed = call_lg_discover_h2.delBed,
			threads = "1",
			mem_gb = "8"
	}
	call call_merge_lg_ins_h2 {
			input:
				hap = "h2",
				svtype = "ins",
				bed = call_lg_discover_h2.insBed,
				threads = "1",
				mem_gb = "8"
	}
	call call_merge_lg_inv_h2 {
		input:
			hap = "h2",
			svtype = "inv",
			bed = call_lg_discover_h2.invBed,
			threads = "1",
			mem_gb = "8"
	}
	call call_inv_cluster_indel_h1 {
		input:
			hap = "h1",
			vartype="indel",
			bed = call_cigar_merge_h1.insdelBedMerge,
			threads = "4",
			mem_gb = "16"
	}
	call call_inv_cluster_snv_h1 {
		input:
			hap = "h1",
			vartype="snv",
			bed = call_cigar_merge_h1.insdelBedMerge,
			threads = "4",
			mem_gb = "16"
	}
	call call_inv_cluster_sv_h1 {
		input:
			hap = "h1",
			vartype="sv",
			bed = call_cigar_merge_h1.insdelBedMerge,
			threads = "4",
			mem_gb = "16"
	}
	call call_inv_cluster_indel_h2 {
		input:
			hap = "h2",
			vartype="indel",
			bed = call_cigar_merge_h2.insdelBedMerge,
			threads = "4",
			mem_gb = "16"
	}
	call call_inv_cluster_snv_h2 {
		input:
			hap = "h2",
			vartype="snv",
			bed = call_cigar_merge_h2.insdelBedMerge,
			threads = "4",
			mem_gb = "16"
	}
	call call_inv_cluster_sv_h2 {
		input:
			hap = "h2",
			vartype="sv",
			bed = call_cigar_merge_h2.insdelBedMerge
			threads = "4",
			mem_gb = "16"
	}
	call call_inv_flag_insdel_cluster_h1 {
		input:
			hap = "h1",
			vartype="indel",
			bed = call_cigar_merge_h1.insdelBedMerge,
			threads = "4",
			mem_gb = "16"
	}
	call call_inv_flag_insdel_cluster_h2 {
		input:
			hap = "h2",
			vartype="indel",
			bed = call_cigar_merge_h2.insdelBedMerge,
			threads = "4",
			mem_gb = "16"
	}
	call call_inv_flag_insdel_cluster_sv_h1 {
		input:
			hap = "h1",
			vartype="sv",
			bed = call_cigar_merge_h1.insdelBedMerge,
			threads = "4",
			mem_gb = "16"
	}
	call call_inv_flag_insdel_cluster_sv_h2 {
		input:
			hap = "h2",
			vartype="sv",
			bed = call_cigar_merge_h2.insdelBedMerge,
			threads = "4",
			mem_gb = "16"
	}
	call call_mappable_bed_h1 {
		input:
			hap = "h1",
			trimBed = align_cut_tig_overlap_h1.trimBed,
			delBed = call_merge_lg_del_h1.mergeBed,
			insBed = call_merge_lg_ins_h1.mergeBed,
			invBed = call_merge_lg_inv_h1.mergeBed,
			threads = "1",
			mem_gb = "8"
	}
	call call_mappable_bed_h2 {
		input:
			hap = "h2",
			trimBed = align_cut_tig_overlap_h2.trimBed,
			delBed = call_merge_lg_del_h2.mergeBed,
			insBed = call_merge_lg_ins_h2.mergeBed,
			invBed = call_merge_lg_inv_h2.mergeBed,
			threads = "1",
			mem_gb = "8"
	}
	call call_inv_merge_flagged_loci_h1 {
		input:
			hap = "h1"
			indelFlag = call_inv_flag_insdel_cluster_indel_h1.bed,
			svFlag = call_inv_flag_insdel_cluster_sv_h1.bed,
			snvCluster = call_inv_cluster_snv_h1.bed,
			svCluster = call_inv_cluster_sv_h1.bed,
			indelCluster = call_inv_cluster_indel_h1.bed,
			threads = "1",
			mem_gb = "8"
	}
	call call_inv_merge_flagged_loci_h2 {
		input:
			hap = "h2"
			indelFlag = call_inv_flag_insdel_cluster_indel_h2.bed,
			svFlag = call_inv_flag_insdel_cluster_sv_h2.bed,
			snvCluster = call_inv_cluster_snv_h2.bed,
			svCluster = call_inv_cluster_sv_h2.bed,
			indelCluster = call_inv_cluster_indel_h2.bed,
			threads = "1",
			mem_gb = "8"
	}
	scatter(i in range(60)) {
	   call call_inv_batch_h1 {
	   		input:
		   		hap = "h1",
		   		trimBed = align_cut_tig_overlap_h1.trimBed,
		   		flag = call_inv_merge_flagged_loci_h1.bed,
		   		asmFai = align_get_tig_fa_h1.asmFai,
		   		asmGz = align_get_tig_fa_h1.asmGz,
		   		batch = i,
					threads = "8",
					mem_gb = "8"
	   }
	}
	scatter(i in range(60)) {
	   call call_inv_batch_h2 {
	   		input:
		   		hap = "h2",
		   		trimBed = align_cut_tig_overlap_h2.trimBed,
		   		flag = call_inv_merge_flagged_loci_h2.bed,
		   		asmFai = align_get_tig_fa_h2.asmFai,
		   		asmGz = align_get_tig_fa_h2.asmGz,
		   		batch = i,
					threads = "8",
					mem_gb = "8"
	   }
	}
	call call_inv_batch_merge_h1 {
		input:
			hap = "h1"
			invBed = call_inv_batch_h1.bed,
			threads = "1",
			mem_gb = "8"
	}
	call call_inv_batch_merge_h2 {
		input:
			hap = "h2"
			invBed = call_inv_batch_h2.bed,
			threads = "1",
			mem_gb = "8"
	}
	call call_integrate_sources_h1 {
			input:
				hap = "h1"
				preInvSvBed = call_cigar_merge.insdelBedMerge,
				preInvSnvBed = call_cigar_merge.snvBedMerge,
	   		delBed = call_merge_lg_del.mergeBed,
	   		insBed = call_merge_lg_ins.mergeBed,
	   		invBed = call_merge_lg_inv.mergeBed,
	   		invBatch = call_inv_batch_merge.bed,
	   		threads = "2",
	   		mem_gb = "32"
	}
	call call_integrate_sources_h2 {
			input:
				hap = "h2"
				preInvSvBed = call_cigar_merge.insdelBedMerge,
				preInvSnvBed = call_cigar_merge.snvBedMerge,
	   		delBed = call_merge_lg_del.mergeBed,
	   		insBed = call_merge_lg_ins.mergeBed,
	   		invBed = call_merge_lg_inv.mergeBed,
	   		invBatch = call_inv_batch_merge.bed,
	   		threads = "2",
	   		mem_gb = "32"
	}
	scatter(chrom in chroms) {
	   call call_merge_haplotypes_chrom_svindel_ins {
	   		input:
	   			svtype = "svindel_ins",
		   		insBed_h1 = call_integrate_sources_h1.insBed,
		   		insBed_h2 = call_integrate_sources_h2.insBed,
		   		callable_h1 = call_mappable_bed_h2.bed,
		   		callable_h2 = call_mappable_bed_h1.bed,
		   		chrom = chrom,
		   		threads = "8",
		   		mem_gb = "12"
	   }
	}
	scatter(chrom in chroms) {
	   call call_merge_haplotypes_chrom_svindel_del {
	   		input:
	   			svtype = "svindel_del",
		   		delBed_h1 = call_integrate_sources_h1.delBed,
		   		delBed_h2 = call_integrate_sources_h2.delBed,
		   		callable_h1 = call_mappable_bed_h2.bed,
		   		callable_h2 = call_mappable_bed_h1.bed,
		   		chrom = chrom,
		   		threads = "8",
		   		mem_gb = "12"
	   }
	}
	scatter(chrom in chroms) {
	   call call_merge_haplotypes_chrom_svinv {
	   		input:
	   			svtype = "sv_inv",
		   		invBed_h1 = call_integrate_sources_h1.invBed,
		   		invBed_h2 = call_integrate_sources_h2.invBed,
		   		callable_h1 = call_mappable_bed_h2.bed,
		   		callable_h2 = call_mappable_bed_h1.bed,
		   		chrom = chrom,
		   		threads = "8",
		   		mem_gb = "12"
	   }
	}
	scatter(chrom in chroms) {
	   call call_merge_haplotypes_chrom_snv {
	   		input:
	   			svtype = "snv_snv",
		   		snvBed_h1 = call_integrate_sources_h1.snvBed,
		   		snvBed_h2 = call_integrate_sources_h2.snvBed,
		   		callable_h1 = call_mappable_bed_h2.bed,
		   		callable_h2 = call_mappable_bed_h1.bed,
		   		chrom = chrom,
		   		threads = "8",
		   		mem_gb = "12"
	   }
	}
	call call_merge_haplotypes_snv {
		input:
			svtype = "snv_snv",
			bed = call_merge_haplotypes_snv.bed,
			threads = "12",
			mem_gb = "24"
	}
	call call_merge_haplotypes_inv {
		input:
			svtype = "sv_inv",
			bed = call_merge_haplotypes_inv.bed,
			threads = "12",
			mem_gb = "24"
	}
	call call_merge_haplotypes_sv_indel_ins {
		input:
			svtype = "sv_indel_ins",
			bed = call_merge_haplotypes_svindel_ins.bed,
			threads = "12",
			mem_gb = "24"
	}
	call call_merge_haplotypes_sv_indel_del {
		input:
			svtype = "sv_indel_del",
			bed = call_merge_haplotypes_svindel_del.bed,
			threads = "12",
			mem_gb = "24"
	}
	call call_final_bed {
		input:
			invBed = call_merge_haplotypes_inv.bed,
			insBed = call_merge_haplotypes_svindel_ins.bed,
			delBed = call_merge_haplotypes_svindel_del.bed,
			snvBed = call_merge_haplotypes_snv.bed,
			threads = "4",
			mem_gb = "16"
	}
	output {
		File snvBed = call_final_bed.snvBed,
		File invBed = call_final_bed.invBed,
		File svInsBed = call_final_bed.svInsBed,
		File svDelBed = call_final_bed.svDelBed,
		File indelInsBed = call_final_bed.indelInsBed,
		File indelDelBed = call_final_bed.indelDelBed,
		File invFasta = call_final_bed.invFasta,
		File svInsFasta = call_final_bed.svInsFasta,
		File svDelFasta = call_final_bed.svDelFasta,
		File indelInsFasta = call_final_bed.indelInsFasta,
		File indelDelFasta = call_final_bed.indelDelFasta
	}
}

task align_ref {
	input {
		File ref
		String threads
		String mem_gb
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

task align_get_tig_fa_h1 {
	input {
		File asm
		String sample
		String hap
		String threads
		String mem_gb
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

task align_get_tig_fa_h2 {
	input {
		File asm
		String sample
		String hap
		String threads
		String mem_gb
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
		String mem_gb
	}
	command {
		snakemake --threads ${threads} data/ref/n_gap.bed.gz
	}
	output {
		File gaps = "data/ref/n_gap.bed.gz"
	}
}

task align_map_h1 {
	input {
		String hap
		String sample
		File refGz
		File asmGz
		String threads
		String mem_gb
	}
	command {
		snakemake --threads ${threads} temp/${sample}/align/pre-cut/aligned_tig_${hap}.sam.gz
	}
	output {
		File  samGz = "temp/${sample}/align/pre-cut/aligned_tig_${hap}.sam.gz"
	}
}

task align_map_h2 {
	input {
		String hap
		String sample
		File refGz
		File asmGz
		String threads
		String mem_gb
	}
	command {
		snakemake --threads ${threads} temp/${sample}/align/pre-cut/aligned_tig_${hap}.sam.gz
	}
	output {
		File  samGz = "temp/${sample}/align/pre-cut/aligned_tig_${hap}.sam.gz"
	}
}

task align_get_read_bed_h1 {
	input {
		File refGz
		File samGz
		String hap
		String threads
		String mem_gb
	}
	command {
		snakemake --threads ${threads} results/${sample}/align/pre-cut/aligned_tig_${hap}.bed.gz
	}
	output {
		File bedGz = results/${sample}/align/pre-cut/aligned_tig_${hap}.bed.gz
		File headerGz = results/${sample}/align/pre-cut/aligned_tig_${hap}.headers.gz
	}
}

task align_get_read_bed_h2 {
	input {
		File refGz
		File samGz
		String hap
		String threads
		String mem_gb
	}
	command {
		snakemake --threads ${threads} results/${sample}/align/pre-cut/aligned_tig_${hap}.bed.gz
	}
	output {
		File bedGz = results/${sample}/align/pre-cut/aligned_tig_${hap}.bed.gz
		File headerGz = results/${sample}/align/pre-cut/aligned_tig_${hap}.headers.gz
	}
}

task align_cut_tig_overlap_h1 {
	input {
		File asmFai
		File bedGz
		String threads
		String mem_gb
	}
	command {
		snakemake --threads ${threads} results/${sample}/align/aligned_tig_${hap}.bed.gz
	}
	output {
		File trimBed = "results/${sample}/align/aligned_tig_${hap}.bed.gz"
	}
}

task align_cut_tig_overlap_h2 {
	input {
		File asmFai
		File bedGz
		String threads
		String mem_gb
	}
	command {
		snakemake --threads ${threads} results/${sample}/align/aligned_tig_${hap}.bed.gz
	}
	output {
		File trimBed = "results/${sample}/align/aligned_tig_${hap}.bed.gz"
	}
}

task call_cigar_h1 {
	input {
		File trimBed,
		File asmGz,
		String batch,
		String threads
		String mem_gb,
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

task call_cigar_h2 {
	input {
		File trimBed,
		File asmGz,
		String batch,
		String threads
		String mem_gb,
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

task call_lg_split_h1 {
	input {
		File trimBed
		String hap
		String threads
		String mem_gb
	}
	command {
		snakemake --threads ${threads} temp/${sample}/lg_sv/batch_${hap}.tsv.gz
	}
	output {
		File batch = "temp/${sample}/lg_sv/batch_${hap}.tsv.gz"
	}
}

task call_lg_split_h2 {
	input {
		File trimBed
		String hap
		String threads
		String mem_gb
	}
	command {
		snakemake --threads ${threads} temp/${sample}/lg_sv/batch_${hap}.tsv.gz
	}
	output {
		File batch = "temp/${sample}/lg_sv/batch_${hap}.tsv.gz"
	}
}

task call_lg_discover_h1 {
	input {
		File trimBed
		File batch
		File asmGz
		File asmFai
		File gaps
		String hap
		String threads
		String mem_gb
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

task call_lg_discover_h2 {
	input {
		File trimBed
		File batch
		File asmGz
		File asmFai
		File gaps
		String hap
		String threads
		String mem_gb
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

task call_cigar_merge_h1 {
	input {
		String hap
		File insdelBatch
		File snvBatch
		String threads
		String mem_gb
	}
	command {
		snakemake --threads ${threads} temp/${sample}/cigar/pre_inv/svindel_insdel_${hap}.bed.gz temp/${sample}/cigar/pre_inv/snv_snv_${hap}.bed.gz
	}
	output {
		File insdelBedMerge = "temp/${sample}/cigar/pre_inv/svindel_insdel_${hap}.bed.gz"
		File snvBedMerge = "temp/${sample}/cigar/pre_inv/snv_snv_${hap}.bed.gz"
	}
}

task call_cigar_merge_h2 {
	input {
		String hap
		File insdelBatch
		File snvBatch
		String threads
		String mem_gb
	}
	command {
		snakemake --threads ${threads} temp/${sample}/cigar/pre_inv/svindel_insdel_${hap}.bed.gz temp/${sample}/cigar/pre_inv/snv_snv_${hap}.bed.gz
	}
	output {
		File insdelBedMerge = "temp/${sample}/cigar/pre_inv/svindel_insdel_${hap}.bed.gz"
		File snvBedMerge = "temp/${sample}/cigar/pre_inv/snv_snv_${hap}.bed.gz"
	}
}

task call_merge_lg_del_h1 {
	input {
		File bed
		String hap
		String threads
		String mem_gb
		String svtype
	}
	command {
		snakemake --threads ${threads} temp/${sample}/lg_sv/sv_${svtype}_${hap}.bed.gz
	}
	output {
		File mergeBed = "temp/${sample}/lg_sv/sv_${svtype}_${hap}.bed.gz"
	}
}

task call_merge_lg_ins_h1 {
	input {
		File bed
		String hap
		String threads
		String mem_gb
		String svtype
	}
	command {
		snakemake --threads ${threads} temp/${sample}/lg_sv/sv_${svtype}_${hap}.bed.gz
	}
	output {
		File mergeBed = "temp/${sample}/lg_sv/sv_${svtype}_${hap}.bed.gz"
	}
}

task call_merge_lg_inv_h1 {
	input {
		File bed
		String hap
		String threads
		String mem_gb
		String svtype
	}
	command {
		snakemake --threads ${threads} temp/${sample}/lg_sv/sv_${svtype}_${hap}.bed.gz
	}
	output {
		File mergeBed = "temp/${sample}/lg_sv/sv_${svtype}_${hap}.bed.gz"
	}
}

task call_merge_lg_del_h2 {
	input {
		File bed
		String hap
		String threads
		String mem_gb
		String svtype
	}
	command {
		snakemake --threads ${threads} temp/${sample}/lg_sv/sv_${svtype}_${hap}.bed.gz
	}
	output {
		File mergeBed = "temp/${sample}/lg_sv/sv_${svtype}_${hap}.bed.gz"
	}
}

task call_merge_lg_ins_h2 {
	input {
		File bed
		String hap
		String threads
		String mem_gb
		String svtype
	}
	command {
		snakemake --threads ${threads} temp/${sample}/lg_sv/sv_${svtype}_${hap}.bed.gz
	}
	output {
		File mergeBed = "temp/${sample}/lg_sv/sv_${svtype}_${hap}.bed.gz"
	}
}

task call_merge_lg_inv_h2 {
	input {
		File bed
		String hap
		String threads
		String mem_gb
		String svtype
	}
	command {
		snakemake --threads ${threads} temp/${sample}/lg_sv/sv_${svtype}_${hap}.bed.gz
	}
	output {
		File mergeBed = "temp/${sample}/lg_sv/sv_${svtype}_${hap}.bed.gz"
	}
}

task call_inv_flag_insdel_cluster_indel_h1 {
	input {
		File bed
		String hap
		String threads
		String mem_gb
		String vartype
	}
	command {
		snakemake --threads ${threads} temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz
	}
	output {
		File bed = "temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz"
	}
}

task call_inv_flag_insdel_cluster_sv_h1 {
	input {
		File bed
		String hap
		String threads
		String mem_gb
		String vartype
	}
	command {
		snakemake --threads ${threads} temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz
	}
	output {
		File bed = "temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz"
	}
}

task call_inv_flag_insdel_cluster_indel_h2 {
	input {
		File bed
		String hap
		String threads
		String mem_gb
		String vartype
	}
	command {
		snakemake --threads ${threads} temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz
	}
	output {
		File bed = "temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz"
	}
}

task call_inv_flag_insdel_cluster_sv_h2 {
	input {
		File bed
		String hap
		String threads
		String mem_gb
		String vartype
	}
	command {
		snakemake --threads ${threads} temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz
	}
	output {
		File bed = "temp/${sample}/inv_caller/flag/insdel_${vartype}_${hap}.bed.gz"
	}
}

task call_inv_cluster_sv_h1 {
	input {
		File bed
		String threads
		String mem_gb
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

task call_inv_cluster_indel_h1 {
	input {
		File bed
		String threads
		String mem_gb
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

task call_inv_cluster_snv_h1 {
	input {
		File bed
		String threads
		String mem_gb
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

task call_inv_cluster_sv_h2 {
	input {
		File bed
		String threads
		String mem_gb
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

task call_inv_cluster_indel_h2 {
	input {
		File bed
		String threads
		String mem_gb
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

task call_inv_cluster_snv_h2 {
	input {
		File bed
		String threads
		String mem_gb
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

task call_mappable_bed_h1 {
	input {
		File delBed
		File insBed
		File invBed
		File trimBed
		String hap
		String threads
		String mem_gb
	}
	command {
		snakemake --threads ${threads} results/${sample}/callable/callable_regions_${hap}_500.bed.gz
	}
	output {
		File bed = "results/${sample}/callable/callable_regions_${hap}_500.bed.gz"
	}
}

task call_mappable_bed_h2 {
	input {
		File delBed
		File insBed
		File invBed
		File trimBed
		String hap
		String threads
		String mem_gb
	}
	command {
		snakemake --threads ${threads} results/${sample}/callable/callable_regions_${hap}_500.bed.gz
	}
	output {
		File bed = "results/${sample}/callable/callable_regions_${hap}_500.bed.gz"
	}
}

task call_inv_merge_flagged_loci_h1 {
	input {
   		String hap
   		File indelFlag
   		File svFlag
   		File snvCluster
   		File svCluster
   		File indelCluster
   		String threads
		String mem_gb
   }
   command {
   	snakemake --threads results/${sample}/inv_caller/flagged_regions_${hap}.bed.gz
   }
   output {
   	File bed = "results/${sample}/inv_caller/flagged_regions_${hap}.bed.gz"
   }
}

task call_inv_merge_flagged_loci_h2 {
	input {
   		String hap
   		File indelFlag
   		File svFlag
   		File snvCluster
   		File svCluster
   		File indelCluster
   		String threads
		String mem_gb
   }
   command {
   	snakemake --threads results/${sample}/inv_caller/flagged_regions_${hap}.bed.gz
   }
   output {
   	File bed = "results/${sample}/inv_caller/flagged_regions_${hap}.bed.gz"
   }
}

task call_inv_batch_h1 {
	input {
		String hap
		File trimBed
		File flag
		File asmFai
		File asmGz
		String batch
		String threads
		String mem_gb
	}
	command {
		snakemake --threads ${threads} temp/${sample}/inv_caller/batch/${hap}/inv_call_${batch}.bed.gz
	}
	output {
		bed = "temp/${sample}/inv_caller/batch/${hap}/inv_call_${batch}.bed.gz"
	}
}

task call_inv_batch_h2 {
	input {
		String hap
		File trimBed
		File flag
		File asmFai
		File asmGz
		String batch
		String threads
		String mem_gb
	}
	command {
		snakemake --threads ${threads} temp/${sample}/inv_caller/batch/${hap}/inv_call_${batch}.bed.gz
	}
	output {
		bed = "temp/${sample}/inv_caller/batch/${hap}/inv_call_${batch}.bed.gz"
	}
}

task call_inv_batch_merge_h1 {
	input {
		String hap
		File invBed
		String threads
		String mem_gb
	}
	command {
		snakemake --threads ${threads} temp/${sample}/inv_caller/sv_inv_${hap}.bed.gz
	}
	output {
		File bed = "temp/${sample}/inv_caller/sv_inv_${hap}.bed.gz"
	}
}

task call_inv_batch_merge_h2 {
	input {
		String hap
		File invBed
		String threads
		String mem_gb
	}
	command {
		snakemake --threads ${threads} temp/${sample}/inv_caller/sv_inv_${hap}.bed.gz
	}
	output {
		File bed = "temp/${sample}/inv_caller/sv_inv_${hap}.bed.gz"
	}
}

task call_integrate_sources_h1 {
	input {
		String hap
		File preInvSvBed
		File preInvSnvBed
		File delBed
		File insBed
		File invBed
		File invBatch
		String threads
		String mem_gb
	}
	command {
		snakemake --threads ${threads} temp/${sample}/bed/integrated/${hap}/svindel_ins.bed.gz temp/${sample}/bed/integrated/${hap}/svindel_del.bed.gz temp/${sample}/bed/integrated/${hap}/snv_snv.bed.gz temp/${sample}/bed/integrated/${hap}/svindel_inv.bed.gz
	}
	output {
		File insBed = "temp/${sample}/bed/integrated/${hap}/svindel_ins.bed.gz"
		File delBed = "temp/${sample}/bed/integrated/${hap}/svindel_del.bed.gz"
		File snvBed = "temp/${sample}/bed/integrated/${hap}/snv_snv.bed.gz"
		File invBed = "temp/${sample}/bed/integrated/${hap}/svindel_inv.bed.gz"
	}
}

task call_integrate_sources_h2 {
	input {
		String hap
		File preInvSvBed
		File preInvSnvBed
		File delBed
		File insBed
		File invBed
		File invBatch
		String threads
		String mem_gb
	}
	command {
		snakemake --threads ${threads} temp/${sample}/bed/integrated/${hap}/svindel_ins.bed.gz temp/${sample}/bed/integrated/${hap}/svindel_del.bed.gz temp/${sample}/bed/integrated/${hap}/snv_snv.bed.gz temp/${sample}/bed/integrated/${hap}/svindel_inv.bed.gz
	}
	output {
		File insBed = "temp/${sample}/bed/integrated/${hap}/svindel_ins.bed.gz"
		File delBed = "temp/${sample}/bed/integrated/${hap}/svindel_del.bed.gz"
		File snvBed = "temp/${sample}/bed/integrated/${hap}/snv_snv.bed.gz"
		File invBed = "temp/${sample}/bed/integrated/${hap}/svindel_inv.bed.gz"
	}
}

task call_merge_haplotypes_chrom_svindel_del {
	input {
		String svtype
		File delBed_h1
		File delBed_h2
		File callable_h1
		File callable_h2
		String chrom
		String threads
		String mem_gb
	}
	command {
		snakemake --threads ${threads} temp/${sample}/bed/bychrom/${svtype}/${chrom}.bed.gz
	}
	output {
		File bed = "temp/${sample}/bed/bychrom/${svtype}/${chrom}.bed.gz"
	}
}

task call_merge_haplotypes_chrom_svindel_ins {
	input {
		String svtype
		File insBed_h1
		File insBed_h2
		File callable_h1
		File callable_h2
		String chrom
		String threads
		String mem_gb
	}
	command {
		snakemake --threads ${threads} temp/${sample}/bed/bychrom/${svtype}/${chrom}.bed.gz
	}
	output {
		File bed = "temp/${sample}/bed/bychrom/${svtype}/${chrom}.bed.gz"
	}
}

task call_merge_haplotypes_chrom_svinv {
	input {
		String svtype
		File invBed_h1
		File invBed_h2
		File callable_h1
		File callable_h2
		String chrom
		String threads
		String mem_gb
	}
	command {
		snakemake --threads ${threads} temp/${sample}/bed/bychrom/${svtype}/${chrom}.bed.gz
	}
	output {
		File bed = "temp/${sample}/bed/bychrom/${svtype}/${chrom}.bed.gz"
	}
}

task call_merge_haplotypes_chrom_snv {
	input {
		String svtype
		File snvBed_h1
		File snvBed_h2
		File callable_h1
		File callable_h2
		String chrom
		String threads
		String mem_gb
	}
	command {
		snakemake --threads ${threads} temp/${sample}/bed/bychrom/${svtype}/${chrom}.bed.gz
	}
	output {
		File bed = "temp/${sample}/bed/bychrom/${svtype}/${chrom}.bed.gz"
	}
}

task call_merge_haplotypes_inv {
	input {
		String svtype
		File bed
		String threads
		String mem_gb
	}
	command {
		snakemake --threads ${threads} temp/${sample}/bed/merged/${svtype}.bed.gz
	}
	output {
		File bed = "temp/${sample}/bed/merged/${svtype}.bed.gz"
	}
}

task call_merge_haplotypes_svindel_del {
	input {
		String svtype
		File bed
		String threads
		String mem_gb
	}
	command {
		snakemake --threads ${threads} temp/${sample}/bed/merged/${svtype}.bed.gz
	}
	output {
		File bed = "temp/${sample}/bed/merged/${svtype}.bed.gz"
	}
}

task call_merge_haplotypes_svindel_ins {
	input {
		String svtype
		File bed
		String threads
		String mem_gb
	}
	command {
		snakemake --threads ${threads} temp/${sample}/bed/merged/${svtype}.bed.gz
	}
	output {
		File bed = "temp/${sample}/bed/merged/${svtype}.bed.gz"
	}
}

task call_merge_haplotypes_snv {
	input {
		String svtype
		File bed
		String threads
		String mem_gb
	}
	command {
		snakemake --threads ${threads} temp/${sample}/bed/merged/${svtype}.bed.gz
	}
	output {
		File bed = "temp/${sample}/bed/merged/${svtype}.bed.gz"
	}
}

task call_final_bed {
	input {
		File invBed
		File insBed
		File delBed
		File snvBed
		String threads
		String mem_gb
	}
	command {
		snakemake --threads ${threads} results/${sample}/bed/snv_snv.bed.gz results/${sample}/bed/indel_ins.bed.gz results/${sample}/bed/indel_del.bed.gz results/${sample}/bed/sv_ins.bed.gz results/${sample}/bed/sv_del.bed.gz results/${sample}/bed/sv_inv.bed.gz results/${sample}/bed/fa/indel_ins.fa.gz results/${sample}/bed/fa/indel_del.fa.gz results/${sample}/bed/fa/sv_ins.fa.gz results/${sample}/bed/fa/sv_del.fa.gz results/${sample}/bed/fa/sv_inv.fa.gz
	}
	output {
		File snvBed = "results/${sample}/bed/snv_snv.bed.gz"
		File indelInsBed = "results/${sample}/bed/indel_ins.bed.gz"
		File indelDelBed = "results/${sample}/bed/indel_del.bed.gz"
		File svInsBed = "results/${sample}/bed/sv_ins.bed.gz"
		File svDelBed = "results/${sample}/bed/sv_del.bed.gz"
		File inBed = "results/${sample}/bed/sv_inv.bed.gz"
		File indelInsFasta = "results/${sample}/bed/fa/indel_ins.fa.gz"
		File indelDelFasta = "results/${sample}/bed/fa/indel_del.fa.gz"
		File svInsFasta = "results/${sample}/bed/fa/sv_ins.fa.gz"
		File svDelFasta = "results/${sample}/bed/fa/sv_del.fa.gz"
		File invFasta "results/${sample}/bed/fa/sv_inv.fa.gz"
	}
}
