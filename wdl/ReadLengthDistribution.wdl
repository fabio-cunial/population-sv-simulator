version 1.0


# Studies how read length distribution affects SV discovery. 
#
workflow ReadLengthDistribution {
    input {
        Int n_trios
        Array[String] trio_fasta_child
        Array[String] trio_fai_child
        Array[String] trio_vcfs_child
        Array[String] trio_vcfs_parent1
        Array[String] trio_vcfs_parent2
        File reference_fa
        File reference_fai
        File reference_mmi
        File reference_tandem_repeats
        Array[Int] coverages
        String bucket_dir
        Int n_cpus_single_sample
        Int n_cpus_joint
        Int use_pbsv
        Int use_sniffles1
        Int use_sniffles2
        Int use_hifiasm
        Int use_pav
        Int use_paftools
        Int delete_bucket_dir
        Int keep_assemblies
    }
    parameter_meta {
        reference_mmi: "The minimap2 index of <reference_fa>."
        n_haplotypes: "Total number of trios in the arrays."
        coverages: "Of one diploid individual, sorted in increasing order."
        bucket_dir: "URL of a directory in a bucket. Stores both temporary files and output files."
        use_pbsv: "1=yes, 0=no."
        delete_bucket_dir: "Erases <bucket_dir> before running the analysis (1=yes, 0=no)."
        keep_assemblies: "Stores two assemblies per individual and per configuration in <bucket_dir>. 1=yes, 0=no."
    }
    if (delete_bucket_dir == 1) {
        call DeleteBucketDir {
            input:
                bucket_dir = bucket_dir
        }
    }
    scatter(i range(n_trios)) {
        call ProcessChunk {
            input:
                 id_from = i,
                 chunk_size = GetChunks.chunk_size,
                 reference_fa = reference_fa,
                 reference_fai = reference_fai,
                 reference_mmi = reference_mmi,
                 reference_tandem_repeats = reference_tandem_repeats,
                 haplotype2variants_file = haplotype2variants_file,
                 variants_file = variants_file,
                 coverages = coverages,
                 length_min = length_min,
                 length_max = length_max,
                 length_stdev = length_stdev,
                 lengths = length_means,
                 bucket_dir = bucket_dir,
                 n_cpus = n_cpus_single_sample,
                 use_pbsv = use_pbsv,
                 use_sniffles1 = use_sniffles1,
                 use_sniffles2 = use_sniffles2,
                 use_hifiasm = use_hifiasm,
                 use_pav = use_pav,
                 use_paftools = use_paftools,
                 keep_assemblies = keep_assemblies
        }
    }
    call JointCalling.JointCalling {
        input:
            bucket_dir = bucket_dir,
            coverages = coverages,
            lengths = length_means,
            n_individuals = n_haplotypes/2,
            reference_fa = reference_fa,
            use_pbsv = use_pbsv,
            use_sniffles2 = use_sniffles2,
            force_sequentiality = ProcessChunk.force_sequentiality
    }
    output {
    }
}


task DeleteBucketDir {
    input { 
        String bucket_dir
    }    
    command <<<
        set -euxo pipefail
        GSUTIL_DELAY_S="60"
        
        TEST=$(gsutil -q stat ~{bucket_dir} && echo 0 || echo 1)
        if [ ${TEST} -eq 0 ]; then
            while ; do
                TEST=$(gsutil -m rm -rf ~{bucket_dir} && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]
                    echo "Error deleting <~{bucket_dir}>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        fi
        echo "1" > force_sequentiality.txt
    >>>
    output {
        Int force_sequentiality = read_int("force_sequentiality.txt")
    }
    runtime {
        docker: "fcunial/simulation"
    }
}


# Given the VCFs created by a given SV caller on a trio at the maximum
# possible coverage, the task returns just the child calls that occur also in
# some parent.
#
# Remark: genotype information is ignored.
#
task CreateTruthVCF {
    input {
        Int child_id
        String vcf_child
        String vcf_parent1
        String vcf_parent2
        File reference_fa
        String bucket_dir
    }
    String docker_dir = "/simulation"
    String work_dir = "/cromwell_root/trios"
    Int ram_size_gb = ceil(size(vcf_child, "GB") + size(vcf_parent1, "GB") + size(vcf_parent2, "GB") + size(reference_fa, "GB"))*2
    Int disk_size_gb = ram_size_gb*2
    command <<<
        set -euxo pipefail
        cd ~{work_dir}
        TRUVARI_BENCH_FLAGS=" "  # Default settings for now
        BCFTOOLS_MERGE_FLAGS="--force-samples --merge none"
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        echo "Running <CreateTruthVCF> on ${N_THREADS} cores on the following node:"
        lscpu
        cat /proc/meminfo
        df -h
        
        
        # Creating the truth VCF
        TEST=$(gsutil -q stat ~{bucket_dir}/child~{child_id}/truth.vcf.gz && echo 0 || echo 1)
        if [ ${TEST} -eq 0 ]; then
            while : ; do
                TEST=$(gsutil cp ~{bucket_dir}/child~{child_id}/truth.vcf.gz . && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading file <~{bucket_dir}/child~{child_id}/truth.vcf.gz>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        else
            bcftools sort --output-type z --output child.vcf.gz ~{vcf_child}
            tabix child.vcf.gz
            bcftools sort --output-type z --output parent1.vcf.gz ~{vcf_parent1}
            tabix parent1.vcf.gz
            bcftools sort --output-type z --output parent2.vcf.gz ~{vcf_parent2}
            tabix parent2.vcf.gz
            ${TIME_COMMAND} truvari bench ${TRUVARI_BENCH_FLAGS} --prog --base parent1.vcf.gz --comp child.vcf.gz --reference ~{reference_fa} --output output_parent1/
            mv output_parent1/tp-call.vcf in_parent1.vcf
            rm -rf output_parent1/
            ${TIME_COMMAND} truvari bench ${TRUVARI_BENCH_FLAGS} --prog --base parent2.vcf.gz --comp child.vcf.gz --reference ~{reference_fa} --output output_parent2/
            mv output_parent2/tp-call.vcf in_parent2.vcf
            rm -rf output_parent2/
            bcftools sort --output-type z --output in_parent1.vcf.gz in_parent1.vcf
            tabix in_parent1.vcf.gz
            bcftools sort --output-type z --output in_parent2.vcf.gz in_parent2.vcf
            tabix in_parent2.vcf
            ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} ${BCFTOOLS_MERGE_FLAGS} in_parent1.vcf.gz in_parent2.vcf.gz --output-type z --output truth.vcf.gz
            while : ; do
                TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp truth.vcf.gz ~{bucket_dir}/child~{child_id}/truth.vcf.gz && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error uploading file <~{bucket_dir}/child~{child_id}/truth.vcf.gz>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        fi
        
        
        
        
        
        
    >>>
    output {
        File vcf_child_truth = work_dir + "/out.vcf.gz"
    }
    runtime {
        docker: "fcunial/simulation"
        cpu: 4  # Arbitrary
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 3
    }
}
























# ----------->Given a chunk of simulated haplotype IDs, the task uses every pair of
# consecutive haplotype IDs in the chunk to build a diploid individual and to
# call variants on its read alignments, where reads are simulated at several
# settings of coverage and mean read length. Every individual is processed
# sequentially.
#
# Remark: the task checkpoints by storing files in <bucket_dir>: (1) a
# checkpoint file for each chunk, that contains lines <individual length
# coverage>, is used to skip through the three nested loops of the simulation;
# (2) the files generated by each step are used to skip that step if they
# already exist.
#
# File size analysis for 20x coverage of one haplotype:
# FILE     | PROGRAM   | SIZE 
# svsig.gz | pbsv      | 15 MB
# .snf     | sniffles2 | 2.5 MB
# .vcf     | any       | 8 MB 
#
task ProcessChild {
    input {
        Int id_from
        Int chunk_size
        File reference_fa
        File reference_fai
        File reference_mmi
        File reference_tandem_repeats
        File haplotype2variants_file
        File variants_file
        Array[Int] coverages
        Int length_min
        Int length_max
        Int length_stdev
        Array[Int] lengths
        String bucket_dir
        Int n_cpus
        Int use_pbsv
        Int use_sniffles1
        Int use_sniffles2
        Int use_hifiasm
        Int use_pav
        Int use_paftools
        Int keep_assemblies
    }
    parameter_meta {
        id_from: "First simulated haplotype to process (zero-based, inclusive)."
        reference_fa: "Contains just chr1, and the header is '>chr1'."
        reference_mmi: "The minimap2 index of <reference_fa>."
        reference_tandem_repeats: "Contains just chr1."
        haplotype2variants_file: "Created by the population simulator."
        variants_file: "Created by the population simulator."
        coverages: "Of one haplotype, sorted in increasing order."
        lengths: "Mean read length values, in no particular order."
        n_cpus: "Used only for setting the runtime."
        use_pbsv: "1=yes, 0=no."
        keep_assemblies: "Stores two assemblies per individual and per configuration in the remote bucket. 1=yes, 0=no."
    }
    
    Int min_coverage = coverages[0]
    Int last = length(coverages) - 1
    Int max_coverage = coverages[last]
    Int ram_size_gb = 4 + ceil((size(reference_fa, "GB")*max_coverage*2) * 3)
    # *3: from loosely rounding up PAV; +4: to leave some free RAM to the OS.
    Int disk_size_gb = ram_size_gb*2 + ceil( size(reference_fa, "GB")*3 + size(reference_mmi, "GB") + size(reference_tandem_repeats, "GB") + size(haplotype2variants_file, "GB") + size(variants_file, "GB") )
    String docker_dir = "/simulation"
    String work_dir = "/cromwell_root/simulation"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        ID_TO=$(( ~{id_from} + ~{chunk_size} - 1 ))
        GSUTIL_DELAY_S="600"
        
        
        
        
        # Initializing files for sampling reads
        while read FASTQ_FILE; do
            java -cp ~{docker_dir} Fastq2LengthHistogram ${FASTQ_FILE} ~{bin_length} ~{max_read_length} ${FASTQ_FILE}.histogram ${FASTQ_FILE}.max
            echo "Local maxima:"
            cat ${FASTQ_FILE}.max
            echo "Histogram:"
            cat ${FASTQ_FILE}.histogram
            while : ; do
                TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp ${FASTQ_FILE}.histogram ${FASTQ_FILE}.max ~{bucket_dir}/histograms/ && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error uploading file <~{bucket_dir}/histograms/${FASTQ_FILE}.histogram>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        done < ~{flowcells_file}
        java -cp ~{docker_dir} BuildReadLengthBuckets ~{flowcells_file} ~{bin_length} ~{max_read_length} bucket_
        
        
        
        
        
        
        
        
        CHECKPOINT_FILE="checkpoint_i~{id_from}_i${ID_TO}.txt"
        TEST=$(gsutil -q stat ~{bucket_dir}/checkpoints/${CHECKPOINT_FILE} && echo 0 || echo 1)
        if [ ${TEST} -eq 0 ]; then
            while : ; do
                TEST=$(gsutil cp ~{bucket_dir}/checkpoints/${CHECKPOINT_FILE} . && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading file <~{bucket_dir}/checkpoints/${CHECKPOINT_FILE}>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        else
            echo "0 0 0" > ${CHECKPOINT_FILE}
        fi
        COVERAGES=~{sep='-' coverages}
        LENGTHS=~{sep='-' lengths}
        LENGTHS=$(echo ${LENGTHS} | tr '-' ' ')
        for ID1 in $(seq ~{id_from} 2 ${ID_TO}); do
            CHECKPOINT_INDIVIDUAL=$(tail -n1 ${CHECKPOINT_FILE} | awk '{ print $1 }')
            if [ ${ID1} -lt ${CHECKPOINT_INDIVIDUAL} ]; then
                continue
            fi
            ID2=$(( ${ID1} + 1 ))
            java -cp ~{docker_dir} -Xmx10g PrintHaplotypes ${ID1} ${ID2} ~{reference_fa} ~{haplotype2variants_file} ~{variants_file} .
            for LENGTH in ${LENGTHS}; do
                CHECKPOINT_LENGTH=$(tail -n1 ${CHECKPOINT_FILE} | awk '{ print $2 }')
                if [ ${ID1} -eq ${CHECKPOINT_INDIVIDUAL} -a ${LENGTH} -lt ${CHECKPOINT_LENGTH} ]; then
                    continue
                fi
                bash ~{docker_dir}/haplotype2reads.sh ${ID1} ${ID2} ~{length_min} ~{length_max} ${LENGTH} ~{length_stdev} ~{max_coverage} ~{bucket_dir} ~{work_dir}
                READS_FILE="reads_i${ID1}_i${ID2}_l${LENGTH}_c~{max_coverage}.fa"
                bash ~{docker_dir}/reads2svs.sh ${READS_FILE} ${ID1} ${LENGTH} ~{min_coverage} ~{max_coverage} ${COVERAGES} ~{reference_fa} ~{reference_fai} ~{reference_mmi} ~{reference_tandem_repeats} ${CHECKPOINT_FILE} ~{bucket_dir} ~{use_pbsv} ~{use_sniffles1} ~{use_sniffles2} ~{use_hifiasm} ~{use_pav} ~{use_paftools} ~{work_dir} ~{docker_dir} ~{keep_assemblies}
            done
            rm -f haplotype_${ID1}.fa haplotype_${ID2}.fa
        done
        echo "1" > force_sequentiality.txt
    >>>
    
    output {
        Int force_sequentiality = read_int(work_dir + "/force_sequentiality.txt")
    }
    runtime {
        docker: "fcunial/simulation"
        cpu: n_cpus
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 3
    }
}
