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
























# Assume that a trio child has high coverage, and that the read length
# distribution of the union of all its flowcells is the mixture of two
# Gaussians. The procedure merges all the flowcells of the child, it creates
# several bimodal distributions with a different weight assigned to the left
# Gaussian, it creates a read set (at a smaller fixed coverage) by sampling
# from every such distribution, and it runs SV callers on every read set.
#
# Remark: the task checkpoints at the file level, by storing files in
# <bucket_dir> and by skipping their creation if they already exist.
#
task ProcessChild {
    input {
        Int child_id
        File flowcells_list
        Int flowcells_size_gb
        Int bin_length
        Int max_read_length
        Array[Float] left_weights
        Float target_coverage_haploid
        File reference_fa
        File reference_fai
        File reference_mmi
        File reference_tandem_repeats
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
        flowcells_list: "Contains the address of every flowcell of the child (one per line)."
        flowcells_size_gb: "Upper bound on the size of the union of all flowcells in the list. Used just for setting the runtime."
        bin_length: "Of the read length histogram"
        max_read_length: "Of the read length histogram"
        left_weights: "Weights to assign to the left Gaussian of a bimodal distribution"
        target_coverage_haploid: "Desired coverage of a readset generated by sampling the bimodal"
        n_cpus: "Used only for setting the runtime"
        use_pbsv: "1=yes, 0=no."
        keep_assemblies: "Stores two assemblies per individual and per configuration in the remote bucket. 1=yes, 0=no."
    }
    
    Int ram_size_gb = 4 + flowcells_size_gb*3
    # *3: from loosely rounding up PAV; +4: to leave some free RAM to the OS.
    Int disk_size_gb = ram_size_gb*2 + ceil( size(reference_fa, "GB")*3 + size(reference_mmi, "GB") + size(reference_tandem_repeats, "GB") )
    String docker_dir = "/simulation"
    String work_dir = "/cromwell_root/simulation"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        GSUTIL_DELAY_S="600"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        REFERENCE_LENGTH=$(cut -f 2 ~{reference_fai} | awk '{s+=$1} END {print s}')
        READ_GROUP="@RG\tID:movie\tSM:~{child_id}"
        
        # Building the files needed for sampling reads from the distribution
        TEST=$(gsutil -q stat ~{bucket_dir}/~{child_id}/bins/bin_0 && echo 0 || echo 1)
        if [ ${TEST} -eq 0 ]; then
            while : ; do
                TEST=$(gsutil cp "~{bucket_dir}/~{child_id}/bins/*" . && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading files from <~{bucket_dir}/~{child_id}/bins/>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        else
            while read FASTQ_FILE; do
                java -cp ~{docker_dir} Fastq2LengthHistogram ${FASTQ_FILE} ~{bin_length} ~{max_read_length} ${FASTQ_FILE}.histogram ${FASTQ_FILE}.max
                echo "File: ${FASTQ_FILE}"
                echo "Local maxima:"
                cat ${FASTQ_FILE}.max
                echo "Histogram:"
                cat ${FASTQ_FILE}.histogram
                while : ; do
                    TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp ${FASTQ_FILE}.histogram ${FASTQ_FILE}.max ~{bucket_dir}/~{child_id}/histograms/ && echo 0 || echo 1)
                    if [ ${TEST} -eq 1 ]; then
                        echo "Error uploading file <~{bucket_dir}/~{child_id}/histograms/${FASTQ_FILE}.histogram>. Trying again..."
                        sleep ${GSUTIL_DELAY_S}
                    else
                        break
                    fi
                done
            done < ~{flowcells_list}
            java -cp ~{docker_dir} BuildReadLengthBins ~{flowcells_list} ~{bin_length} ~{max_read_length} bin_
            while : ; do
                TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp "bin_*" ~{bucket_dir}/~{child_id}/bins/ && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error uploading bin files to <~{bucket_dir}/~{child_id}/bins/>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        done
        
        # Aborting if the distribution does not have exactly two local maxima
        N_MAXIMA=0
        MEAN_LEFT=0
        STD_LEFT=0
        MEAN_RIGHT=0
        STD_RIGHT=0
        i=0
        while read line; do
            if [ $i -eq 0 ]; then
                N_MAXIMA=${line}
                if [ ${N_MAXIMA} -ne 2 ]; then
                    break
                fi
            elif [ $i -eq 1 ]; then
                MEAN_LEFT=${line}
            elif [ $i -eq 2 ]; then
                STD_LEFT=${line}
            elif [ $i -eq 3 ]; then
                MEAN_RIGHT=${line}
            elif [ $i -eq 4 ]; then
                STD_RIGHT=${line}
            fi
            i=$(( $i + 1 ))
        done < bin_.stats
        if [ ${N_MAXIMA} -ne 2 ]; then
            echo "The read length distribution of the union of all flowcells has ${N_MAXIMA} local maxima. Aborting."
            echo "Statistics:"
            cat bin_.stats
            exit 0
        fi
        
        # Creating reads with multiple length distributions
        WEIGHTS=~{sep='-' left_weights}
        WEIGHTS=$(echo ${WEIGHTS} | tr '-' ' ')
        for WEIGHT_LEFT in ${WEIGHTS}; do
            TEST=$(gsutil -q stat ~{bucket_dir}/~{child_id}/reads_w${WEIGHT_LEFT}/reads.fastq && echo 0 || echo 1)
            if [ ${TEST} -eq 0 ]; then
                while : ; do
                    TEST=$(gsutil cp "~{bucket_dir}/~{child_id}/reads_w${WEIGHT_LEFT}/reads.fastq" . && echo 0 || echo 1)
                    if [ ${TEST} -eq 1 ]; then
                        echo "Error downloading <~{bucket_dir}/~{child_id}/reads_w${WEIGHT_LEFT}/reads.fastq>. Trying again..."
                        sleep ${GSUTIL_DELAY_S}
                    else
                        break
                    fi
                done
                while : ; do
                    TEST=$(gsutil cp "~{bucket_dir}/~{child_id}/reads_w${WEIGHT_LEFT}/reads.bam" . && echo 0 || echo 1)
                    if [ ${TEST} -eq 1 ]; then
                        echo "Error downloading <~{bucket_dir}/~{child_id}/reads_w${WEIGHT_LEFT}/reads.bam>. Trying again..."
                        sleep ${GSUTIL_DELAY_S}
                    else
                        break
                    fi
                done
            else
                GENOME_LENGTH_HAPLOID=$(cut -f 2 ~{reference_fai} | awk '{s+=$1} END {print s}')
                java -Xmx4G -cp ~{docker_dir} SampleReadsFromLengthBins ${MEAN_LEFT} ${STD_LEFT} ${MEAN_RIGHT} ${STD_RIGHT} ${WEIGHT_LEFT} ~{bin_length} ~{max_read_length} bin_ ${GENOME_LENGTH_HAPLOID} ~{target_coverage_haploid} 2000000000 reads.fastq
                if [ $? -ne 0 ]; then
                    rm -rf reads.fastq
                    touch reads.fastq
                    rm -rf reads.bam
                    touch reads.bam
                else 
                	${TIME_COMMAND} ${MINIMAP_COMMAND} -R ${READ_GROUP} ~{reference_mmi} reads.fastq > reads.sam
                	${TIME_COMMAND} samtools calmd -@ ${N_THREADS} -b reads.sam ~{reference_fa} > reads.1.bam
                	rm -f reads.sam
                	${TIME_COMMAND} samtools sort -@ ${N_THREADS} reads.1.bam > reads.bam
                	rm -f reads.1.bam
                fi
                while : ; do
                    TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp reads.fastq ~{bucket_dir}/~{child_id}/reads_w${WEIGHT_LEFT}/reads.fastq && echo 0 || echo 1)
                    if [ ${TEST} -eq 1 ]; then
                        echo "Error uploading file <~{bucket_dir}/~{child_id}/reads_w${WEIGHT_LEFT}/reads.fastq>. Trying again..."
                        sleep ${GSUTIL_DELAY_S}
                    else
                        break
                    fi
                fi
                while : ; do
                    TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp reads.bam ~{bucket_dir}/~{child_id}/reads_w${WEIGHT_LEFT}/reads.bam && echo 0 || echo 1)
                    if [ ${TEST} -eq 1 ]; then
                        echo "Error uploading file <~{bucket_dir}/~{child_id}/reads_w${WEIGHT_LEFT}/reads.bam>. Trying again..."
                        sleep ${GSUTIL_DELAY_S}
                    else
                        break
                    fi
                fi
            done
            COVERAGE=$( sed -n '2~4p' reads.fastq | wc -c )
            COVERAGE=$(( ${COVERAGE} / (2*${REFERENCE_LENGTH}) ))  # 1 haplotype
            bash ~{docker_dir}/reads2svs_impl.sh ~{child_id} reads.bam reads.fastq ${COVERAGE} ~{child_id} ${N_THREADS} ~{reference_fa} ~{reference_fai} ~{reference_tandem_repeats} ~{bucket_dir}/~{child_id}/reads_w${WEIGHT_LEFT} ~{use_pbsv} ~{use_sniffles1} ~{use_sniffles2} ~{use_hifiasm} ~{use_pav} ~{use_paftools} ~{keep_assemblies} ~{work_dir} ~{docker_dir}
            rm -f reads.fastq reads.bam
        done
    >>>
    
    output {
    }
    runtime {
        docker: "fcunial/simulation"
        cpu: n_cpus
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 3
    }
}
