version 1.0


# Studies how a bimodal read length distribution affects SV discovery, by
# sampling reads from high-coverage trio children.
#
workflow ReadLengthDistribution {
    input {
        String bucket_dir
        File children_ids
        Int bin_length
        Int max_read_length
        Int answer_question1
        Array[Float] question1_left_weights
        Float question1_target_coverage_one_haplotype
        Int answer_question2
        Array[Int] question2_left_coverages
        Int question2_min_coverage
        Int question2_max_coverage
        File reference_fa
        File reference_fai
        File reference_tandem_repeats
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
        children_ids: "The trio children to be processed. Every child X is assumed to be described by files <X.fastqs> and <X.parents> in <bucket_dir/trios_info>."
        bin_length: "Of the read length histogram"
        max_read_length: "Of the read length histogram"
        answer_question1: "Answers the question: Given a fixed coverage, how does performance vary with different bimodal read length distributions?"
        question1_left_weights: "Weights to assign to the left Gaussian of a bimodal distribution"
        question1_target_coverage_one_haplotype: "Desired coverage of a readset generated by sampling the bimodal"
        answer_question2: "Answers the question: Given a fixed coverage around the right mode, how does performance vary with increasing coverage on the left mode?"
        question2_left_coverages: "Coverages of the left mode. Of one haplotype."
        question2_max_coverage: "The max of <question2_left_coverages>. Used also as the fixed coverage of the right mode."
        n_cpus: "Used only for setting the runtime"
        use_pbsv: "1=yes, 0=no."
        keep_assemblies: "Stores two assemblies per individual and per configuration in the remote bucket. 1=yes, 0=no."
    }
    
    Int flowcells_size_gb = 30*3  # Assuming 30x coverage per individual
    scatter(child in read_lines(children_ids)) {
        call ProcessTrioChild {
            input:
                child_id = child,
                flowcells_size_gb = flowcells_size_gb,
                bin_length = bin_length,
                max_read_length = max_read_length,
                answer_question1 = answer_question1,
                question1_left_weights = question1_left_weights,
                question1_target_coverage_one_haplotype = question1_target_coverage_one_haplotype,
                answer_question2 = answer_question2,
                question2_left_coverages = question2_left_coverages,
                question2_min_coverage = question2_min_coverage,
                question2_max_coverage = question2_max_coverage,
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                reference_tandem_repeats = reference_tandem_repeats,
                bucket_dir = bucket_dir,
                n_cpus = n_cpus,
                use_pbsv = use_pbsv,
                use_sniffles1 = use_sniffles1,
                use_sniffles2 = use_sniffles2,
                use_hifiasm = use_hifiasm,
                use_pav = use_pav,
                use_paftools = use_paftools,
                keep_assemblies = keep_assemblies
        }
    }
    output {
    }
}


# Assume that a trio child has high coverage, and that the read length
# distribution of the union of all its flowcells is the mixture of at least two
# Gaussians. The procedure merges all the flowcells of the child; it creates
# several bimodal distributions using the last two Gaussians, assigning a
# different weight to the last-but-one Gaussian; it creates a read set (at a
# small fixed coverage) by sampling from every such bimodal distribution; and
# it runs SV callers on every such read set.
#
# Remark: the task checkpoints at the file level, by storing files in
# <bucket_dir/child_id> and by skipping their creation if they already exist.
#
task ProcessTrioChild {
    input {
        String child_id
        Int flowcells_size_gb
        Int bin_length
        Int max_read_length
        Int answer_question1
        Array[Float] question1_left_weights
        Float question1_target_coverage_one_haplotype
        Int answer_question2
        Array[Int] question2_left_coverages
        Int question2_min_coverage
        Int question2_max_coverage
        File reference_fa
        File reference_fai
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
        flowcells_size_gb: "Upper bound on the size of the union of all flowcells of the child. Used just for setting the runtime."
        bin_length: "Of the read length histogram"
        max_read_length: "Of the read length histogram"
        answer_question1: "Answers the question: Given a fixed coverage, how does performance vary with different bimodal read length distributions?"
        question1_left_weights: "Weights to assign to the left Gaussian of a bimodal distribution"
        question1_target_coverage_one_haplotype: "Desired coverage of a readset generated by sampling the bimodal"
        answer_question2: "Answers the question: Given a fixed coverage around the right mode, how does performance vary with increasing coverage on the left mode?"
        question2_left_coverages: "Coverages of the left mode. Of one haplotype."
        question2_max_coverage: "The max of <question2_left_coverages>. Used also as the fixed coverage of the right mode."
        n_cpus: "Used only for setting the runtime"
        use_pbsv: "1=yes, 0=no."
        keep_assemblies: "Stores two assemblies per individual and per configuration in the remote bucket. 1=yes, 0=no."
    }
    
    Int ram_size_gb = 4 + flowcells_size_gb*3
    # *3: from loosely rounding up PAV; +4: to leave some free RAM to the OS.
    Int ram_size_gb_effective = ram_size_gb - 4
    Int disk_size_gb = ram_size_gb*2 + flowcells_size_gb*8 + ceil(size(reference_fa, "GB") + size(reference_tandem_repeats, "GB"))
    String docker_dir = "/simulation"
    String work_dir = "/cromwell_root/simulation"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        TIME_COMMAND="/usr/bin/time --verbose"
        READ_GROUP="@RG\tID:movie\tSM:~{child_id}"
        MINIMAP_COMMAND="minimap2 -t ${N_THREADS} -aYx map-hifi --eqx"
        GENOME_LENGTH_HAPLOID=$(cut -f 2 ~{reference_fai} | awk '{s+=$1} END {printf "%lu", s}')
        echo "Running <ProcessTrioChild> on ${N_THREADS} cores on the following node:"
        cat /proc/meminfo
        df -h

        # Building the files needed for sampling reads from the distribution
        TEST=$(gsutil -q stat ~{bucket_dir}/~{child_id}/bins/bin_0.bin && echo 0 || echo 1)
        if [ ${TEST} -eq 0 ]; then
            while : ; do
                TEST=$(gsutil -m cp "~{bucket_dir}/~{child_id}/bins/*" . && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading files from <~{bucket_dir}/~{child_id}/bins/>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        else
            while : ; do
                TEST=$(gsutil cp ~{bucket_dir}/trios_info/~{child_id}.fastqs . && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading file <~{bucket_dir}/trios_info/~{child_id}.fastqs>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
            rm -f list.txt
            while read FASTQ_FILE; do
                while : ; do
                    TEST=$(gsutil cp ${FASTQ_FILE} . && echo 0 || echo 1)
                    if [ ${TEST} -eq 1 ]; then
                        echo "Error downloading file <${FASTQ_FILE}>. Trying again..."
                        sleep ${GSUTIL_DELAY_S}
                    else
                        break
                    fi
                done
                FASTQ_FILE_LOCAL=$(basename ${FASTQ_FILE})
                gunzip ${FASTQ_FILE_LOCAL}
                FASTQ_FILE_LOCAL=${FASTQ_FILE_LOCAL%.gz}
                echo ${FASTQ_FILE_LOCAL} >> list.txt
                java -cp ~{docker_dir} Fastq2LengthHistogram ${FASTQ_FILE_LOCAL} ~{bin_length} ~{max_read_length} ${FASTQ_FILE_LOCAL}.histogram ${FASTQ_FILE_LOCAL}.max
                echo "File: ${FASTQ_FILE_LOCAL}"
                echo "Local maxima:"
                cat ${FASTQ_FILE_LOCAL}.max
                echo "Histogram:"
                cat ${FASTQ_FILE_LOCAL}.histogram
                while : ; do
                    TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp ${FASTQ_FILE_LOCAL}.histogram ${FASTQ_FILE_LOCAL}.max ~{bucket_dir}/~{child_id}/histograms/ && echo 0 || echo 1)
                    if [ ${TEST} -eq 1 ]; then
                        echo "Error uploading file <~{bucket_dir}/~{child_id}/histograms/${FASTQ_FILE_LOCAL}.histogram>. Trying again..."
                        sleep ${GSUTIL_DELAY_S}
                    else
                        break
                    fi
                done
            done < ~{child_id}.fastqs
            java -cp ~{docker_dir} -Xmx~{ram_size_gb_effective}G BuildReadLengthBins list.txt ~{bin_length} ~{max_read_length} ${GENOME_LENGTH_HAPLOID} bin_
            rm -f *.fq
            while : ; do
                TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp "bin_*" ~{bucket_dir}/~{child_id}/bins/ && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error uploading bin files to <~{bucket_dir}/~{child_id}/bins/>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        fi

        # Aborting if the distribution does not have at least two local maxima
        N_MAXIMA=0
        MEAN_LEFT=0
        STD_LEFT=0
        MEAN_RIGHT=0
        STD_RIGHT=0
        i=0
        while read line; do
            if [ $i -eq 0 ]; then
                N_MAXIMA=${line}
                if [ ${N_MAXIMA} -lt 2 ]; then
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
        if [ ${N_MAXIMA} -lt 2 ]; then
            echo "The read length distribution of the union of all flowcells has ${N_MAXIMA}<2 local maxima. Aborting."
            echo "Statistics:"
            cat bin_.stats
            exit 0
        fi

        # Creating readsets with a fixed coverage but multiple length
        # distributions
        if [ ~{answer_question1} -eq 1 ]; then
            WEIGHTS=~{sep='-' question1_left_weights}
            WEIGHTS=$(echo ${WEIGHTS} | tr '-' ' ')
            for WEIGHT_LEFT in ${WEIGHTS}; do
                EXISTS="0"
                TEST1=$(gsutil -q stat ~{bucket_dir}/~{child_id}/reads_w${WEIGHT_LEFT}/reads.fastq && echo 0 || echo 1)
                TEST2=$(gsutil -q stat ~{bucket_dir}/~{child_id}/reads_w${WEIGHT_LEFT}/reads.bam && echo 0 || echo 1)
                if [ ${TEST1} -eq 0 -a ${TEST2} -eq 0 ]; then
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
                        TEST=$(gsutil cp "~{bucket_dir}/~{child_id}/reads_w${WEIGHT_LEFT}/reads.bam*" . && echo 0 || echo 1)
                        if [ ${TEST} -eq 1 ]; then
                            echo "Error downloading <~{bucket_dir}/~{child_id}/reads_w${WEIGHT_LEFT}/reads.bam>. Trying again..."
                            sleep ${GSUTIL_DELAY_S}
                        else
                            break
                        fi
                    done
                    if [ -s reads.fastq -a -s reads.bam ]; then
                        EXISTS="1"
                    fi
                fi
                if [ ${EXISTS} -eq 0 ]; then
                    TEST=$(java -Xms~{ram_size_gb_effective}G -Xmx~{ram_size_gb_effective}G -cp ~{docker_dir}:~{docker_dir}/commons-math3.jar SampleReadsFromLengthBins ${MEAN_LEFT} ${STD_LEFT} ${MEAN_RIGHT} ${STD_RIGHT} ${WEIGHT_LEFT} ~{bin_length} ~{max_read_length} bin_ ${GENOME_LENGTH_HAPLOID} ~{question1_target_coverage_one_haplotype} ~{flowcells_size_gb} reads.fastq && echo 0 || echo 1)
                    if [ ${TEST} -eq 1 ]; then
                        rm -f reads.fastq.histogram; touch reads.fastq.histogram 
                        rm -f reads.fastq.max; touch reads.fastq.max
                        rm -f reads.fastq; touch reads.fastq
                        rm -f reads.bam; touch reads.bam
                        rm -f reads.bam.bai; touch reads.bam.bai
                    fi
                    while : ; do
                        TEST2=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp "reads.fastq*" ~{bucket_dir}/~{child_id}/reads_w${WEIGHT_LEFT}/ && echo 0 || echo 1)
                        if [ ${TEST2} -eq 1 ]; then
                            echo "Error uploading file <~{bucket_dir}/~{child_id}/reads_w${WEIGHT_LEFT}/reads.fastq>. Trying again..."
                            sleep ${GSUTIL_DELAY_S}
                        else
                            break
                        fi
                    done
                    if [ ${TEST} -eq 0 ]; then
                        ${TIME_COMMAND} ${MINIMAP_COMMAND} -R ${READ_GROUP} ~{reference_fa} reads.fastq > reads.sam
                        ${TIME_COMMAND} samtools sort -@ ${N_THREADS} --output-fmt BAM reads.sam > reads.1.bam
                        rm -f reads.sam
                        ${TIME_COMMAND} samtools calmd -@ ${N_THREADS} -b reads.1.bam ~{reference_fa} > reads.bam
                        rm -f reads.1.bam
                        samtools index -@ ${N_THREADS} reads.bam
                    fi
                    while : ; do
                        TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp "reads.bam*" ~{bucket_dir}/~{child_id}/reads_w${WEIGHT_LEFT}/ && echo 0 || echo 1)
                        if [ ${TEST} -eq 1 ]; then
                            echo "Error uploading file <~{bucket_dir}/~{child_id}/reads_w${WEIGHT_LEFT}/reads.bam>. Trying again..."
                            sleep ${GSUTIL_DELAY_S}
                        else
                            break
                        fi
                    done
                fi
                if [ -s reads.fastq -a -s reads.bam ]; then
                    if [ ! -f reads.bam.bai ]; then
                        samtools index -@ ${N_THREADS} reads.bam
                    fi
                    COVERAGE=$( sed -n '2~4p' reads.fastq | wc -c )
                    COVERAGE=$(( ${COVERAGE} / (2*${GENOME_LENGTH_HAPLOID}) ))  # 1 haplotype
                    bash ~{docker_dir}/reads2svs_impl.sh ~{child_id} reads.bam reads.fastq ${COVERAGE} ~{child_id} ${N_THREADS} ~{reference_fa} ~{reference_fai} ~{reference_tandem_repeats} ~{bucket_dir}/~{child_id}/reads_w${WEIGHT_LEFT} ~{use_pbsv} ~{use_sniffles1} ~{use_sniffles2} ~{use_hifiasm} ~{use_pav} ~{use_paftools} ~{keep_assemblies} ~{work_dir} ~{docker_dir}
                else
                    echo "Empty FASTQ or BAM file for weight ${WEIGHT_LEFT}. Stopping."
                fi
                rm -f reads.fastq reads.bam reads.bai
            done
        fi
        
        # Creating readsets with a fixed right coverage but multiple left
        # coverages.
        if [ ~{answer_question2} -eq 1 ]; then
            EXISTS="0"
            TEST1=$(gsutil -q stat ~{bucket_dir}/~{child_id}/reads_maxCoverage_right.fastq && echo 0 || echo 1)
            TEST2=$(gsutil -q stat ~{bucket_dir}/~{child_id}/reads_maxCoverage_left.fastq && echo 0 || echo 1)
            if [ ${TEST1} -eq 0 -a ${TEST2} -eq 0 ]; then
                while : ; do
                    TEST=$(gsutil -m cp ~{bucket_dir}/~{child_id}/reads_maxCoverage_right.fastq . && echo 0 || echo 1)
                    if [ ${TEST} -eq 1 ]; then
                        echo "Error downloading <~{bucket_dir}/~{child_id}/reads_maxCoverage_right.fastq>. Trying again..."
                        sleep ${GSUTIL_DELAY_S}
                    else
                        break
                    fi
                done
                while : ; do
                    TEST=$(gsutil -m cp ~{bucket_dir}/~{child_id}/reads_maxCoverage_left.fastq . && echo 0 || echo 1)
                    if [ ${TEST} -eq 1 ]; then
                        echo "Error downloading <~{bucket_dir}/~{child_id}/reads_maxCoverage_left.fastq>. Trying again..."
                        sleep ${GSUTIL_DELAY_S}
                    else
                        break
                    fi
                done
                if [ -s reads_maxCoverage_right.fastq -a -s reads_maxCoverage_left.fastq ]; then
                    EXISTS="1"
                fi
            fi
            if [ ${EXISTS} -eq 0 ]; then
                TEST=$(java -Xms~{ram_size_gb_effective}G -Xmx~{ram_size_gb_effective}G -cp ~{docker_dir}:~{docker_dir}/commons-math3.jar SampleReadsFromLengthBins ${MEAN_LEFT} ${STD_LEFT} ${MEAN_RIGHT} ${STD_RIGHT} 0 ~{bin_length} ~{max_read_length} bin_ ${GENOME_LENGTH_HAPLOID} ~{question2_max_coverage} ~{flowcells_size_gb} reads_maxCoverage_right.fastq && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    rm -f reads_maxCoverage_right.fastq; touch reads_maxCoverage_right.fastq
                fi
                while : ; do
                    TEST2=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp "reads_maxCoverage_right.fastq" ~{bucket_dir}/~{child_id}/ && echo 0 || echo 1)
                    if [ ${TEST2} -eq 1 ]; then
                        echo "Error uploading file <~{bucket_dir}/~{child_id}/reads_maxCoverage_right.fastq>. Trying again..."
                        sleep ${GSUTIL_DELAY_S}
                    else
                        break
                    fi
                done
                TEST=$(java -Xms~{ram_size_gb_effective}G -Xmx~{ram_size_gb_effective}G -cp ~{docker_dir}:~{docker_dir}/commons-math3.jar SampleReadsFromLengthBins ${MEAN_LEFT} ${STD_LEFT} ${MEAN_RIGHT} ${STD_RIGHT} 1 ~{bin_length} ~{max_read_length} bin_ ${GENOME_LENGTH_HAPLOID} ~{question2_max_coverage} ~{flowcells_size_gb} reads_maxCoverage_left.fastq && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    rm -f reads_maxCoverage_left.fastq; touch reads_maxCoverage_left.fastq
                fi
                while : ; do
                    TEST2=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp "reads_maxCoverage_left.fastq" ~{bucket_dir}/~{child_id}/ && echo 0 || echo 1)
                    if [ ${TEST2} -eq 1 ]; then
                        echo "Error uploading file <~{bucket_dir}/~{child_id}/reads_maxCoverage_left.fastq>. Trying again..."
                        sleep ${GSUTIL_DELAY_S}
                    else
                        break
                    fi
                done
            fi
            if [ -s reads_maxCoverage_left.fastq -a -s reads_maxCoverage_right.fastq ]; then
                bash ~{docker_dir}/readLengthDistribution_impl.sh reads_maxCoverage_left.fastq reads_maxCoverage_right.fastq ~{child_id} ~{question2_min_coverage} ~{question2_max_coverage} ~{question2_left_coverages} ~{reference_fa} ~{reference_fai} ~{reference_tandem_repeats} ~{bucket_dir}/~{child_id} ~{use_pbsv} ~{use_sniffles1} ~{use_sniffles2} ~{use_hifiasm} ~{use_pav} ~{use_paftools} ~{keep_assemblies} ~{work_dir} ~{docker_dir}
            else
                echo "Empty left or right read sets for weight ${WEIGHT_LEFT}. Aborting."
            fi
            rm -f reads_maxCoverage_*.fastq
        fi
    >>>
    
    output {
    }
    runtime {
        docker: "fcunial/simulation"
        cpu: n_cpus
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        bootDiskSizeGb: 15
        preemptible: 0
    }
}
