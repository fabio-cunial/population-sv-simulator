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
        Array[Float] left_weights
        Float target_coverage_one_haplotype
        File reference_fa
        File reference_fai
        File reference_mmi
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
        left_weights: "Weights to be assigned to the left Gaussian of a bimodal distribution"
        target_coverage_one_haplotype: "Desired coverage of a readset generated by sampling the bimodal"
        n_cpus: "Used only for setting the runtime"
    }
    
    Int flowcells_size_gb = 30*3
    # Assuming 30x coverage per individual
    scatter(child in read_lines(children_ids)) {
        call ProcessTrioChild {
            input:
                child_id = child,
                flowcells_size_gb = flowcells_size_gb,
                bin_length = bin_length,
                max_read_length = max_read_length,
                left_weights = left_weights,
                target_coverage_one_haplotype = target_coverage_one_haplotype,
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                reference_mmi = reference_mmi,
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
# distribution of the union of all its flowcells is the mixture of two
# Gaussians. The procedure merges all the flowcells of the child, it creates
# several bimodal distributions with a different weight assigned to the left
# Gaussian, it creates a read set (at a smaller fixed coverage) by sampling
# from every such distribution, and it runs SV callers on every such read set.
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
        Array[Float] left_weights
        Float target_coverage_one_haplotype
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
        flowcells_size_gb: "Upper bound on the size of the union of all flowcells in the list. Used just for setting the runtime."
        bin_length: "Of the read length histogram"
        max_read_length: "Of the read length histogram"
        left_weights: "Weights to assign to the left Gaussian of a bimodal distribution"
        target_coverage_one_haplotype: "Desired coverage of a readset generated by sampling the bimodal"
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
        MINIMAP_COMMAND="minimap2 -t ${N_THREADS} -aYx map-hifi --eqx"

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
            java -cp ~{docker_dir} BuildReadLengthBins list.txt ~{bin_length} ~{max_read_length} bin_
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
                java -Xmx4G -cp ~{docker_dir} SampleReadsFromLengthBins ${MEAN_LEFT} ${STD_LEFT} ${MEAN_RIGHT} ${STD_RIGHT} ${WEIGHT_LEFT} ~{bin_length} ~{max_read_length} bin_ ${GENOME_LENGTH_HAPLOID} ~{target_coverage_one_haplotype} 2000000000 reads.fastq
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
            fi
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
