version 1.0


# Given a set of trio children sequenced at high coverage, and such that their
# read length distribution is bimodal, the workflow creates new readsets by
# sampling from the right mode at increasing coverages, and studies how coverage
# affects performance of several SV discovery tools.
#
workflow TriosCoverageEffect {
    input {
        String bucket_dir
        File children_ids
        Int bin_length
        Int max_read_length
        Array[Float] right_coverages
        Float min_right_coverage
        Float max_right_coverage
        Float coverage_quantum
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
        right_coverages: "Coverages of the right mode. Of each haplotype."
        max_right_coverage: "The max of <right_coverages>."
        coverage_quantum: "Of each haplotype. Assumed to be a float <=1.0."
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
                right_coverages = right_coverages,
                min_right_coverage = min_right_coverage,
                max_right_coverage = max_right_coverage,
                coverage_quantum = coverage_quantum,
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
# distribution of the union of all its flowcells is a Gaussian or a mixture of
# Gaussians. The procedure merges all the flowcells of the child; it creates
# several readsets by sampling from the rightmost Gaussian at several
# coverages; and it runs SV callers on every such readset.
#
# Remark: the task checkpoints at the file level, by storing files in
# <bucket_dir/child_id/coverage_effect> and by skipping their creation if they
# already exist.
#
task ProcessTrioChild {
    input {
        String child_id
        Int flowcells_size_gb
        Int bin_length
        Int max_read_length
        Array[Float] right_coverages
        Float min_right_coverage
        Float max_right_coverage
        Float coverage_quantum
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
        right_coverages: "Coverages of the right mode. Of each haplotype."
        max_right_coverage: "The max of <right_coverages>."
        coverage_quantum: "Of each haplotype. Assumed to be a float <=1.0."
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

        # Building the bin files needed for sampling reads from the distribution
        TEST=$(gsutil -q stat ~{bucket_dir}/~{child_id}/coverage_effect/bins/bin_0.bin && echo 0 || echo 1)
        if [ ${TEST} -eq 0 ]; then
            while : ; do
                TEST=$(gsutil -m cp "~{bucket_dir}/~{child_id}/coverage_effect/bins/*" . && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading files from <~{bucket_dir}/~{child_id}/coverage_effect/bins/>. Trying again..."
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
                    TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp ${FASTQ_FILE_LOCAL}.histogram ${FASTQ_FILE_LOCAL}.max ~{bucket_dir}/~{child_id}/coverage_effect/histograms/ && echo 0 || echo 1)
                    if [ ${TEST} -eq 1 ]; then
                        echo "Error uploading file <~{bucket_dir}/~{child_id}/coverage_effect/histograms/${FASTQ_FILE_LOCAL}.histogram>. Trying again..."
                        sleep ${GSUTIL_DELAY_S}
                    else
                        break
                    fi
                done
            done < ~{child_id}.fastqs
            java -cp ~{docker_dir} -Xmx~{ram_size_gb_effective}G BuildReadLengthBins list.txt ~{bin_length} ~{max_read_length} ${GENOME_LENGTH_HAPLOID} bin_
            head -n 40 bin_0.bin
            rm -f *.fq
            while : ; do
                TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp "bin_*" ~{bucket_dir}/~{child_id}/coverage_effect/bins/ && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error uploading bin files to <~{bucket_dir}/~{child_id}/coverage_effect/bins/>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        fi

        # Computing properties of the distribution
        N_MAXIMA=0; MEAN_LEFT=0; STD_LEFT=0; MEAN_RIGHT=0; STD_RIGHT=0
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
        
        # Creating readsets with no left coverage and multiple right coverages
        RIGHT_COVERAGES=~{sep='-' right_coverages}
        TEST=$(gsutil -q stat ~{bucket_dir}/~{child_id}/coverage_effect/long_coverage_~{child_id}/reads.fastq && echo 0 || echo 1)
        if [ ${TEST} -eq 0 ]; then
            while : ; do
                TEST=$(gsutil cp ~{bucket_dir}/~{child_id}/coverage_effect/long_coverage_~{child_id}/reads.fastq . && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading file <~{bucket_dir}/~{child_id}/coverage_effect/long_coverage_~{child_id}/reads.fastq>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        else            
            # Sampling from the right distribution only                
            java -cp ~{docker_dir}:~{docker_dir}/commons-math3.jar SampleReadsFromLengthBins ${MEAN_RIGHT} ${STD_RIGHT} ${MEAN_LEFT} ${STD_LEFT} 1 ~{bin_length} ~{max_read_length} bin_ ${GENOME_LENGTH_HAPLOID} ~{max_right_coverage} reads.fastq
            rm -f bin_*.bin
            java -Xmx~{ram_size_gb_effective}G -cp ~{docker_dir} Fastq2LengthHistogram reads.fastq ~{bin_length} ~{max_read_length} reads.fastq.histogram reads.fastq.max
            COVERAGE_EACH_HAPLOTYPE=$( sed -n '2~4p' reads.fastq | wc -c )
            echo "scale=8; ${COVERAGE_EACH_HAPLOTYPE} / (2.0*${GENOME_LENGTH_HAPLOID})" | bc > reads.fastq.coverage
            while : ; do
                TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp "reads.fastq*" "~{bucket_dir}/~{child_id}/coverage_effect/long_coverage_~{child_id}/" && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error uploading file <~{bucket_dir}/~{child_id}/coverage_effect/long_coverage_~{child_id}/reads.fastq>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        fi
        bash ~{docker_dir}/readLengthDistribution_impl.sh reads.fastq -1 ~{child_id} ~{min_right_coverage} ~{max_right_coverage} ${RIGHT_COVERAGES} ~{reference_fa} ~{reference_fai} ~{reference_tandem_repeats} ~{bucket_dir}/~{child_id}/coverage_effect ~{use_pbsv} ~{use_sniffles1} ~{use_sniffles2} ~{use_hifiasm} ~{use_pav} ~{use_paftools} ~{keep_assemblies} ~{work_dir} ~{docker_dir} ~{coverage_quantum} ~{bin_length} ~{max_read_length} ${GENOME_LENGTH_HAPLOID}
        rm -f reads.fastq
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
