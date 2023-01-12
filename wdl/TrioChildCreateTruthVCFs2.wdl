version 1.0

import "TriosTruthVCFUtils.wdl"


# Runs several callers on a child of a trio and on its parents, using, for each
# individual, a readset whose length histogram is identical to the full readset,
# but all the modes except the one centered on the rightmost peak are removed. 
# Finally, keeps only the variants of a child that occur also in some parent.
#
workflow TrioChildCreateTruthVCFs2 {
    input {
        String child_id
        String bucket_dir
        Int bin_length
        Int max_read_length
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
        n_cpus: "For the bottleneck part of creating VCFs from large read sets"
    }
    
    Int flowcells_size_gb = 30*3  # Assuming 30x coverage per individual
    Int vcf_size_gb = 1  # Arbitrary upper bound
    call TriosTruthVCFUtils.Child2Family {
        input: 
            child_id = child_id,
            bucket_dir = bucket_dir
    }
    scatter(individual_id in Child2Family.individuals) {
        call CreateLongReadsVCFs {
            input:
                child_id = child_id,
                individual_id = individual_id,
                bin_length = bin_length,
                max_read_length = max_read_length,
                flowcells_size_gb = flowcells_size_gb,
                bucket_dir = bucket_dir,
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                reference_tandem_repeats = reference_tandem_repeats,
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
    call TriosTruthVCFUtils.IntersectVCFs {
        input:
            child_id = child_id,
            bucket_dir = bucket_dir,
            prefix = "long_coverage",
            reference_fa = reference_fa,
            vcf_size_gb = vcf_size_gb,
            force_sequentiality = CreateLongReadsVCFs.force_sequentiality
    }
    
    output {
    }
}


# Given a trio child or parent, the task creates a VCF for each caller using
# a readset whose length histogram is identical to the full readset, but
# all the modes except the one centered on the rightmost peak are removed.
#
task CreateLongReadsVCFs {
    input {
        String child_id
        String individual_id
        Int bin_length
        Int max_read_length
        Int flowcells_size_gb
        String bucket_dir
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
        child_id: "The ID of a trio child"
        individual_id: "The ID of the child or of one of its parents"
        flowcells_size_gb: "Upper bound on the size of the union of all flowcells of the individual. Used just for setting the runtime."
    }
    
    Int ram_size_gb = 4 + flowcells_size_gb*3
    # *3: from loosely rounding up PAV; +4: to leave some free RAM to the OS.
    Int ram_size_gb_effective = ram_size_gb - 4
    Int disk_size_gb = ceil(ram_size_gb*4 + size(reference_fa, "GB")*5 + size(reference_tandem_repeats, "GB"))
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
        READ_GROUP="@RG\tID:movie\tSM:~{individual_id}"
        MINIMAP_COMMAND="minimap2 -t ${N_THREADS} -aYx map-hifi --eqx"
        GENOME_LENGTH_HAPLOID=$(cut -f 2 ~{reference_fai} | awk '{s+=$1} END {printf "%lu", s}')
        
        echo 1 > ~{work_dir}/force_sequentiality.txt
        TEST=$(gsutil cp ~{bucket_dir}/~{child_id}/long_coverage_~{individual_id}/reads.fastq . && echo 0 || echo 1)
        if [ ${TEST} -eq 0 ]; then
            while : ; do
                TEST=$(gsutil cp ~{bucket_dir}/~{child_id}/long_coverage_~{individual_id}/reads.fastq . && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading file <~{bucket_dir}/~{child_id}/long_coverage_~{individual_id}/reads.fastq>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        else
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
                    TEST=$(gsutil cp ~{bucket_dir}/trios_info/~{individual_id}.fastqs . && echo 0 || echo 1)
                    if [ ${TEST} -eq 1 ]; then
                        echo "Error downloading file <~{bucket_dir}/trios_info/~{individual_id}.fastqs>. Trying again..."
                        sleep ${GSUTIL_DELAY_S}
                    else
                        break
                    fi
                done
                rm -f reads.fastq list.txt
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
                done < ~{individual_id}.fastqs
                java -cp ~{docker_dir} BuildReadLengthBins list.txt ~{bin_length} ~{max_read_length} ${GENOME_LENGTH_HAPLOID} bin_
                rm -f *.fastq
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
                # Using all reads
                cat bin_* > reads.fastq
            else
                # Cleaning the distribution
                java -cp ~{docker_dir} DeleteLeftMode ${MEAN_RIGHT} bin_ ~{bin_length} ~{max_read_length} reads.fastq
            fi
            rm -f bin_*
            while : ; do
                TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp reads.fastq "~{bucket_dir}/~{child_id}/long_coverage_~{individual_id}/reads.fastq" && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error uploading file <~{bucket_dir}/~{child_id}/long_coverage_~{individual_id}/reads.fastq>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        fi
        TEST=$(gsutil cp ~{bucket_dir}/~{child_id}/long_coverage_~{individual_id}/reads.bam . && echo 0 || echo 1)
        if [ ${TEST} -eq 0 ]; then
            while : ; do
                TEST=$(gsutil cp "~{bucket_dir}/~{child_id}/long_coverage_~{individual_id}/reads.bam*" . && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading file <~{bucket_dir}/~{child_id}/long_coverage_~{individual_id}/reads.bam>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        else
            ${TIME_COMMAND} ${MINIMAP_COMMAND} -R ${READ_GROUP} ~{reference_fa} reads.fastq > reads.sam
            ${TIME_COMMAND} samtools sort -@ ${N_THREADS} --output-fmt BAM reads.sam > reads.1.bam
            rm -f reads.sam
            ${TIME_COMMAND} samtools calmd -@ ${N_THREADS} -b reads.1.bam ~{reference_fa} > reads.bam
            rm -f reads.1.bam
            samtools index -@ ${N_THREADS} reads.bam
            while : ; do
                TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp reads.bam reads.bam.bai ~{bucket_dir}/~{child_id}/long_coverage_~{individual_id}/ && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error uploading file <~{bucket_dir}/~{child_id}/long_coverage_~{individual_id}/reads.bam>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        fi
        if [ ! -f reads.bam.bai ]; then
            samtools index -@ ${N_THREADS} reads.bam
        fi
        COVERAGE=$( sed -n '2~4p' reads.fastq | wc -c )
        COVERAGE=$(echo "scale=8; ${COVERAGE} / (2.0*${GENOME_LENGTH_HAPLOID})" | bc)  # Of each haplotype
        bash ~{docker_dir}/reads2svs_impl.sh ~{individual_id} reads.bam reads.fastq ${COVERAGE} ~{individual_id} ${N_THREADS} ~{reference_fa} ~{reference_fai} ~{reference_tandem_repeats} ~{bucket_dir}/~{child_id}/long_coverage_~{individual_id} ~{use_pbsv} ~{use_sniffles1} ~{use_sniffles2} ~{use_hifiasm} ~{use_pav} ~{use_paftools} ~{keep_assemblies} ~{work_dir} ~{docker_dir}
        rm -f reads.bam reads.bai reads.fastq
    >>>
    
    output {
        Int force_sequentiality = read_int("simulation/force_sequentiality.txt")
    }
    runtime {
        docker: "fcunial/simulation"
        cpu: n_cpus
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        bootDiskSizeGb: 15
        preemptible: 3
    }
}
