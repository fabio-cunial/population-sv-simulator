version 1.0


# Runs several callers on a child of a trio and on its parents, at max 
# coverage, and keeps only the variants of a child that occur also in some
# parent.
#
workflow TrioChildCreateTruthVCFs {
    input {
        String child_id
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
    
    Int flowcells_size_gb = 30*3  # Assuming 30x coverage per individual
    Int vcf_size_gb = 1  # Arbitrary upper bound
    call Child2Family {
        input: 
            child_id = child_id,
            bucket_dir = bucket_dir
    }
    scatter(individual_id in Child2Family.individuals) {
        call CreateFullCoverageVCFs {
            input:
                child_id = child_id,
                individual_id = individual_id,
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
    call IntersectVCFs {
        input:
            child_id = child_id,
            bucket_dir = bucket_dir,
            reference_fa = reference_fa,
            vcf_size_gb = vcf_size_gb,
            force_sequentiality = CreateFullCoverageVCFs.force_sequentiality
    }
    
    output {
    }
}


# Returns a WDL array with the IDs of a child and of its parents.
#
task Child2Family {
    input {
        String child_id
        String bucket_dir
    }
    command <<<
        while : ; do
            TEST=$(gsutil cp ~{bucket_dir}/trios_info/~{child_id}.parents . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading file <~{bucket_dir}/trios_info/~{child_id}.parents>. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        echo ~{child_id} >> ~{child_id}.parents
    >>>
    output {
        Array[String] individuals = read_lines(child_id + ".parents")
    }
    runtime {
        docker: "ubuntu:latest"
    }
}


# Given a trio child or parent, the task creates a VCF for each caller using
# all the available reads.
#
task CreateFullCoverageVCFs {
    input {
        String child_id
        String individual_id
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
    Int disk_size_gb = ram_size_gb*2 + ceil( size(reference_fa, "GB")*5 + size(reference_tandem_repeats, "GB") )
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
        
        echo 1 > force_sequentiality.txt
        TEST=$(gsutil cp ~{bucket_dir}/~{child_id}/full_coverage_~{individual_id}/reads.fastq . && echo 0 || echo 1)
        if [ ${TEST} -eq 0 ]; then
            while : ; do
                TEST=$(gsutil cp ~{bucket_dir}/~{child_id}/full_coverage_~{individual_id}/reads.fastq . && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading file <~{bucket_dir}/~{child_id}/full_coverage_~{individual_id}/reads.fastq>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
            while : ; do
                TEST=$(gsutil cp ~{bucket_dir}/~{child_id}/full_coverage_~{individual_id}/reads.bam . && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading file <~{bucket_dir}/~{child_id}/full_coverage_~{individual_id}/reads.bam>. Trying again..."
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
            rm -f reads.fastq
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
                cat ${FASTQ_FILE_LOCAL} >> reads.fastq
                rm -f ${FASTQ_FILE_LOCAL}
            done < ${INDIVIDUAL}.fastqs
            ${TIME_COMMAND} ${MINIMAP_COMMAND} -R ${READ_GROUP} ~{reference_fa} reads.fastq > reads.sam
            ${TIME_COMMAND} samtools sort -@ ${N_THREADS} --output-fmt BAM reads.sam > reads.1.bam
            rm -f reads.sam
            ${TIME_COMMAND} samtools calmd -@ ${N_THREADS} -b reads.1.sam ~{reference_fa} > reads.bam
            rm -f reads.1.bam
            while : ; do
                TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp reads.bam reads.fastq ~{bucket_dir}/~{child_id}/full_coverage_~{individual_id}/ && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error uploading file <~{bucket_dir}/~{child_id}/full_coverage_~{individual_id}/reads.bam>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        fi
        COVERAGE=$( sed -n '2~4p' reads.fastq | wc -c )
        COVERAGE=$(( ${COVERAGE} / (2*${GENOME_LENGTH_HAPLOID}) ))  # 1 haplotype
        bash ~{docker_dir}/reads2svs_impl.sh ~{individual_id} reads.bam reads.fastq ${COVERAGE} ~{individual_id} ${N_THREADS} ~{reference_fa} ~{reference_fai} ~{reference_tandem_repeats} ~{bucket_dir}/~{child_id}/full_coverage_~{individual_id} ~{use_pbsv} ~{use_sniffles1} ~{use_sniffles2} ~{use_hifiasm} ~{use_pav} ~{use_paftools} ~{keep_assemblies} ~{work_dir} ~{docker_dir}
        rm -f reads.bam reads.fastq
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


# Given the VCFs created by a given SV caller on a trio at max coverage, the
# task returns just the child calls that occur also in some parent.
#
task IntersectVCFs {
    input {
        String child_id
        String bucket_dir
        File reference_fa
        Int vcf_size_gb
        Array[Int] force_sequentiality
    }
    parameter_meta {
        vcf_size_gb: "Upper bound on the size of a VCF. Used just to set the runtime."
        force_sequentiality: "Used just for enforcing an order on tasks"
    }
    
    String docker_dir = "/simulation"
    String work_dir = "/cromwell_root/trios"
    Int ram_size_gb = vcf_size_gb*2 + size(reference_fa, "GB")*2
    Int disk_size_gb = vcf_size_gb*3 + size(reference_fa, "GB")
    
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

        TEST=$(gsutil -q stat ~{bucket_dir}/~{child_id}/truth.vcf.gz && echo 0 || echo 1)
        if [ ${TEST} -eq 1 ]; then
            while : ; do
                TEST=$(gsutil -m cp "~{bucket_dir}/~{child_id}/full_coverage_~{child_id}/*.vcf" . && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading files <~{bucket_dir}/~{child_id}/full_coverage_~{child_id}/*.vcf>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
            while : ; do
                TEST=$(gsutil cp ~{bucket_dir}/trios_info/~{child_id}.parents . && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading file <~{bucket_dir}/trios_info/~{child_id}.parents>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
            PARENT1=""; PARENT2=""; i=0
            while read PARENT_ID; do
                while : ; do
                    TEST=$(gsutil -m cp "~{bucket_dir}/~{child_id}/full_coverage_${PARENT_ID}/*.vcf" . && echo 0 || echo 1)
                    if [ ${TEST} -eq 1 ]; then
                        echo "Error downloading files <~{bucket_dir}/~{child_id}/full_coverage_${PARENT_ID}/*.vcf>. Trying again..."
                        sleep ${GSUTIL_DELAY_S}
                    else
                        break
                    fi
                done
                if [ $i -eq 0 ]; then
                    PARENT1=${PARENT_ID}
                else
                    PARENT2=${PARENT_ID}
                fi
                i=$(( $i + 1 ))
            done < ~{child_id}.parents
            for FILE in $(find . -type f -name "*_~{child_id}.vcf"); do
                FILE=$(basename ${FILE})
                CALLER=${FILE%_~{child_id}.vcf}
                bcftools sort --output-type z --output child.vcf.gz ${FILE}
                tabix child.vcf.gz
                bcftools sort --output-type z --output parent1.vcf.gz ${CALLER}_${PARENT1}.vcf
                tabix parent1.vcf.gz
                bcftools sort --output-type z --output parent2.vcf.gz ${CALLER}_${PARENT2}.vcf
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
                rm -f in_parent1.vcf.gz in_parent2.vcf.gz
                while : ; do
                    TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp truth.vcf.gz ~{bucket_dir}/~{child_id}/ && echo 0 || echo 1)
                    if [ ${TEST} -eq 1 ]; then
                        echo "Error uploading file <~{bucket_dir}/~{child_id}/truth.vcf.gz>. Trying again..."
                        sleep ${GSUTIL_DELAY_S}
                    else
                        break
                    fi
                done
            done
        fi
    >>>
    output {
    }
    runtime {
        docker: "fcunial/simulation"
        cpu: 8  # Arbitrary
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 3
    }
}