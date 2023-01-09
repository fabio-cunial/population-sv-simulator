version 1.0


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
        docker: "fcunial/simulation"
        bootDiskSizeGb: 15
    }
}


# Given the VCFs created by a given SV caller on a trio, the task returns just
# the child calls that occur also in some parent.
#
task IntersectVCFs {
    input {
        String child_id
        String bucket_dir
        String prefix
        File reference_fa
        Int vcf_size_gb
        Array[Int]? force_sequentiality
    }
    parameter_meta {
        prefix: "VCFs are assumed to be in directories <bucket_dir>/<child_id>/<prefix>_<X>/, where <X> is the ID of the child or parent."
        vcf_size_gb: "Upper bound on the size of a VCF. Used just to set the runtime."
        force_sequentiality: "Used just for enforcing an order on tasks"
    }
    
    String docker_dir = "/simulation"
    String work_dir = "/cromwell_root/simulation"
    Int ram_size_gb = ceil(vcf_size_gb*2 + size(reference_fa, "GB")*2)
    Int disk_size_gb = ceil(vcf_size_gb*3 + size(reference_fa, "GB"))
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TRUVARI_BENCH_FLAGS=" "  # Default settings for now
        BCFTOOLS_MERGE_FLAGS="--force-samples --merge none"
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))

        while : ; do
            TEST=$(gsutil -m cp "~{bucket_dir}/~{child_id}/~{prefix}_~{child_id}/vcfs/*.vcf" . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading files <~{bucket_dir}/~{child_id}/~{prefix}_~{child_id}/vcfs/*.vcf>. Trying again..."
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
                TEST=$(gsutil -m cp "~{bucket_dir}/~{child_id}/~{prefix}_${PARENT_ID}/vcfs/*.vcf" . && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading files <~{bucket_dir}/~{child_id}/~{prefix}_${PARENT_ID}/vcfs/*.vcf>. Trying again..."
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
        
        # Fixing VCF format issues by some callers
        for FILE in $(find . -type f -name "*.vcf"); do
            bcftools view -h ${FILE} > header.txt
            tail -n 1 header.txt > columns.txt
            N_LINES=$(wc -l < header.txt)
            head -n $(( ${N_LINES} - 1 )) header.txt > newFile.vcf
            echo "##FILTER=<ID=STRANDBIAS,Description=\"STRANDBIAS\">" >> newFile.vcf
            cat columns.txt >> newFile.vcf
            bcftools view -H ${FILE} >> newFile.vcf
            rm -f ${FILE}
            mv newFile.vcf ${FILE}
            rm -f header.txt columns.txt
        done
        
        # Comparing VCFs
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
            tabix in_parent2.vcf.gz
            ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} ${BCFTOOLS_MERGE_FLAGS} in_parent1.vcf.gz in_parent2.vcf.gz --output-type z --output truth.vcf.gz
            rm -f in_parent1.vcf.gz in_parent2.vcf.gz
            while : ; do
                TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp truth.vcf.gz ~{bucket_dir}/~{child_id}/~{prefix}_${CALLER}_truth.vcf.gz && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error uploading file <~{bucket_dir}/~{child_id}/~{prefix}_${CALLER}_truth.vcf.gz>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        done
    >>>
    output {
    }
    runtime {
        docker: "fcunial/simulation"
        cpu: 8  # Arbitrary
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        bootDiskSizeGb: 15
        preemptible: 3
    }
}
