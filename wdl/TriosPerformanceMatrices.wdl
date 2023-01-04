version 1.0

import "AnnotateVCFs.wdl"


# The program compares the VCF file for each tuple (trio child, left weight,
# caller) to its corresponding "truth" file created by
# <TrioChildCreateTruthVCFs>. SVs of different lengths are analyzed separately.
# The results are collected in a matrix for every child and caller, whose lines
# follow the format <leftWeight, svLength, metric>.
#
workflow TriosPerformanceMatrices {
    input {
        File children_ids
        String bucket_dir
        Array[String] callers
        Array[Float] left_weights
        Array[Int] sv_lengths
        Int only_pass
        File reference_fa
    }
    parameter_meta {
        children_ids: "The trio children to be processed"
        only_pass: "(0/1) Use only variants with FILTER=PASS."
    }    

    scatter(child_id in read_lines(children_ids)) {
        call ChildPerformanceMatrices { 
            input:
                child_id = child_id,
                bucket_dir = bucket_dir,
                callers = callers,
                left_weights = left_weights,
                sv_lengths = sv_lengths,
                only_pass = only_pass,
                reference_fa = reference_fa
        }
    }
    output {
    }
}


task ChildPerformanceMatrices {
    input {
        String child_id
        String bucket_dir
        Array[String] callers
        Array[Float] left_weights
        Array[Int] sv_lengths
        Int only_pass
        File reference_fa
    }
    parameter_meta {
        only_pass: "Use only calls with FILTER=PASS (0/1)."
    }
    
    String docker_dir = "/simulation"
    String work_dir = "/cromwell_root/simulation"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        TRUVARI_BENCH_FLAGS=" "  # Using the default settings for now
        
        CALLERS=~{sep='-' callers}
        CALLERS=$(echo ${CALLERS} | tr '-' ' ')
        WEIGHTS=~{sep='-' left_weights}
        WEIGHTS=$(echo ${WEIGHTS} | tr '-' ' ')
        SV_LENGTHS=~{sep='-' sv_lengths}
        SV_LENGTHS=$(echo ${SV_LENGTHS} | tr '-' ' ')
        for caller in ${CALLERS}; do
            while : ; do
                TEST=$(gsutil -m cp "~{bucket_dir}/~{child_id}/${caller}_truth.vcf.gz" . && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading file <~{bucket_dir}/~{child_id}/${caller}_truth.vcf.gz>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
            FILTER_STRING_1="((SVLEN>0 && SVLEN>=${sv_length}) || (SVLEN<0 && SVLEN<=-${sv_length}))"
            if [ ~{only_pass} -eq 1 ]; then
                FILTER_STRING_1="(${FILTER_STRING_1}) && FILTER=\"PASS\""
            fi
            PREVIOUS_SV_LENGTH="0"
            for sv_length in ${SV_LENGTHS}; do
                FILTER_STRING_2="((SVLEN>0 && SVLEN>${PREVIOUS_SV_LENGTH} && SVLEN<=${sv_length}) || (SVLEN<0 && SVLEN<-${PREVIOUS_SV_LENGTH} && SVLEN>=-${sv_length}))"
                bcftools filter --threads 0 --include "${FILTER_STRING_1}" --output-type v ${caller}_truth.vcf.gz | bcftools sort --output-type z --output truth1_${sv_length}.vcf.gz
                tabix truth1_${sv_length}.vcf.gz
                bcftools filter --threads 0 --include "${FILTER_STRING_2}" --output-type v ${caller}_truth.vcf.gz | bcftools sort --output-type z --output truth2_${sv_length}.vcf.gz
                tabix truth2_${sv_length}.vcf.gz
                PREVIOUS_SV_LENGTH=${sv_length}
            done
            TP_MATRIX_1="${caller}_matrix1_tp.txt"
            FP_MATRIX_1="${caller}_matrix1_fp.txt"
            FN_MATRIX_1="${caller}_matrix1_fn.txt"
            PRECISION_MATRIX_1="${caller}_matrix1_precision.txt"
            RECALL_MATRIX_1="${caller}_matrix1_recall.txt"
            F1_MATRIX_1="${caller}_matrix1_f1.txt"
            TP_MATRIX_2="${caller}_matrix2_tp.txt"
            FP_MATRIX_2="${caller}_matrix2_fp.txt"
            FN_MATRIX_2="${caller}_matrix2_fn.txt"
            PRECISION_MATRIX_2="${caller}_matrix2_precision.txt"
            RECALL_MATRIX_2="${caller}_matrix2_recall.txt"
            F1_MATRIX_2="${caller}_matrix2_f1.txt"
            touch ${TP_MATRIX_1} ${FP_MATRIX_1} ${FN_MATRIX_1} ${PRECISION_MATRIX_1} ${RECALL_MATRIX_1} ${F1_MATRIX_1}
            touch ${TP_MATRIX_2} ${FP_MATRIX_2} ${FN_MATRIX_2} ${PRECISION_MATRIX_2} ${RECALL_MATRIX_2} ${F1_MATRIX_2}
            for weight in ${WEIGHTS}; do
                TEST=$(gsutil -q stat "~{bucket_dir}/~{child_id}/reads_${weight}/${caller}_~{child_id}.vcf" . && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    continue
                fi
                while : ; do
                    TEST=$(gsutil cp "~{bucket_dir}/~{child_id}/reads_${weight}/${caller}_~{child_id}.vcf" . && echo 0 || echo 1)
                    if [ ${TEST} -eq 1 ]; then
                        echo "Error downloading file <~{bucket_dir}/~{child_id}/reads_${weight}/${caller}_~{child_id}.vcf>. Trying again..."
                        sleep ${GSUTIL_DELAY_S}
                    else
                        break
                    fi
                done
                PREVIOUS_SV_LENGTH="0"
                for sv_length in ${SV_LENGTHS}; do
                    # Cumulative length
                    echo -n "${weight},${sv_length}," >> ${TP_MATRIX_1}
                    echo -n "${weight},${sv_length}," >> ${FP_MATRIX_1}
                    echo -n "${weight},${sv_length}," >> ${FN_MATRIX_1}
                    echo -n "${weight},${sv_length}," >> ${PRECISION_MATRIX_1}
                    echo -n "${weight},${sv_length}," >> ${RECALL_MATRIX_1}
                    echo -n "${weight},${sv_length}," >> ${F1_MATRIX_1}
                    bcftools filter --threads 0 --include "${FILTER_STRING_1}" --output-type v ${caller}_~{child_id}.vcf | bcftools sort --output-type z --output measured.vcf.gz
                    tabix measured.vcf.gz
                    truvari bench ${TRUVARI_BENCH_FLAGS} --prog --base truth1_${sv_length}.vcf.gz --comp measured.vcf.gz --reference ~{reference_fa} --output output_dir/
                    rm -f measured.vcf.gz
                    grep "\"TP-call\":" output_dir/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${TP_MATRIX_1}
                    grep "\"FP\":" output_dir/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${FP_MATRIX_1}
                    grep "\"FN\":" output_dir/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${FN_MATRIX_1}
                    grep "\"precision\":" output_dir/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${PRECISION_MATRIX_1}
                    grep "\"recall\":" output_dir/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${RECALL_MATRIX_1}
                    grep "\"f1\":" output_dir/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${F1_MATRIX_1}
                    rm -rf output_dir/
                    echo "" >> ${TP_MATRIX_1}; echo "" >> ${FP_MATRIX_1}; echo "" >> ${FN_MATRIX_1}; echo "" >> ${PRECISION_MATRIX_1}; echo "" >> ${RECALL_MATRIX_1}; echo "" >> ${F1_MATRIX_1}
                    # Length bins
                    FILTER_STRING_2="((SVLEN>0 && SVLEN>${PREVIOUS_SV_LENGTH} && SVLEN<=${sv_length}) || (SVLEN<0 && SVLEN<-${PREVIOUS_SV_LENGTH} && SVLEN>=-${sv_length}))"
                    if [ ~{only_pass} -eq 1 ]; then
                        FILTER_STRING_2="(${FILTER_STRING_2}) && FILTER=\"PASS\""
                    fi
                    echo -n "${weight},${sv_length}," >> ${TP_MATRIX_2}
                    echo -n "${weight},${sv_length}," >> ${FP_MATRIX_2}
                    echo -n "${weight},${sv_length}," >> ${FN_MATRIX_2}
                    echo -n "${weight},${sv_length}," >> ${PRECISION_MATRIX_2}
                    echo -n "${weight},${sv_length}," >> ${RECALL_MATRIX_2}
                    echo -n "${weight},${sv_length}," >> ${F1_MATRIX_2}
                    bcftools filter --threads 0 --include "${FILTER_STRING_2}" --output-type v ${caller}_~{child_id}.vcf | bcftools sort --output-type z --output measured.vcf.gz
                    tabix measured.vcf.gz
                    truvari bench ${TRUVARI_BENCH_FLAGS} --prog --base truth2_${sv_length}.vcf.gz --comp measured.vcf.gz --reference ~{reference_fa} --output output_dir/
                    rm -f measured.vcf.gz
                    grep "\"TP-call\":" output_dir/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${TP_MATRIX_2}
                    grep "\"FP\":" output_dir/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${FP_MATRIX_2}
                    grep "\"FN\":" output_dir/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${FN_MATRIX_2}
                    grep "\"precision\":" output_dir/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${PRECISION_MATRIX_2}
                    grep "\"recall\":" output_dir/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${RECALL_MATRIX_2}
                    grep "\"f1\":" output_dir/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${F1_MATRIX_2}
                    rm -rf output_dir/                    
                    echo "" >> ${TP_MATRIX_2}; echo "" >> ${FP_MATRIX_2}; echo "" >> ${FN_MATRIX_2}; echo "" >> ${PRECISION_MATRIX_2}; echo "" >> ${RECALL_MATRIX_2}; echo "" >> ${F1_MATRIX_2}
                    PREVIOUS_SV_LENGTH=${sv_length}
                done
            done
            while : ; do
                TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp ${TP_MATRIX_1} ${FP_MATRIX_1} ${FN_MATRIX_1} ${PRECISION_MATRIX_1} ${RECALL_MATRIX_1} ${F1_MATRIX_1} ~{bucket_dir}/~{child_id}/ && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error uploading cumulative length matrices. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
            while : ; do
                TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp ${TP_MATRIX_2} ${FP_MATRIX_2} ${FN_MATRIX_2} ${PRECISION_MATRIX_2} ${RECALL_MATRIX_2} ${F1_MATRIX_2} ~{bucket_dir}/~{child_id}/ && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error uploading binned length matrices. Trying again..."
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
        cpu: 1
        memory: "8GB"  # Arbitrary
        disks: "local-disk 50 HDD"  # Arbitrary
        preemptible: 0
    }
}
