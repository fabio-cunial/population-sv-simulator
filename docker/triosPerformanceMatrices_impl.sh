#!/bin/bash
#
BUCKET_DIR=$1
CHILD_ID=$2
CALLERS=$3
VALUES=$4  # Separated by -
SV_LENGTHS=$5  # Separated by -
TRUTH_VCF_PREFIX=$6  # String before <_caller>
MEASURED_CHARACTER_CODE=$7  # Single character after <reads_>
ONLY_PASS=$8
REFERENCE_FA=$9
WORK_DIR=${10}

set -euxo pipefail
cd ${WORK_DIR}

GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
GSUTIL_DELAY_S="600"
TRUVARI_BENCH_FLAGS=" "  # Using the default settings for now
CALLERS=$(echo ${CALLERS} | tr '-' ' ')
VALUES=$(echo ${VALUES} | tr '-' ' ')
SV_LENGTHS=$(echo ${SV_LENGTHS} | tr '-' ' ')

for caller in ${CALLERS}; do
    while : ; do
        TEST=$(gsutil -m cp "${BUCKET_DIR}/${CHILD_ID}/${TRUTH_VCF_PREFIX}_${caller}_truth.vcf.gz" . && echo 0 || echo 1)
        if [ ${TEST} -eq 1 ]; then
            echo "Error downloading file <${BUCKET_DIR}/${CHILD_ID}/${TRUTH_VCF_PREFIX}_${caller}_truth.vcf.gz>. Trying again..."
            sleep ${GSUTIL_DELAY_S}
        else
            break
        fi
    done
    FILTER_STRING_1="((SVLEN>0 && SVLEN>=${sv_length}) || (SVLEN<0 && SVLEN<=-${sv_length}))"
    if [ ${ONLY_PASS} -eq 1 ]; then
        FILTER_STRING_1="(${FILTER_STRING_1}) && FILTER=\"PASS\""
    fi
    PREVIOUS_SV_LENGTH="0"
    for sv_length in ${SV_LENGTHS}; do
        FILTER_STRING_2="((SVLEN>0 && SVLEN>${PREVIOUS_SV_LENGTH} && SVLEN<=${sv_length}) || (SVLEN<0 && SVLEN<-${PREVIOUS_SV_LENGTH} && SVLEN>=-${sv_length}))"
        bcftools filter --threads 0 --include "${FILTER_STRING_1}" --output-type v ${TRUTH_VCF_PREFIX}_${caller}_truth.vcf.gz | bcftools sort --output-type z --output truth1_${sv_length}.vcf.gz
        tabix truth1_${sv_length}.vcf.gz
        bcftools filter --threads 0 --include "${FILTER_STRING_2}" --output-type v ${TRUTH_VCF_PREFIX}_${caller}_truth.vcf.gz | bcftools sort --output-type z --output truth2_${sv_length}.vcf.gz
        tabix truth2_${sv_length}.vcf.gz
        PREVIOUS_SV_LENGTH=${sv_length}
    done
    TP_MATRIX_1="${TRUTH_VCF_PREFIX}_${caller}_matrix1_tp.txt"
    FP_MATRIX_1="${TRUTH_VCF_PREFIX}_${caller}_matrix1_fp.txt"
    FN_MATRIX_1="${TRUTH_VCF_PREFIX}_${caller}_matrix1_fn.txt"
    PRECISION_MATRIX_1="${TRUTH_VCF_PREFIX}_${caller}_matrix1_precision.txt"
    RECALL_MATRIX_1="${TRUTH_VCF_PREFIX}_${caller}_matrix1_recall.txt"
    F1_MATRIX_1="${TRUTH_VCF_PREFIX}_${caller}_matrix1_f1.txt"
    TP_MATRIX_2="${TRUTH_VCF_PREFIX}_${caller}_matrix2_tp.txt"
    FP_MATRIX_2="${TRUTH_VCF_PREFIX}_${caller}_matrix2_fp.txt"
    FN_MATRIX_2="${TRUTH_VCF_PREFIX}_${caller}_matrix2_fn.txt"
    PRECISION_MATRIX_2="${TRUTH_VCF_PREFIX}_${caller}_matrix2_precision.txt"
    RECALL_MATRIX_2="${TRUTH_VCF_PREFIX}_${caller}_matrix2_recall.txt"
    F1_MATRIX_2="${TRUTH_VCF_PREFIX}_${caller}_matrix2_f1.txt"
    touch ${TP_MATRIX_1} ${FP_MATRIX_1} ${FN_MATRIX_1} ${PRECISION_MATRIX_1} ${RECALL_MATRIX_1} ${F1_MATRIX_1}
    touch ${TP_MATRIX_2} ${FP_MATRIX_2} ${FN_MATRIX_2} ${PRECISION_MATRIX_2} ${RECALL_MATRIX_2} ${F1_MATRIX_2}
    for value in ${VALUES}; do
        TEST=$(gsutil -q stat "${BUCKET_DIR}/${CHILD_ID}/reads_${MEASURED_CHARACTER_CODE}${value}/${caller}_${CHILD_ID}.vcf" . && echo 0 || echo 1)
        if [ ${TEST} -eq 1 ]; then
            continue
        fi
        while : ; do
            TEST=$(gsutil cp "${BUCKET_DIR}/${CHILD_ID}/reads_${MEASURED_CHARACTER_CODE}${value}/${caller}_${CHILD_ID}.vcf" . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading file <${BUCKET_DIR}/${CHILD_ID}/reads_${MEASURED_CHARACTER_CODE}${value}/${caller}_${CHILD_ID}.vcf>. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        PREVIOUS_SV_LENGTH="0"
        for sv_length in ${SV_LENGTHS}; do
            # Cumulative length
            echo -n "${value},${sv_length}," >> ${TP_MATRIX_1}
            echo -n "${value},${sv_length}," >> ${FP_MATRIX_1}
            echo -n "${value},${sv_length}," >> ${FN_MATRIX_1}
            echo -n "${value},${sv_length}," >> ${PRECISION_MATRIX_1}
            echo -n "${value},${sv_length}," >> ${RECALL_MATRIX_1}
            echo -n "${value},${sv_length}," >> ${F1_MATRIX_1}
            bcftools filter --threads 0 --include "${FILTER_STRING_1}" --output-type v ${caller}_${CHILD_ID}.vcf | bcftools sort --output-type z --output measured.vcf.gz
            tabix measured.vcf.gz
            truvari bench ${TRUVARI_BENCH_FLAGS} --prog --base truth1_${sv_length}.vcf.gz --comp measured.vcf.gz --reference ${REFERENCE_FA} --output output_dir/
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
            if [ ${ONLY_PASS} -eq 1 ]; then
                FILTER_STRING_2="(${FILTER_STRING_2}) && FILTER=\"PASS\""
            fi
            echo -n "${value},${sv_length}," >> ${TP_MATRIX_2}
            echo -n "${value},${sv_length}," >> ${FP_MATRIX_2}
            echo -n "${value},${sv_length}," >> ${FN_MATRIX_2}
            echo -n "${value},${sv_length}," >> ${PRECISION_MATRIX_2}
            echo -n "${value},${sv_length}," >> ${RECALL_MATRIX_2}
            echo -n "${value},${sv_length}," >> ${F1_MATRIX_2}
            bcftools filter --threads 0 --include "${FILTER_STRING_2}" --output-type v ${caller}_${CHILD_ID}.vcf | bcftools sort --output-type z --output measured.vcf.gz
            tabix measured.vcf.gz
            truvari bench ${TRUVARI_BENCH_FLAGS} --prog --base truth2_${sv_length}.vcf.gz --comp measured.vcf.gz --reference ${REFERENCE_FA} --output output_dir/
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
        TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp ${TP_MATRIX_1} ${FP_MATRIX_1} ${FN_MATRIX_1} ${PRECISION_MATRIX_1} ${RECALL_MATRIX_1} ${F1_MATRIX_1} ${BUCKET_DIR}/${CHILD_ID}/ && echo 0 || echo 1)
        if [ ${TEST} -eq 1 ]; then
            echo "Error uploading cumulative length matrices. Trying again..."
            sleep ${GSUTIL_DELAY_S}
        else
            break
        fi
    done
    while : ; do
        TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp ${TP_MATRIX_2} ${FP_MATRIX_2} ${FN_MATRIX_2} ${PRECISION_MATRIX_2} ${RECALL_MATRIX_2} ${F1_MATRIX_2} ${BUCKET_DIR}/${CHILD_ID}/ && echo 0 || echo 1)
        if [ ${TEST} -eq 1 ]; then
            echo "Error uploading binned length matrices. Trying again..."
            sleep ${GSUTIL_DELAY_S}
        else
            break
        fi
    done
done