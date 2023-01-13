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


function updateMatrices() {
    local CALLER=$1
    local FILTER_STRING=$2
    local VALUE=$3
    local SV_LENGTH=$4
    local TP_MATRIX=$5
    local FP_MATRIX=$6
    local FN_MATRIX=$7
    local PRECISION_MATRIX=$8
    local RECALL_MATRIX=$9
    local F1_MATRIX=${10}
    local TRUTH_FILE=${11}
    
    echo -n "${VALUE},${SV_LENGTH}," >> ${TP_MATRIX}
    echo -n "${VALUE},${SV_LENGTH}," >> ${FP_MATRIX}
    echo -n "${VALUE},${SV_LENGTH}," >> ${FN_MATRIX}
    echo -n "${VALUE},${SV_LENGTH}," >> ${PRECISION_MATRIX}
    echo -n "${VALUE},${SV_LENGTH}," >> ${RECALL_MATRIX}
    echo -n "${VALUE},${SV_LENGTH}," >> ${F1_MATRIX}
    rm -f measured.vcf.gz measured.vcf.gz.tbi
    bcftools filter --threads 0 --include "${FILTER_STRING}" --output-type v ${CALLER}_${CHILD_ID}.vcf | bcftools sort --output-type z --output measured.vcf.gz
    tabix measured.vcf.gz
    truvari bench ${TRUVARI_BENCH_FLAGS} --prog --base ${TRUTH_FILE} --comp measured.vcf.gz --reference ${REFERENCE_FA} --output output_dir/
    rm -f measured.vcf.gz measured.vcf.gz.tbi
    grep "\"TP-call\":" output_dir/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${TP_MATRIX}
    grep "\"FP\":" output_dir/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${FP_MATRIX}
    grep "\"FN\":" output_dir/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${FN_MATRIX}
    grep "\"precision\":" output_dir/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${PRECISION_MATRIX}
    grep "\"recall\":" output_dir/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${RECALL_MATRIX}
    grep "\"f1\":" output_dir/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${F1_MATRIX}
    rm -rf output_dir/
    echo "" >> ${TP_MATRIX}; echo "" >> ${FP_MATRIX}; echo "" >> ${FN_MATRIX}; echo "" >> ${PRECISION_MATRIX}; echo "" >> ${RECALL_MATRIX}; echo "" >> ${F1_MATRIX}
}


function uploadMatrices() {
    local TAG=$1
    local TP_MATRIX=$2
    local FP_MATRIX=$3
    local FN_MATRIX=$4
    local PRECISION_MATRIX=$5
    local RECALL_MATRIX=$6
    local F1_MATRIX=$7
    
    while : ; do
        TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp ${TP_MATRIX} ${FP_MATRIX} ${FN_MATRIX} ${PRECISION_MATRIX} ${RECALL_MATRIX} ${F1_MATRIX} ${BUCKET_DIR}/${CHILD_ID}/ && echo 0 || echo 1)
        if [ ${TEST} -eq 1 ]; then
            echo "Error uploading cumulative length matrices (${TAG}). Trying again..."
            sleep ${GSUTIL_DELAY_S}
        else
            break
        fi
    done
}


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
    PREVIOUS_SV_LENGTH="0"
    for sv_length in ${SV_LENGTHS}; do
        FILTER_STRING_1="(SVLEN>0 && SVLEN>=${sv_length}) || (SVLEN<0 && SVLEN<=-${sv_length})"
        FILTER_STRING_2="(SVLEN>0 && SVLEN>${PREVIOUS_SV_LENGTH} && SVLEN<=${sv_length}) || (SVLEN<0 && SVLEN<-${PREVIOUS_SV_LENGTH} && SVLEN>=-${sv_length})"
        FILTER_STRING_3="(SVLEN>0 && SVLEN<=${sv_length}) || (SVLEN<0 && SVLEN>=-${sv_length})"
        if [ ${ONLY_PASS} -eq 1 ]; then
            FILTER_STRING_1="(${FILTER_STRING_1}) && FILTER=\"PASS\""
            FILTER_STRING_2="(${FILTER_STRING_2}) && FILTER=\"PASS\""
            FILTER_STRING_3="(${FILTER_STRING_3}) && FILTER=\"PASS\""
        fi
        bcftools filter --threads 0 --include "${FILTER_STRING_1}" --output-type v ${TRUTH_VCF_PREFIX}_${caller}_truth.vcf.gz | bcftools sort --output-type z --output truth1_${sv_length}.vcf.gz
        tabix truth1_${sv_length}.vcf.gz
        bcftools filter --threads 0 --include "${FILTER_STRING_2}" --output-type v ${TRUTH_VCF_PREFIX}_${caller}_truth.vcf.gz | bcftools sort --output-type z --output truth2_${sv_length}.vcf.gz
        tabix truth2_${sv_length}.vcf.gz
        bcftools filter --threads 0 --include "${FILTER_STRING_3}" --output-type v ${TRUTH_VCF_PREFIX}_${caller}_truth.vcf.gz | bcftools sort --output-type z --output truth3_${sv_length}.vcf.gz
        tabix truth3_${sv_length}.vcf.gz
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
    TP_MATRIX_3="${TRUTH_VCF_PREFIX}_${caller}_matrix3_tp.txt"
    FP_MATRIX_3="${TRUTH_VCF_PREFIX}_${caller}_matrix3_fp.txt"
    FN_MATRIX_3="${TRUTH_VCF_PREFIX}_${caller}_matrix3_fn.txt"
    PRECISION_MATRIX_3="${TRUTH_VCF_PREFIX}_${caller}_matrix3_precision.txt"
    RECALL_MATRIX_3="${TRUTH_VCF_PREFIX}_${caller}_matrix3_recall.txt"
    F1_MATRIX_3="${TRUTH_VCF_PREFIX}_${caller}_matrix3_f1.txt"
    touch ${TP_MATRIX_1} ${FP_MATRIX_1} ${FN_MATRIX_1} ${PRECISION_MATRIX_1} ${RECALL_MATRIX_1} ${F1_MATRIX_1}
    touch ${TP_MATRIX_2} ${FP_MATRIX_2} ${FN_MATRIX_2} ${PRECISION_MATRIX_2} ${RECALL_MATRIX_2} ${F1_MATRIX_2}
    touch ${TP_MATRIX_3} ${FP_MATRIX_3} ${FN_MATRIX_3} ${PRECISION_MATRIX_3} ${RECALL_MATRIX_3} ${F1_MATRIX_3}
    for value in ${VALUES}; do
        TEST=$(gsutil -q stat "${BUCKET_DIR}/${CHILD_ID}/reads_${MEASURED_CHARACTER_CODE}${value}/vcfs/${caller}_${CHILD_ID}.vcf" && echo 0 || echo 1)
        if [ ${TEST} -eq 1 ]; then
            continue
        fi
        while : ; do
            TEST=$(gsutil cp "${BUCKET_DIR}/${CHILD_ID}/reads_${MEASURED_CHARACTER_CODE}${value}/vcfs/${caller}_${CHILD_ID}.vcf" . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading file <${BUCKET_DIR}/${CHILD_ID}/reads_${MEASURED_CHARACTER_CODE}${value}/vcfs/${caller}_${CHILD_ID}.vcf>. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        PREVIOUS_SV_LENGTH="0"
        for sv_length in ${SV_LENGTHS}; do
            FILTER_STRING_1="(SVLEN>0 && SVLEN>=${sv_length}) || (SVLEN<0 && SVLEN<=-${sv_length})"
            FILTER_STRING_2="(SVLEN>0 && SVLEN>${PREVIOUS_SV_LENGTH} && SVLEN<=${sv_length}) || (SVLEN<0 && SVLEN<-${PREVIOUS_SV_LENGTH} && SVLEN>=-${sv_length})"
            FILTER_STRING_3="(SVLEN>0 && SVLEN<=${sv_length}) || (SVLEN<0 && SVLEN>=-${sv_length})"
            if [ ${ONLY_PASS} -eq 1 ]; then
                FILTER_STRING_1="(${FILTER_STRING_1}) && FILTER=\"PASS\""
                FILTER_STRING_2="(${FILTER_STRING_2}) && FILTER=\"PASS\""
                FILTER_STRING_3="(${FILTER_STRING_3}) && FILTER=\"PASS\""
            fi
            updateMatrices ${caller} "${FILTER_STRING_1}" ${value} ${sv_length} ${TP_MATRIX_1} ${FP_MATRIX_1} ${FN_MATRIX_1} ${PRECISION_MATRIX_1} ${RECALL_MATRIX_1} ${F1_MATRIX_1} truth1_${sv_length}.vcf.gz
            updateMatrices ${caller} "${FILTER_STRING_2}" ${value} ${sv_length} ${TP_MATRIX_2} ${FP_MATRIX_2} ${FN_MATRIX_2} ${PRECISION_MATRIX_2} ${RECALL_MATRIX_2} ${F1_MATRIX_2} truth2_${sv_length}.vcf.gz
            updateMatrices ${caller} "${FILTER_STRING_3}" ${value} ${sv_length} ${TP_MATRIX_3} ${FP_MATRIX_3} ${FN_MATRIX_3} ${PRECISION_MATRIX_3} ${RECALL_MATRIX_3} ${F1_MATRIX_3} truth3_${sv_length}.vcf.gz
            PREVIOUS_SV_LENGTH=${sv_length}
        done
    done
    uploadMatrices long ${TP_MATRIX_1} ${FP_MATRIX_1} ${FN_MATRIX_1} ${PRECISION_MATRIX_1} ${RECALL_MATRIX_1} ${F1_MATRIX_1}
    uploadMatrices binned ${TP_MATRIX_2} ${FP_MATRIX_2} ${FN_MATRIX_2} ${PRECISION_MATRIX_2} ${RECALL_MATRIX_2} ${F1_MATRIX_2}
    uploadMatrices short ${TP_MATRIX_3} ${FP_MATRIX_3} ${FN_MATRIX_3} ${PRECISION_MATRIX_3} ${RECALL_MATRIX_3} ${F1_MATRIX_3}
done