#!/bin/bash
#
# Given a set of rules for filtering VCFs, the program: 
#
# 1. Compares every filtered VCF file in <measured_vcfs/> to the corresponding
# file in <ground_truth_vcfs/>, in parallel over all the cores of the machine,
# collecting per-individual statistics.
#
# 2. Merges all filtered VCF files in <measured_vcfs/>, and compares the result
# to the filtered set of distinct ground-truth variants in the population.
#
# 3. Compares a filtered joint-calling file to the filtered set of distinct
# ground-truth variants in the population.
#
FILTER_STRING=$1  # Spaces are replaced with '+'. It equals '+' if empty.
JOINT_CALLING_FILE=$2  # Relative path or "null".
REFERENCE_FA=$3
N_THREADS=$4
WORK_DIR=$5  # Absolute path
CONFIGURATION_ID=$6  # Used just for storing the all-individuals merged VCF
BUCKET_DIR_ALLINDIVIDUALS_VCFS=$7

# Per-individual matrices. The script appends many values (one value per
# individual) to a single row of each matrix.
TP_MATRIX=$8
FP_MATRIX=$9
FN_MATRIX=${10}
PRECISION_MATRIX=${11}
RECALL_MATRIX=${12}
F1_MATRIX=${13}

# All-individuals matrices. The script appends one aggregated value to a single
# row of each matrix.
TP_MATRIX_MERGE=${14}
FP_MATRIX_MERGE=${15}
FN_MATRIX_MERGE=${16}
PRECISION_MATRIX_MERGE=${17}
RECALL_MATRIX_MERGE=${18}
F1_MATRIX_MERGE=${19}

# Joint-calling matrices (ignored if there is no joint calling file). The
# script appends one aggregated value to a single row of each matrix.
TP_MATRIX_JOINT=${20}
FP_MATRIX_JOINT=${21}
FN_MATRIX_JOINT=${22}
PRECISION_MATRIX_JOINT=${23}
RECALL_MATRIX_JOINT=${24}
F1_MATRIX_JOINT=${25}


set -euxo pipefail
cd ${WORK_DIR}
GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
TIME_COMMAND="/usr/bin/time --verbose"
if [ ${FILTER_STRING} = "+" ]; then
    FILTER_STRING=""
else
    FILTER_STRING=$(echo ${FILTER_STRING} | tr '+' ' ')
fi


# Filters every VCF file in a chunk of files using $FILTER_STRING$, filters the
# corresponding ground truth as well, and compares the two, sequentially.
#
# Remark: ground truth VCFs must be indexed every time, since we are filtering
# them by different criteria every time.
#
function processChunk() {
    local CHUNK_ID=$1
    
    local TP_MATRIX="tp_${CHUNK_ID}.txt"
    local FP_MATRIX="fp_${CHUNK_ID}.txt"
    local FN_MATRIX="fn_${CHUNK_ID}.txt"
    local PRECISION_MATRIX="precision_${CHUNK_ID}.txt"
    local RECALL_MATRIX="recall_${CHUNK_ID}.txt"
    local F1_MATRIX="f1_${CHUNK_ID}.txt"
    while read VCF_FILE; do
        MEASURED_FILE="${VCF_FILE%.vcf}_filtered.vcf"
        TRUE_FILE="true_${CHUNK_ID}.vcf"
        OUTPUT_DIR="ouput_${CHUNK_ID}"
        if [ ${#FILTER_STRING} -ne 0 ]; then
            bcftools filter --include '${FILTER_STRING}' ${VCF_FILE} > ${MEASURED_FILE}
        else 
            cp ${VCF_FILE} ${MEASURED_FILE}
        fi
        bgzip -@ 1 ${MEASURED_FILE}
        tabix ${MEASURED_FILE}.gz
        ID=$(basename ${VCF_FILE} .vcf)
        ID=${ID#${CALLER}_i}
        ID=${ID%_i*_l${READ_LENGTH}_c${COVERAGE}}
        if [ ${#FILTER_STRING} -ne 0 ]; then
            bcftools filter --include '${FILTER_STRING}' ground_truth_vcfs/groundTruth_individual_${ID}.vcf > ${TRUE_FILE}
        else
            cp ground_truth_vcfs/groundTruth_individual_${ID}.vcf ${TRUE_FILE}
        fi
        bgzip -@ 1 ${TRUE_FILE}
        tabix ${TRUE_FILE}.gz
        ${TIME_COMMAND} truvari bench -b ${TRUE_FILE}.gz -c ${MEASURED_FILE}.gz -f ${REFERENCE_FA} -o ${OUTPUT_DIR}/
        rm -f ${TRUE_FILE}*
        grep "\"TP-call\":" ${OUTPUT_DIR}/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${TP_MATRIX}
        grep "\"FP\":" ${OUTPUT_DIR}/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${FP_MATRIX}
        grep "\"FN\":" ${OUTPUT_DIR}/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${FN_MATRIX}
        grep "\"precision\":" ${OUTPUT_DIR}/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${PRECISION_MATRIX}
        grep "\"recall\":" ${OUTPUT_DIR}/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${RECALL_MATRIX}
        grep "\"f1\":" ${OUTPUT_DIR}/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${F1_MATRIX}
        rm -rf ${OUTPUT_DIR}/
    done < chunk-${CHUNK_ID}
}


# 1. Per-individual matrices
find measured_vcfs/ -maxdepth 1 -name '*_i*_i*_l*_c*.vcf' > files.txt
split -d -n ${N_THREADS} files.txt chunk-
rm -f files.txt
for CHUNK_FILE in $(ls chunk-*); do
    processChunk ${CHUNK_FILE#chunk-} &
done
wait
for CHUNK_FILE in $(ls chunk-*); do
    PREFIX=${CHUNK_FILE#chunk-}
    cat tp_${PREFIX}.txt >> ${TP_MATRIX}
    cat fp_${PREFIX}.txt >> ${FP_MATRIX}
    cat fn_${PREFIX}.txt >> ${FN_MATRIX}
    cat precision_${PREFIX}.txt >> ${PRECISION_MATRIX}
    cat recall_${PREFIX}.txt >> ${RECALL_MATRIX}
    cat f1_${PREFIX}.txt >> ${F1_MATRIX}
done
rm -f chunk-* tp_*.txt fp_*.txt fn_*.txt precision_*.txt recall_*.txt f1_*.txt

# 2. All-individuals matrices: merging filtered calls over all individuals, and
# comparing the merge to the filtered set of distinct variants in the ground
# truth.
rm -f files.txt
find measured_vcfs/ -maxdepth 1 -name '*_filtered.vcf.gz' > files.txt
${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --output-type z --merge none --file-list file.txt --output merge.vcf.gz 
truvari collapse --threads ${N_THREADS} --keep common --minhaplen 40 --sizemin 40 --input merge.vcf.gz --output truvari_merge.vcf --collapsed-output truvari_collapsed.vcf --reference ${REFERENCE_FA}
bgzip -@ ${N_THREADS} truvari_merge.vcf
tabix truvari_merge.vcf.gz
gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp truvari_merge.vcf.gz ${BUCKET_DIR_ALLINDIVIDUALS_VCFS}/${CONFIGURATION_ID}_mergedIndividuals.vcf.gz
if [ ${#FILTER_STRING} -ne 0 ]; then
    ${TIME_COMMAND} bcftools filter --include '${FILTER_STRING}' ground_truth_vcfs/groundTruth_joint.vcf > true_filtered.vcf
else 
    cp ground_truth_vcfs/groundTruth_joint.vcf true_filtered.vcf
fi
bgzip -@ ${N_THREADS} true_filtered.vcf
tabix true_filtered.vcf.gz
${TIME_COMMAND} truvari bench -b true_filtered.vcf.gz -c truvari_merge.vcf.gz -f ${REFERENCE_FA} -o ouput/
grep "\"TP-call\":" output/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${TP_MATRIX_MERGE}
grep "\"FP\":" output/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${FP_MATRIX_MERGE}
grep "\"FN\":" output/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${FN_MATRIX_MERGE}
grep "\"precision\":" output/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${PRECISION_MATRIX_MERGE}
grep "\"recall\":" output/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${RECALL_MATRIX_MERGE}
grep "\"f1\":" output/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${F1_MATRIX_MERGE}
rm -rf truvari_merge.vcf* ouput/

# 3. Joint-calling matrices (if any): filtering the joint calling file with
# $FILTER_STRING$, and comparing it to the filtered set of distinct variants in
# the ground truth.
if [ ${JOINT_CALLING_FILE} != null ]; then
    if [ ${#FILTER_STRING} -ne 0 ]; then
        ${TIME_COMMAND} bcftools filter --include '${FILTER_STRING}' ${JOINT_CALLING_FILE} > joint_filtered.vcf
    else 
        cp ${JOINT_CALLING_FILE} joint_filtered.vcf
    fi
    bgzip -@ ${N_THREADS} joint_filtered.vcf
    tabix joint_filtered.vcf.gz
    ${TIME_COMMAND} truvari bench -b true_filtered.vcf.gz -c joint_filtered.vcf.gz -f ${REFERENCE_FA} -o ouput/
    grep "\"TP-call\":" output/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${TP_MATRIX_JOINT}
    grep "\"FP\":" output/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${FP_MATRIX_JOINT}
    grep "\"FN\":" output/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${FN_MATRIX_JOINT}
    grep "\"precision\":" output/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${PRECISION_MATRIX_JOINT}
    grep "\"recall\":" output/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${RECALL_MATRIX_JOINT}
    grep "\"f1\":" output/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${F1_MATRIX_JOINT}
    rm -rf joint_filtered.vcf* ouput/
fi
rm -rf true_filtered.vcf*
