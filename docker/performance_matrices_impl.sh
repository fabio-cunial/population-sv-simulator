#!/bin/bash
#
# Given a set of rules for filtering VCFs, the program: 
#
# 1. Compares every filtered VCF file in <experimental_vcfs/> to the
# corresponding file in <ground_truth_vcfs/>, in parallel over all the cores of
# the machine, collecting per-individual statistics.
#
# 2. Merges all filtered VCF files in <experimental_vcfs/>, and compares the
# result to the filtered set of distinct ground-truth variants in the
# population.
#
# 3. Compares a filtered joint-calling file to the filtered set of distinct
# ground-truth variants in the population.
#
# Remark: genotypes are not considered. BNDs are not interpreted (this might
# penalize some callers).
#
# Resource analysis. Intel Xeon, 2.30GHz, 16 physical cores.
#
# TASK                    | TIME   | RAM    | CORES | COMMENT
# bcftools merge          | 1 s    | 100 MB |   1   |
# truvari collapse        | 2 s    | 100 MB |  3-4  |
# bcftools sort merged    | 1 s    | 10 MB  |   1   |
# bcftools filter truth   | 1 s    | 20 MB  |   1   |
# truvari bench merged    | 1 s    | 100 MB |  4.5  |
#
FILTER_STRING=$1  # Spaces are replaced with '+'. Equals '+' if empty.
FILTER_STRING_FREQUENCY=$2  # Spaces are replaced with '+'. Equals '+' if empty.
JOINT_CALLING_FILE=$3  # Relative path or "null".
REFERENCE_FA_ADDRESS=$4  # Address in a bucket
REFERENCE_FAI_ADDRESS=$5  # Address in a bucket
N_THREADS=$6
WORK_DIR=$7  # Absolute path
CONFIGURATION_ID=$8  # Used just for storing the all-individuals merged VCF
CALLER=$9
READ_LENGTH=${10}
COVERAGE=${11}
BUCKET_DIR_ALLINDIVIDUALS_VCFS=${12}

# Per-individual matrices. The script appends many values (one value per
# individual) to a single row of each matrix.
TP_MATRIX=${13}
FP_MATRIX=${14}
FN_MATRIX=${15}
PRECISION_MATRIX=${16}
RECALL_MATRIX=${17}
F1_MATRIX=${18}

# All-individuals matrices. The script appends one aggregated value to a single
# row of each matrix.
TP_MATRIX_MERGE=${19}
FP_MATRIX_MERGE=${20}
FN_MATRIX_MERGE=${21}
PRECISION_MATRIX_MERGE=${22}
RECALL_MATRIX_MERGE=${23}
F1_MATRIX_MERGE=${24}

# Joint-calling matrices (ignored if there is no joint calling file). The
# script appends one aggregated value to a single row of each matrix.
TP_MATRIX_JOINT=${25}
FP_MATRIX_JOINT=${26}
FN_MATRIX_JOINT=${27}
PRECISION_MATRIX_JOINT=${28}
RECALL_MATRIX_JOINT=${29}
F1_MATRIX_JOINT=${30}


set -euxo pipefail
cd ${WORK_DIR}
GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
GSUTIL_DELAY_S="600"
TIME_COMMAND="/usr/bin/time --verbose"
if [ ${FILTER_STRING} = "+" ]; then
    FILTER_STRING=""
else
    FILTER_STRING=$(echo ${FILTER_STRING} | tr '+' ' ')
fi
if [ ${FILTER_STRING_FREQUENCY} = "+" ]; then
    FILTER_STRING_FREQUENCY=""
else
    FILTER_STRING_FREQUENCY=$(echo ${FILTER_STRING_FREQUENCY} | tr '+' ' ')
fi
FILTER_STRING_TRUTH=${FILTER_STRING}
if [ ${#FILTER_STRING_FREQUENCY} -ne 0  ]; then
    if [ ${#FILTER_STRING_TRUTH} -ne 0  ]; then
        FILTER_STRING_TRUTH="${FILTER_STRING_TRUTH} && ${FILTER_STRING_FREQUENCY}"
    else
        FILTER_STRING_TRUTH=${FILTER_STRING_FREQUENCY}
    fi
fi
BCFTOOLS_MERGE_FLAGS="--force-samples --merge none"
TRUVARI_BENCH_FLAGS=" "  # Default settings for now
TRUVARI_COLLAPSE_FLAGS="--keep common --sizemin 40"
rm -f new_headers.txt
echo "##INFO=<ID=REPEATS_START,Number=1,Type=Integer,Description=\"Repetitive context around the first position\">" >> new_headers.txt
echo "##INFO=<ID=REPEATS_END,Number=1,Type=Integer,Description=\"Repetitive context around the last position\">" >> new_headers.txt
echo "##INFO=<ID=REPEATS_FRACTION,Number=1,Type=Float,Description=\"Fraction of the SV covered by repeats of any type\">" >> new_headers.txt


# Ensures that a VCF file contains all and only the basic columns.
#
function removeVCFColumns() {
    local INPUT_FILE=$1
    
    bcftools view --threads 0 -h ${INPUT_FILE} > ${INPUT_FILE}_header.txt
    N_ROWS=$(wc -l < ${INPUT_FILE}_header.txt)
    head -n $((${N_ROWS} - 1)) ${INPUT_FILE}_header.txt > ${INPUT_FILE}_new.vcf
    rm -f ${INPUT_FILE}_header.txt
    echo "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" >> ${INPUT_FILE}_new.vcf
    echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE" >> ${INPUT_FILE}_new.vcf
    bcftools view --threads 0 -H ${INPUT_FILE} | awk '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tGT\t.\n", $1, $2, $3, $4, $5, $6, $7, $8)}' >> ${INPUT_FILE}_new.vcf
    rm -f ${INPUT_FILE}
    cat ${INPUT_FILE}_new.vcf | sed '/contig=<ID=0>/d' | sed '/bcftools_/d' > ${INPUT_FILE}
    rm -f ${INPUT_FILE}_new.vcf
}


# Appends <new_headers.txt> to a VCF file. This is needed by <bcftools filter>.
#
function reheader() {
    local INPUT_FILE=$1
    
    bcftools view -h ${INPUT_FILE} > ${INPUT_FILE}_headers_old.txt
    N_LINES=$(wc -l < ${INPUT_FILE}_headers_old.txt)
    head -n $((${N_LINES} - 1)) ${INPUT_FILE}_headers_old.txt > ${INPUT_FILE}_headers_new.txt
    cat new_headers.txt >> ${INPUT_FILE}_headers_new.txt
    tail -n 1 ${INPUT_FILE}_headers_old.txt >> ${INPUT_FILE}_headers_new.txt
    bcftools reheader -h ${INPUT_FILE}_headers_new.txt ${INPUT_FILE} | sed '/contig=<ID=0>/d' | sed '/bcftools_/d' > ${INPUT_FILE}.newHeaders
    rm -f ${INPUT_FILE}
    mv ${INPUT_FILE}.newHeaders ${INPUT_FILE}
    rm -f ${INPUT_FILE}_headers_old.txt ${INPUT_FILE}_headers_new.txt
}


# Filters every VCF file in a chunk of files using $FILTER_STRING$, filters the
# corresponding ground truth as well, and compares the two, sequentially.
#
# Remark: we sort experimental VCFs since there is no guarantee they are sorted.
# Ground truth VCFs are assumed to be sorted.
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
        MEASURED_FILE="${VCF_FILE%.vcf}_filtered.vcf.gz"
        TRUE_FILE="true_${CHUNK_ID}.vcf.gz"
        OUTPUT_DIR="output_${CHUNK_ID}"
        if [ ${#FILTER_STRING} -ne 0 ]; then
            reheader ${VCF_FILE}
            bcftools filter --threads 0 --include "${FILTER_STRING}" --output-type v ${VCF_FILE} | bcftools sort --output-type z --output ${MEASURED_FILE} 
        else
            bcftools sort --output-type z --output ${MEASURED_FILE} ${VCF_FILE}
        fi
        tabix ${MEASURED_FILE}
        ID=$(basename ${VCF_FILE} .vcf)
        ID=${ID#${CALLER}_i}
        ID=${ID%_i*_l${READ_LENGTH}_c${COVERAGE}_annotated}
        if [ ${#FILTER_STRING_TRUTH} -ne 0 ]; then
            bcftools filter --threads 0 --include "${FILTER_STRING_TRUTH}" --output-type z --output ${TRUE_FILE} ground_truth_vcfs/groundTruth_individual_${ID}.vcf
        else
            bgzip --threads 1 --stdout ground_truth_vcfs/groundTruth_individual_${ID}.vcf > ${TRUE_FILE}
        fi
        tabix ${TRUE_FILE}
        truvari bench ${TRUVARI_BENCH_FLAGS} --prog --base ${TRUE_FILE} --comp ${MEASURED_FILE} --reference reference.fa --output ${OUTPUT_DIR}/
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


# Main program
# 1. Per-individual matrices
find experimental_vcfs/ -maxdepth 1 -name '*_i*_i*_l*_c*_annotated.vcf' > tmp.txt
shuf tmp.txt > files.txt  # For better balancing
rm -f tmp.txt
N_LINES_PER_CHUNK="0"
N_LINES=$(wc -l < files.txt)
if [ ${N_THREADS} -ge ${N_LINES} ]; then
    N_LINES_PER_CHUNK="1"
else
    N_LINES_PER_CHUNK=$(( ${N_LINES} / ${N_THREADS} ))
fi
split -d -l ${N_LINES_PER_CHUNK} files.txt chunk-
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
#
# Remark: we sort the experimental VCF to be absolutely sure it is correct
# (this might be unnecessary). The ground truth VCF is assumed to be sorted.
#
TEST=$(gsutil -q stat ${BUCKET_DIR_ALLINDIVIDUALS_VCFS}/${CONFIGURATION_ID}_mergedIndividuals.vcf.gz && echo 0 || echo 1)
if [ ${TEST} -eq 0 ]; then
    while : ; do
        TEST=$(gsutil cp ${BUCKET_DIR_ALLINDIVIDUALS_VCFS}/${CONFIGURATION_ID}_mergedIndividuals.vcf.gz truvari_merge.vcf.gz && echo 0 || echo 1)
        if [ ${TEST} -eq 1 ]; then
            echo "Error downloading file <${BUCKET_DIR_ALLINDIVIDUALS_VCFS}/${CONFIGURATION_ID}_mergedIndividuals.vcf.gz>. Trying again..."
            sleep ${GSUTIL_DELAY_S}
        else
            break
        fi
    done
else
    rm -f files.txt
    find experimental_vcfs/ -maxdepth 1 -name '*_filtered.vcf.gz' > files.txt
    ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} ${BCFTOOLS_MERGE_FLAGS} --file-list files.txt --output-type z --output merge.vcf.gz
    tabix merge.vcf.gz
    ${TIME_COMMAND} truvari collapse --threads ${N_THREADS} ${TRUVARI_COLLAPSE_FLAGS} --input merge.vcf.gz --output truvari_merge.vcf --collapsed-output truvari_collapsed.vcf --reference reference.fa
    removeVCFColumns truvari_merge.vcf
    ${TIME_COMMAND} bcftools sort --output-type z --output truvari_merge.vcf.gz truvari_merge.vcf
    while : ; do
        TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp truvari_merge.vcf.gz ${BUCKET_DIR_ALLINDIVIDUALS_VCFS}/${CONFIGURATION_ID}_mergedIndividuals.vcf.gz && echo 0 || echo 1)
        if [ ${TEST} -eq 1 ]; then
            echo "Error uploading file <truvari_merge.vcf.gz>. Trying again..."
            sleep ${GSUTIL_DELAY_S}
        else
            break
        fi
    done
fi
tabix truvari_merge.vcf.gz
if [ ${#FILTER_STRING_TRUTH} -ne 0 ]; then
    ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include "${FILTER_STRING_TRUTH}" --output-type v ground_truth_vcfs/groundTruth_joint.vcf.gz | sed '/contig=<ID=0>/d' | sed '/bcftools_/d' > true_filtered.vcf
    bgzip --threads ${N_THREADS} true_filtered.vcf
else 
    cp ground_truth_vcfs/groundTruth_joint.vcf.gz true_filtered.vcf.gz
fi
tabix true_filtered.vcf.gz
${TIME_COMMAND} truvari bench ${TRUVARI_BENCH_FLAGS} --prog --base true_filtered.vcf.gz --comp truvari_merge.vcf.gz --reference reference.fa --output output/    
grep "\"TP-call\":" output/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${TP_MATRIX_MERGE}
grep "\"FP\":" output/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${FP_MATRIX_MERGE}
grep "\"FN\":" output/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${FN_MATRIX_MERGE}
grep "\"precision\":" output/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${PRECISION_MATRIX_MERGE}
grep "\"recall\":" output/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${RECALL_MATRIX_MERGE}
grep "\"f1\":" output/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${F1_MATRIX_MERGE}
rm -rf truvari_merge.vcf* output/

# 3. Joint-calling matrices (if any): filtering the joint calling file and
# comparing it to the filtered set of distinct variants in the ground truth.
#
# Remark: we sort the experimental VCF since there is no guarantee that it is
# sorted. The ground truth VCF is assumed to be sorted.
#
if [ ${JOINT_CALLING_FILE} != "null" ]; then
    removeVCFColumns ${JOINT_CALLING_FILE}
    if [ ${#FILTER_STRING} -ne 0 ]; then
        reheader ${JOINT_CALLING_FILE}
        bcftools filter --threads ${N_THREADS} --include "${FILTER_STRING}" --output-type v ${JOINT_CALLING_FILE} | bcftools sort --output-type z --output joint_filtered.vcf.gz
    else 
        bcftools sort --output-type z --output joint_filtered.vcf.gz ${JOINT_CALLING_FILE}
    fi
    tabix joint_filtered.vcf.gz
    ${TIME_COMMAND} truvari bench ${TRUVARI_BENCH_FLAGS} --base true_filtered.vcf.gz --comp joint_filtered.vcf.gz --reference reference.fa --output output/
    grep "\"TP-call\":" output/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${TP_MATRIX_JOINT}
    grep "\"FP\":" output/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${FP_MATRIX_JOINT}
    grep "\"FN\":" output/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${FN_MATRIX_JOINT}
    grep "\"precision\":" output/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${PRECISION_MATRIX_JOINT}
    grep "\"recall\":" output/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${RECALL_MATRIX_JOINT}
    grep "\"f1\":" output/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${F1_MATRIX_JOINT}
    rm -rf joint_filtered.vcf* output/
fi
rm -rf true_filtered.vcf*
