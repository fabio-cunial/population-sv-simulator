#!/bin/bash
#
# Runs Truvari between every VCF file in <measured_vcfs/> and the corresponding
# file in <ground_truth_vcfs/>, in parallel over all the cores of the machine.
#
FILTER_STRING=$1  # Spaces are replaced with '+'. It equals '+' if empty.
REFERENCE_FA=$2
N_THREADS=$3
WORK_DIR=$4  # Absolute path
TP_MATRIX=$5
FP_MATRIX=$6
FN_MATRIX=$7
PRECISION_MATRIX=$8
RECALL_MATRIX=$9
F1_MATRIX=${10}

set -euxo pipefail
cd ${WORK_DIR}
TIME_COMMAND="/usr/bin/time --verbose"
if [ ${FILTER_STRING} = "+" ]; then
    FILTER_STRING=""
else
    FILTER_STRING=$(echo ${FILTER_STRING} | tr '+' ' ')
fi


# Handles a chunk of individuals sequentially.
function processChunk() {
    local CHUNK_ID=$1
    
    local TP_MATRIX="${CHUNK_ID}_tp.txt"
    local FP_MATRIX="${CHUNK_ID}_fp.txt"
    local FN_MATRIX="${CHUNK_ID}_fn.txt"
    local PRECISION_MATRIX="${CHUNK_ID}_precision.txt"
    local RECALL_MATRIX="${CHUNK_ID}_recall.txt"
    local F1_MATRIX="${CHUNK_ID}_f1.txt"
    while read VCF_FILE; do
        if [ ${#FILTER_STRING} -ne 0 ]; then
            bcftools filter --include '${FILTER_STRING}' ${VCF_FILE} > measured.vcf
        else 
            cp ${VCF_FILE} measured.vcf
        fi
        bgzip -@ 1 measured.vcf
        tabix measured.vcf.gz
        ID=$(basename ${VCF_FILE} .vcf)
        ID=${ID#${CALLER}_i}
        ID=${ID%_i*_l${READ_LENGTH}_c${COVERAGE}}
        if [ ${#FILTER_STRING} -ne 0 ]; then
            bcftools filter --include '${FILTER_STRING}' ground_truth_vcfs/groundTruth_individual_${ID}.vcf > true.vcf
        else
            cp ground_truth_vcfs/groundTruth_individual_${ID}.vcf true.vcf
        fi
        bgzip -@ 1 true.vcf
        tabix true.vcf.gz
        ${TIME_COMMAND} truvari bench -b true.vcf.gz -c measured.vcf.gz -f ${REFERENCE_FA} -o ouput/
        grep "\"TP-call\":" output/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${TP_MATRIX}
        grep "\"FP\":" output/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${FP_MATRIX}
        grep "\"FN\":" output/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${FN_MATRIX}
        grep "\"precision\":" output/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${PRECISION_MATRIX}
        grep "\"recall\":" output/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${RECALL_MATRIX}
        grep "\"f1\":" output/summary.txt | awk 'BEGIN {ORS=""} {print $2}' >> ${F1_MATRIX}
        rm -rf measured.vcf* true.vcf* ouput/
    done < chunk-${CHUNK_ID}
}


# Main program
find measured_vcfs/ -maxdepth 1 -name '*.vcf' > files.txt
split -d -n ${N_THREADS} files.txt chunk-
rm -f files.txt
for CHUNK_FILE in $(ls chunk-*); do
    processChunk ${CHUNK_FILE#chunk-} &
done
wait
for CHUNK_FILE in $(ls chunk-*); do
    PREFIX=${CHUNK_FILE#chunk-}
    cat ${PREFIX}_tp.txt >> ${TP_MATRIX}
    cat ${PREFIX}_fp.txt >> ${FP_MATRIX}
    cat ${PREFIX}_fn.txt >> ${FN_MATRIX}
    cat ${PREFIX}_precision.txt >> ${PRECISION_MATRIX}
    cat ${PREFIX}_recall.txt >> ${RECALL_MATRIX}
    cat ${PREFIX}_f1.txt >> ${F1_MATRIX}
done
rm -f chunk-*.txt *_tp.txt *_fp.txt *_fn.txt *_precision.txt *_recall.txt *_f1.txt







bcftools merge --threads ${N_THREADS} -m none --file-list file.txt | bgzip > merge.vcf.gz
truvari collapse --threads ${N_THREADS} --keep common --minhaplen 40 --sizemin 40 -i merge.vcf.gz -o truvari_merge.vcf -c truvari_collapsed.vcf --reference ${REFERENCE_FA} 