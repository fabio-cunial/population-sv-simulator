#!/bin/bash
#
REPEATMASKER_FILE="../buildModel/input/GRCh37/repeats/hg19.fa.out.cleaned"
REPEATMASKER_FILE_NROWS="$(wc -l < ${REPEATMASKER_FILE})"
TRF_FILE="../buildModel/input/GRCh37/repeats/hg19.trf.bed.cleaned"
TRF_FILE_NROWS="$(wc -l < ${TRF_FILE})"
SEGDUPS_FILE="../buildModel/input/GRCh37/repeats/segdups.txt"
SEGDUPS_FILE_NROWS="$(( $(wc -l < ${SEGDUPS_FILE}) - 1))"
POPULATION_SV_FILE="./output/groundTruth_joint.vcf"
ONLY_PASS="1"
OUTPUT_DIR="./output/check"
VERBOSE="1"

rm -rf ${OUTPUT_DIR}; mkdir -p ${OUTPUT_DIR}
java -cp ../src BuildModel ${REPEATMASKER_FILE} ${REPEATMASKER_FILE_NROWS} ${TRF_FILE} ${TRF_FILE_NROWS} ${SEGDUPS_FILE} ${SEGDUPS_FILE_NROWS} ${POPULATION_SV_FILE} ${ONLY_PASS} ${OUTPUT_DIR} ${VERBOSE}
