#!/bin/bash
REPEATMASKER_FILE="../build_model/input/repeats/hg19.fa.out.cleaned"
REPEATMASKER_FILE_NROWS="$(wc -l < ${REPEATMASKER_FILE})"
TRF_FILE="../build_model/input/repeats/hg19.trf.bed.cleaned"
TRF_FILE_NROWS="$(wc -l < ${TRF_FILE})"
POPULATION_SV_FILE="./output/simulatedVariants.vcf"
ONLY_PASS="1"
OUTPUT_DIR="./output/check"
VERBOSE="1"

rm -rf ${OUTPUT_DIR}; mkdir -p ${OUTPUT_DIR}
java -cp ./src BuildModel ${REPEATMASKER_FILE} ${REPEATMASKER_FILE_NROWS} ${TRF_FILE} ${TRF_FILE_NROWS} ${POPULATION_SV_FILE} ${ONLY_PASS} ${OUTPUT_DIR} ${VERBOSE}
