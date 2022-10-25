#!/bin/bash
#
REPEATMASKER_FILE="./input/GRCh37/repeats/hg19.fa.out.cleaned"
REPEATMASKER_FILE_NROWS="$(wc -l < ${REPEATMASKER_FILE})"
TRF_FILE="./input/GRCh37/repeats/hg19.trf.bed.cleaned"
TRF_FILE_NROWS="$(wc -l < ${TRF_FILE})"
SEGDUPS_FILE="./input/GRCh37/repeats/segdups.txt"
SEGDUPS_FILE_NROWS="$(( $(wc -l < ${SEGDUPS_FILE}) - 1))"
GNOMADSV_FILE="./input/GRCh37/gnomad_v2.1_sv.sites.vcf"
ONLY_PASS="1"
OUTPUT_DIR="./output"
VERBOSE="1"

java -cp ../src BuildModel ${REPEATMASKER_FILE} ${REPEATMASKER_FILE_NROWS} ${TRF_FILE} ${TRF_FILE_NROWS} ${SEGDUPS_FILE} ${SEGDUPS_FILE_NROWS} ${GNOMADSV_FILE} ${ONLY_PASS} ${OUTPUT_DIR} ${VERBOSE}
