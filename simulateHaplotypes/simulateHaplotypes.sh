#!/bin/bash
#
MODEL_DIR="./input"
POPULATION_NAME="ALL"
REFERENCE_FILE="./input/T2T/hs1_chr1.fa"
REPEAT_INTERVALS_FILE="./input/T2T/repeatIntervals.txt"
SEGDUPS_FILE="./input/T2T/segdups.txt"
INSERTION_STRINGS_FILE="./input/insertions/union_170509_refalt_gt49bp.sort.INS.txt"
N_INSERTION_STRINGS="$(wc -l < ${INSERTION_STRINGS_FILE})"
ALLELE_FREQUENCY_THRESHOLD="0.01"
N_VARIANTS_FREQUENT="5000"
N_INDIVIDUALS="300"   #"10000"
OUTPUT_DIR="./output"

java -Xmx10g -cp ../src SimulateHaplotypes ${MODEL_DIR} ${POPULATION_NAME} ${REFERENCE_FILE} ${REPEAT_INTERVALS_FILE} ${SEGDUPS_FILE} ${INSERTION_STRINGS_FILE} ${N_INSERTION_STRINGS} ${ALLELE_FREQUENCY_THRESHOLD} ${N_VARIANTS_FREQUENT} ${N_INDIVIDUALS} ${OUTPUT_DIR}
