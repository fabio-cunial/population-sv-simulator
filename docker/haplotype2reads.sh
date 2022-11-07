#!/bin/bash
#
# Creates a diploid individual by taking two haplotypes with consecutive IDs.
# Runs the read simulator at max coverage on each. Randomly permutes its
# output, to make sure the two haplotypes are mixed when splitting the file
# later.
#
# Resource analysis for 20x coverage of one haplotype (~10 GB). Intel Xeon,
# 2.30GHz, 8 physical cores.
#
# TASK  | TIME   | RAM    | CORES | COMMENT
# pbsim | 7 mins | 250 MB |   1   |
# shuf  | 1.5 m  | 10 GB  |   1   | loads the whole file in memory
# bucket upload/download: ~3 m
#
ID1=$1
ID2=$2
LENGTH_MIN=$3
LENGTH_MAX=$4
LENGTH_MEAN=$5
LENGTH_STDEV=$6
MAX_COVERAGE=$7
BUCKET_DIR=$8  # Root dir of the simulation in the bucket
WORK_DIR=$9  # Absolute path

GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
GSUTIL_DELAY_S="600"
TIME_COMMAND="/usr/bin/time --verbose"
PBSIM1_MODEL="/simulation/pbsim_data/model_qc_ccs"
PBSIM2_MODEL="/simulation/pbsim_data/P6C4.model"
PBSIM1_COMMAND="pbsim --data-type CLR --length-min ${LENGTH_MIN} --length-max ${LENGTH_MAX} --accuracy-min 0.99 --accuracy-max 1.0 --difference-ratio 6:21:73 --model_qc  ${PBSIM1_MODEL} --length-mean ${LENGTH_MEAN} --length-sd ${LENGTH_STDEV} --accuracy-mean 0.995 --accuracy-sd 0.002"
PBSIM2_COMMAND="pbsim                 --length-min ${LENGTH_MIN} --length-max ${LENGTH_MAX} --accuracy-min 0.99 --accuracy-max 1.0 --difference-ratio 6:50:54 --hmm_model ${PBSIM2_MODEL} --length-mean ${LENGTH_MEAN} --length-sd ${LENGTH_STDEV} --accuracy-mean 0.995"

set -euxo pipefail
echo "Running <haplotype2reads.sh> on the following node:"
lscpu
cat /proc/meminfo
cd ${WORK_DIR}

OUTPUT_FILE="reads_i${ID1}_i${ID2}_l${LENGTH_MEAN}_c${MAX_COVERAGE}.fa"
TEST=$(gsutil -q stat ${BUCKET_DIR}/reads/${OUTPUT_FILE} && echo 0 || echo 1)
if [ ${TEST} -eq 0 ]; then
    while : ; do
        TEST=$(gsutil cp ${BUCKET_DIR}/reads/${OUTPUT_FILE} . && echo 0 || echo 1)
        if [ ${TEST} -eq 1 ]; then
            echo "Error downloading file <${BUCKET_DIR}/reads/${OUTPUT_FILE}>. Trying again..."
            sleep ${GSUTIL_DELAY_S}
        else
            break
        fi
    done
else
    ${TIME_COMMAND} ${PBSIM2_COMMAND} --depth ${MAX_COVERAGE} --prefix reads_${ID1} haplotype_${ID1}.fa
    mv reads_${ID1}_0001.fastq reads_${ID1}.fa
    rm -f reads_${ID1}_0001.* 
    ${TIME_COMMAND} ${PBSIM2_COMMAND} --depth ${MAX_COVERAGE} --prefix reads_${ID2} haplotype_${ID2}.fa
    mv reads_${ID2}_0001.fastq reads_${ID2}.fa
    rm -f reads_${ID2}_0001.* 
    cat reads_${ID2}.fa >> reads_${ID1}.fa
    rm -f reads_${ID2}.fa
    ${TIME_COMMAND} shuf reads_${ID1}.fa > reads_${ID1}_${ID2}.1
    rm -f reads_${ID1}.fa
    tr 'Z' '\n' < reads_${ID1}_${ID2}.1 > ${OUTPUT_FILE}
    rm -f reads_${ID1}_${ID2}.1
    while : ; do
        TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp ${OUTPUT_FILE} ${BUCKET_DIR}/reads/ && echo 0 || echo 1)
        if [ ${TEST} -eq 1 ]; then
            echo "Error uploading file <${OUTPUT_FILE}>. Trying again..."
            sleep ${GSUTIL_DELAY_S}
        else
            break
        fi
    done
fi
