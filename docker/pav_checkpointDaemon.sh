#!/bin/bash
#
# Periodically checks for new PAV state files and checkpoints them to a remote
# bucket, in an infinite loop.
#
# Remark: this is not the perfect strategy for checkpointing a Snakemake script:
# it is just a simple heuristic to avoid wasting a lot of time if the machine
# gets preempted while running the script.
#
BUCKET_DIR=$1  # Root dir of the simulation in the bucket
SAMPLE_ID=$2
WORK_DIR=$3  # Absolute path
DOCKER_DIR=$4  # Absolute path

STATE_FILES_LIST="${DOCKER_DIR}/pav_state_files.txt"
IN_BUCKET_FILE="pav_in_bucket.txt"
DELAY_S="60"
GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
SUFFIX="checkpointed"

# <set -u> has to be issued after array construction, otherwise it gives unbound
# variable error on the array.
set -eo pipefail
cd ${WORK_DIR}

STATE_ARRAY=()
while read LINE; do
    STATE_ARRAY+=(${LINE})
done < ${IN_BUCKET_FILE}
echo "STATE_ARRAY=${STATE_ARRAY[*]}"
set -u
rm -f *.${SUFFIX}
while true; do
    tree data/ temp/ results/ || echo ""
    i=0
    while read FILE; do
        echo "Checking if file exists: ${FILE}"
        if [ ${STATE_ARRAY[${i}]} -eq 0 -a -e ${FILE} ]; then
            echo "FILE exists! ${FILE} lastModified=$(stat --format %Y ${FILE}) currentTime=$(date +%s) DELAY_S=${DELAY_S}"
            if [ $(date +%s) -ge $(( $(stat --format %Y ${FILE}) + ${DELAY_S} )) ]; then
                cp ${FILE} $(echo ${FILE} | tr '/' '+').${SUFFIX}
                STATE_ARRAY[${i}]="1"
                echo "STATE_ARRAY[${i}] set to one! FILE=${FILE}"
            fi
        fi
        i=$((${i} + 1))
    done < ${STATE_FILES_LIST}
    while read FILE; do
        FILE_PRIME=$(echo ${FILE} | tr '/' '+').${SUFFIX}
        if [ -e ${FILE_PRIME} ]; then
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp ${FILE_PRIME} ${BUCKET_DIR}/pav/${SAMPLE_ID}/${FILE} && echo 0 || echo 1)
            if [ ${TEST} -eq 0 ]; then
                rm -f ${FILE_PRIME}
            fi
        fi
    done < ${STATE_FILES_LIST}
    sleep ${DELAY_S}
done
