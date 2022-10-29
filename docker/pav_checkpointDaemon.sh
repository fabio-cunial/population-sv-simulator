#!/bin/bash
#
# Periodically checks for new PAV state files and checkpoints them to a remote
# bucket, in an infinite loop.
#
BUCKET_DIR=$1  # Root dir of the simulation in the bucket
SAMPLE_ID=$2
WORK_DIR=$3  # Local directory

STATE_FILES_LIST="pav_state_files.txt"
IN_BUCKET_FILE="pav_in_bucket.txt"
DELAY_S="300"
DELAY_MS=$((${DELAY_S} * 1000))
GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"

set -euxo pipefail

cd ${WORK_DIR}
STATE_ARRAY=()
while read LINE; do
    STATE_ARRAY+=( ${LINE} )
done < ${IN_BUCKET_FILE}
currentTime=$(date +%s)
while true; do
    i=0
    while read FILE; do
        lastModified=$(stat --format %Y ${FILE})
        if [ ${STATE_ARRAY[${i}]} -eq 0 -a -e ${FILE} -a ${currentTime} -ge $((${lastModified} + ${DELAY_MS})) ]; then
            gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp ${FILE} ${BUCKET_DIR}/pav/${SAMPLE_ID}/${FILE}
            STATE_ARRAY[${i}]="1"
        fi
        i=$((${i} + 1))
    done < ${STATE_FILES_LIST}
    sleep ${DELAY_S}s
done
