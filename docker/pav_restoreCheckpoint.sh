#!/bin/bash
#
# Downloads all the PAV state files that were checkpointed to a remote bucket.
#
BUCKET_DIR=$1  # Root dir of the simulation in the bucket
SAMPLE_ID=$2
WORK_DIR=$3  # Absolute path
DOCKER_DIR=$4  # Absolute path

STATE_FILES_LIST="${DOCKER_DIR}/pav_state_files.txt"
IN_BUCKET_FILE="pav_in_bucket.txt"

set -euxo pipefail
cd ${WORK_DIR}

# Initializing local directory structure
mkdir -p data/ref/
mkdir -p temp/sample/align/pre-cut/
mkdir -p results/sample/align/pre-cut/
mkdir -p temp/sample/cigar/batched/
mkdir -p temp/sample/lg_sv/
mkdir -p temp/sample/cigar/pre_inv/
mkdir -p temp/sample/inv_caller/flag/
mkdir -p results/sample/callable/
mkdir -p results/sample/inv_caller/
mkdir -p temp/sample/inv_caller/batch/
mkdir -p temp/sample/bed/integrated/h1/
mkdir -p temp/sample/bed/integrated/h2/
mkdir -p temp/sample/bed/bychrom/svindel_ins/
mkdir -p temp/sample/bed/bychrom/svindel_del/
mkdir -p temp/sample/bed/bychrom/sv_inv/
mkdir -p temp/sample/bed/bychrom/snv_snv/
mkdir -p temp/sample/bed/merged/
mkdir -p results/sample/bed/
mkdir -p results/sample/bed/fa/

# Restoring files
rm -f ${IN_BUCKET_FILE}
while read FILE; do
    TEST=$(gsutil -q stat ${BUCKET_DIR}/pav/${SAMPLE_ID}/${FILE} && echo 0 || echo 1)
    if [ ${TEST} -eq 0 ]; then
        TEST=$(gsutil cp ${BUCKET_DIR}/pav/${SAMPLE_ID}/${FILE} ${FILE} && echo 0 || echo 1)
        if [ ${TEST} -eq 0 ]; then
            echo "1" >> ${IN_BUCKET_FILE}
        else 
            echo "0" >> ${IN_BUCKET_FILE}
        fi
    else
        echo "0" >> ${IN_BUCKET_FILE}
    fi
done < ${STATE_FILES_LIST}
cat ${IN_BUCKET_FILE}
