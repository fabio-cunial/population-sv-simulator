#!/bin/bash
#
# Given a read set for the right mode at a fixed coverage, and a read set for
# the left mode at its max coverage, the program builds several merged read sets
# that all include the full right dataset but have different values of left
# coverage. For every such dataset, the program runs several callers.
#
READS_FILE_LEFT=$1
READS_FILE_RIGHT=$2  # -1=no such file exists.
SAMPLE_ID=$3  # SM field in the .sam file (needed later for joint calling)
MIN_COVERAGE_LEFT=$4  # Min value in $COVERAGES_LEFT$.
MAX_COVERAGE_LEFT=$5    # Max value in $COVERAGES_LEFT$. It is assumed to be equal to the total coverage of $READS_FILE_LEFT$.
COVERAGES_LEFT=$6  # Of each haplotype. Sorted in increasing order. Assumed to be floats and multiples of $COVERAGE_QUANTUM$. String separated by "-".
REFERENCE_FA=$7
REFERENCE_FAI=$8
REFERENCE_TANDEM_REPEATS=$9
BUCKET_DIR=${10}  # Root dir of the simulation in the bucket
USE_PBSV=${11}
USE_SNIFFLES1=${12}
USE_SNIFFLES2=${13}
USE_HIFIASM=${14}
USE_PAV=${15}
USE_PAFTOOLS=${16}
KEEP_ASSEMBLIES=${17}
WORK_DIR=${18}
DOCKER_DIR=${19}
COVERAGE_QUANTUM=${20}  # Of each haplotype. Assumed to be a double <=1.0.
BIN_LENGTH=${21}
MAX_READ_LENGTH=${22}
GENOME_LENGTH_HAPLOID=${23}

GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
GSUTIL_DELAY_S="600"
TIME_COMMAND="/usr/bin/time --verbose"
COVERAGES_LEFT=$(echo ${COVERAGES_LEFT} | tr '-' ' ')
N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
MINIMAP_COMMAND="minimap2 -t ${N_THREADS} -aYx map-hifi --eqx"
READ_GROUP="@RG\tID:movie\tSM:${SAMPLE_ID}"
export LC_NUMERIC="en_US.UTF-8"

set -euxo pipefail
cd ${WORK_DIR}

# Aligning the right reads to the reference (if any).
if [ ${READS_FILE_RIGHT} != "-1" ]; then
    TEST=$(gsutil -q stat ${BUCKET_DIR}/reads_question2_chunks/right.bam && echo 0 || echo 1)
    if [ ${TEST} -eq 0 ]; then
        while : ; do
            TEST=$(gsutil cp ${BUCKET_DIR}/reads_question2_chunks/right.bam . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading <${BUCKET_DIR}/reads_question2_chunks/right.bam>. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
    else
        ${TIME_COMMAND} ${MINIMAP_COMMAND} -R ${READ_GROUP} ${REFERENCE_FA} ${READS_FILE_RIGHT} > right.sam
        ${TIME_COMMAND} samtools sort -@ ${N_THREADS} --output-fmt BAM right.sam > right.1.bam
        rm -f right.sam
        ${TIME_COMMAND} samtools calmd -@ ${N_THREADS} -b right.1.bam ${REFERENCE_FA} > right.bam
        rm -f right.1.bam
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp right.bam ${BUCKET_DIR}/reads_question2_chunks/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading <${BUCKET_DIR}/reads_question2_chunks/right.bam>. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
    fi
fi

# Splitting the left reads into chunks equal to a fractional quantum of
# coverage of each haplotype, and aligning each chunk to the reference in
# isolation.
N_CHUNKS=$(echo "scale=0; ${MAX_COVERAGE_LEFT} / ${COVERAGE_QUANTUM}" | bc)
N_ROWS=$(wc -l < ${READS_FILE_LEFT})
N_ROWS_ONE_QUANTUM=$(( ${N_ROWS} / ${N_CHUNKS} ))
if [ $((${N_ROWS_ONE_QUANTUM} % 4)) -ne 0 ]; then
    # Making sure it is a multiple of 4, to make FASTQ files work.
    N_ROWS_ONE_QUANTUM=$(( (${N_ROWS_ONE_QUANTUM}/4 + 1)*4 ))
fi
rm -f chunk-*
split -d -l ${N_ROWS_ONE_QUANTUM} ${READS_FILE_LEFT} chunk-
rm -f ${READS_FILE_LEFT}
for CHUNK in $( find . -maxdepth 1 -name 'chunk-*' ); do
    FILE_NAME="${BUCKET_DIR}/reads_question2_chunks/$(basename ${CHUNK}).bam"
    TEST=$(gsutil -q stat ${FILE_NAME} && echo 0 || echo 1)
    if [ ${TEST} -eq 0 ]; then
        while : ; do
            TEST=$(gsutil cp ${FILE_NAME} . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading <${FILE_NAME}>. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
    else
        mv ${CHUNK} ${CHUNK}.fastq
        ${TIME_COMMAND} ${MINIMAP_COMMAND} -R ${READ_GROUP} ${REFERENCE_FA} ${CHUNK}.fastq > ${CHUNK}.sam
        mv ${CHUNK}.fastq ${CHUNK}
        tail -n 16 ${CHUNK}
        tail -n 16 ${CHUNK}.sam
        ${TIME_COMMAND} samtools sort -@ ${N_THREADS} --output-fmt BAM ${CHUNK}.sam > ${CHUNK}.1.bam
        rm -f ${CHUNK}.sam
        ${TIME_COMMAND} samtools calmd -@ ${N_THREADS} -b ${CHUNK}.1.bam ${REFERENCE_FA} > ${CHUNK}.bam
        rm -f ${CHUNK}.1.bam
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp ${CHUNK}.bam ${FILE_NAME} && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading <${FILE_NAME}>. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
    fi
done

# Building the BAM and reads file of the smallest coverage
echo "Starting coverage ${MIN_COVERAGE_LEFT} of each haplotype..."
rm -f coverage_${MIN_COVERAGE_LEFT}.bam coverage_${MIN_COVERAGE_LEFT}.fastq
mv ${READS_FILE_RIGHT} coverage_${MIN_COVERAGE_LEFT}.fastq
if [ ${MIN_COVERAGE_LEFT} == "0" -o ${MIN_COVERAGE_LEFT} == "0.0" -o ${MIN_COVERAGE_LEFT} == ".0" ]; then
    if [ ${READS_FILE_RIGHT} != "-1" ]; then
        mv right.bam coverage_${MIN_COVERAGE_LEFT}.bam
    fi
else
    LAST_CHUNK=$(echo "scale=0; ${MIN_COVERAGE_LEFT} / ${COVERAGE_QUANTUM}" | bc)
    LAST_CHUNK=$(( ${LAST_CHUNK}-1 ))
    IDS=""
    for i in $(seq -f "%02g" 0 ${LAST_CHUNK}); do
        IDS="${IDS} chunk-${i}.bam"
        cat chunk-${i} >> coverage_${MIN_COVERAGE_LEFT}.fastq
    done
    if [ ${READS_FILE_RIGHT} != "-1" ]; then
        ${TIME_COMMAND} samtools merge -@ ${N_THREADS} -o coverage_${MIN_COVERAGE_LEFT}.bam right.bam ${IDS}
        rm -f right.bam
    else
        ${TIME_COMMAND} samtools merge -@ ${N_THREADS} -o coverage_${MIN_COVERAGE_LEFT}.bam ${IDS}
    fi
fi
samtools index -@ ${N_THREADS} coverage_${MIN_COVERAGE_LEFT}.bam

# Iterating over coverages
PREVIOUS_COVERAGE="-1"; PREVIOUS_LAST_CHUNK="-1";
for COVERAGE in ${COVERAGES_LEFT}; do
    echo "Starting coverage ${COVERAGE} of each haplotype..."
    TEST1=$(gsutil -q stat ${BUCKET_DIR}/reads_c${COVERAGE}/coverage_${COVERAGE}.fastq && echo 0 || echo 1)
    TEST2=$(gsutil -q stat ${BUCKET_DIR}/reads_c${COVERAGE}/coverage_${COVERAGE}.bam && echo 0 || echo 1)
    if [ ${TEST1} -eq 0 -a ${TEST2} -eq 0 ]; then
        while : ; do
            TEST=$(gsutil cp ${BUCKET_DIR}/reads_c${COVERAGE}/coverage_${COVERAGE}.fastq . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading <${BUCKET_DIR}/reads_c${COVERAGE}/coverage_${COVERAGE}.fastq>. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        while : ; do
            TEST=$(gsutil cp ${BUCKET_DIR}/reads_c${COVERAGE}/coverage_${COVERAGE}.bam* . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading <${BUCKET_DIR}/reads_c${COVERAGE}/coverage_${COVERAGE}.bam>. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
    else
        if [ ${PREVIOUS_COVERAGE} != -1 ]; then
            IDS="coverage_${PREVIOUS_COVERAGE}.bam"
            LAST_CHUNK=$(echo "scale=0; ${COVERAGE} / ${COVERAGE_QUANTUM}" | bc)
            LAST_CHUNK=$(( ${LAST_CHUNK}-1 ))
            for i in $(seq -f "%02g" $((${PREVIOUS_LAST_CHUNK} + 1)) ${LAST_CHUNK} ); do
                IDS="${IDS} chunk-${i}.bam"
            done
            ${TIME_COMMAND} samtools merge -@ ${N_THREADS} -o coverage_${COVERAGE}.bam ${IDS}
            samtools index -@ ${N_THREADS} coverage_${COVERAGE}.bam
            rm -f coverage_${PREVIOUS_COVERAGE}.bam coverage_${PREVIOUS_COVERAGE}.bai
            cp coverage_${PREVIOUS_COVERAGE}.fastq coverage_${COVERAGE}.fastq
            for i in $(seq -f "%02g" $((${PREVIOUS_LAST_CHUNK} + 1)) ${LAST_CHUNK}); do
                cat chunk-${i} >> coverage_${COVERAGE}.fastq
            done
        fi
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp coverage_${COVERAGE}.fastq ${BUCKET_DIR}/reads_c${COVERAGE}/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading <coverage_${COVERAGE}.fastq>. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp coverage_${COVERAGE}.bam coverage_${COVERAGE}.bam.bai ${BUCKET_DIR}/reads_c${COVERAGE}/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading <coverage_${COVERAGE}.bam>. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
    fi
    java -cp ${DOCKER_DIR} Fastq2LengthHistogram coverage_${COVERAGE}.fastq ${BIN_LENGTH} ${MAX_READ_LENGTH} coverage_${COVERAGE}.fastq.histogram coverage_${COVERAGE}.fastq.max
    COVERAGE_EACH_HAPLOTYPE=$( sed -n '2~4p' coverage_${COVERAGE}.fastq | wc -c )
    echo "scale=8; ${COVERAGE_EACH_HAPLOTYPE} / (2.0*${GENOME_LENGTH_HAPLOID})" | bc > coverage_${COVERAGE}.fastq.coverage
    while : ; do
        TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp "coverage_${COVERAGE}.fastq.*" ${BUCKET_DIR}/reads_c${COVERAGE}/ && echo 0 || echo 1)
        if [ ${TEST} -eq 1 ]; then
            echo "Error uploading <coverage_${COVERAGE}.fastq.histogram>. Trying again..."
            sleep ${GSUTIL_DELAY_S}
        else
            break
        fi
    done
    bash ${DOCKER_DIR}/reads2svs_impl.sh ${SAMPLE_ID} coverage_${COVERAGE}.bam coverage_${COVERAGE}.fastq ${COVERAGE} ${SAMPLE_ID} ${N_THREADS} ${REFERENCE_FA} ${REFERENCE_FAI} ${REFERENCE_TANDEM_REPEATS} "${BUCKET_DIR}/reads_c${COVERAGE}" ${USE_PBSV} ${USE_SNIFFLES1} ${USE_SNIFFLES2} ${USE_HIFIASM} ${USE_PAV} ${USE_PAFTOOLS} ${KEEP_ASSEMBLIES} ${WORK_DIR} ${DOCKER_DIR} 1
    rm -f coverage_${PREVIOUS_COVERAGE}.bam coverage_${PREVIOUS_COVERAGE}.fastq
    
    # Next iteration
    PREVIOUS_COVERAGE=${COVERAGE}
    PREVIOUS_LAST_CHUNK=$(echo "scale=0; ${COVERAGE} / ${COVERAGE_QUANTUM}" | bc)
    PREVIOUS_LAST_CHUNK=$(( ${PREVIOUS_LAST_CHUNK}-1 ))
    tree -L 2 || echo ""
    df -h
done
rm -f coverage_* chunk-*
