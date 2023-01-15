#!/bin/bash
#
# Splits the reads of a diploid individual into 1x diploid chunks and aligns
# each chunk to the reference. Then, for each coverage, merges a prefix of the
# chunks and runs the callers on the merged BAM.
#
# Resource analysis for 20x coverage of one chr1 haplotype. Intel Xeon,
# 2.30GHz, 8 physical cores.
#
# TASK                    | TIME   | RAM    | CORES | COMMENT
# minimap2 one chunk      | 30 s   | 3 GB   |   7   |
# samtools sort one chunk | 2 s    | 500 MB |   6   |
# samtools merge          | 30 s   | 15 MB  |   6   |
# pbsv discover           | 3 min  | 300 MB |   1   |
# pbsv call               | 2 min  | 3 GB   |   1   |
# sniffles1               | 4 min  | 100 MB |   4   |
# sniffles2               | 3 min  | 150 MB |   1   |
# hifiasm                 | 45 min | 20 GB  |   8   |
# PAV                     | 1..8 h | 7 GB   | 2...8 | faster with high covg
# paftools                | 2 s    | 300 MB |   1   |
#
READS_FILE=$1
SAMPLE_ID=$2  # SM field in the .sam file (needed later for joint calling)
LENGTH=$3
MIN_COVERAGE=$4  # Min value in $COVERAGES$.
MAX_COVERAGE=$5  # Max value in $COVERAGES$.
COVERAGES=$6  # Of each haplotype. String, separated by "-".
REFERENCE_FA=$7
REFERENCE_FAI=$8
REFERENCE_MMI=$9
REFERENCE_TANDEM_REPEATS=${10}
CHECKPOINT_FILE=${11}
BUCKET_DIR=${12}  # Root dir of the simulation in the bucket
USE_PBSV=${13}
USE_SNIFFLES1=${14}
USE_SNIFFLES2=${15}
USE_HIFIASM=${16}
USE_PAV=${17}
USE_PAFTOOLS=${18}
WORK_DIR=${19}
DOCKER_DIR=${20}
KEEP_ASSEMBLIES=${21}

GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
GSUTIL_DELAY_S="600"
TIME_COMMAND="/usr/bin/time --verbose"
COVERAGES=$(echo ${COVERAGES} | tr '-' ' ')
N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
MINIMAP_COMMAND="minimap2 -t ${N_THREADS} -aYx map-hifi --eqx"
MINIMAP2_ASM_FLAG="asm5"  # asm5/asm10/asm20 for ~0.1/1/5% sequence divergence
PAFTOOLS_MIN_SV_LENGTH="40"
READ_GROUP="@RG\tID:movie\tSM:${SAMPLE_ID}"
REFERENCE_LENGTH=$(cut -f 2 ${REFERENCE_FAI} | awk '{s+=$1} END {print s}')
ID1=${SAMPLE_ID}
ID2=$((${SAMPLE_ID} + 1))

set -euxo pipefail
echo "Running <reads2svs.sh> on ${N_THREADS} cores on the following node:"
lscpu
cat /proc/meminfo
cd ${WORK_DIR}

# Splitting the reads into chunks equal to 1x of a diploid individual, and 
# aligning each chunk to the reference in isolation.
N_ROWS=$(wc -l < ${READS_FILE})
N_ROWS_1X=$(( ${N_ROWS} / ${MAX_COVERAGE} ))
if [ $((${N_ROWS_1X} % 4)) -ne 0 ]; then
    # Making sure it is a multiple of 4, to make FASTQ files work.
    N_ROWS_1X=$(( (${N_ROWS_1X}/4 + 1)*4 ))
fi
rm -f chunk-*
split -d -l ${N_ROWS_1X} ${READS_FILE} chunk-
rm -f ${READS_FILE}
ALIGNMENTS_PREFIX="alignments_i${ID1}_i${ID2}_l${LENGTH}_c${MAX_COVERAGE}"
for CHUNK in $( find . -maxdepth 1 -name 'chunk-*' ); do
    FILE_NAME="${BUCKET_DIR}/alignments/${ALIGNMENTS_PREFIX}_$(basename ${CHUNK}).bam"
    TEST=$(gsutil -q stat ${FILE_NAME} && echo 0 || echo 1)
    if [ ${TEST} -eq 0 ]; then
        TEST=$(gsutil cp ${FILE_NAME} ${CHUNK}.bam && echo 0 || echo 1)
        if [ ${TEST} -eq 1 ]; then
            echo "Error downloading file <${FILE_NAME}>."
        fi
    fi
    if [ ${TEST} -eq 1 ]; then
        mv ${CHUNK} ${CHUNK}.fa
    	${TIME_COMMAND} ${MINIMAP_COMMAND} -R ${READ_GROUP} ${REFERENCE_MMI} ${CHUNK}.fa > ${CHUNK}.sam
        mv ${CHUNK}.fa ${CHUNK}
    	${TIME_COMMAND} samtools calmd -@ ${N_THREADS} -b ${CHUNK}.sam ${REFERENCE_FA} > ${CHUNK}.1.bam
    	rm -f ${CHUNK}.sam
    	${TIME_COMMAND} samtools sort -@ ${N_THREADS} ${CHUNK}.1.bam > ${CHUNK}.bam
    	rm -f ${CHUNK}.1.bam
        gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp ${CHUNK}.bam ${FILE_NAME} || echo "Error uploading file <${CHUNK}.bam> to bucket."
    fi
done

# Building the BAM and reads file of the smallest coverage
echo "Starting coverage ${MIN_COVERAGE}..."
rm -f coverage_${MIN_COVERAGE}.bam coverage_${MIN_COVERAGE}.fa
IDS=""
for i in $(seq -f "%02g" 0 $(( ${MIN_COVERAGE}-1 )) ); do
	IDS="${IDS} chunk-${i}.bam"
    cat chunk-${i} >> coverage_${MIN_COVERAGE}.fa
done
${TIME_COMMAND} samtools merge -@ ${N_THREADS} -o coverage_${MIN_COVERAGE}.bam ${IDS}
samtools index -@ ${N_THREADS} coverage_${MIN_COVERAGE}.bam

# Iterating over coverages
PREVIOUS_COVERAGE="0"
for COVERAGE in ${COVERAGES}; do
    echo "Starting coverage ${COVERAGE}..."
    CHECKPOINT_INDIVIDUAL=$(tail -n1 ${CHECKPOINT_FILE} | awk '{ print $1 }')
    CHECKPOINT_LENGTH=$(tail -n1 ${CHECKPOINT_FILE} | awk '{ print $2 }')
    CHECKPOINT_COVERAGE=$(tail -n1 ${CHECKPOINT_FILE} | awk '{ print $3 }')
    if [ ${PREVIOUS_COVERAGE} -ne 0 ]; then
		IDS="coverage_${PREVIOUS_COVERAGE}.bam"
		for i in $(seq -f "%02g" ${PREVIOUS_COVERAGE} $(( ${COVERAGE}-1 )) ); do
			IDS="${IDS} chunk-${i}.bam"
		done
		${TIME_COMMAND} samtools merge -@ ${N_THREADS} -o coverage_${COVERAGE}.bam ${IDS}
		samtools index -@ ${N_THREADS} coverage_${COVERAGE}.bam
        rm -f coverage_${PREVIOUS_COVERAGE}.bam coverage_${PREVIOUS_COVERAGE}.bai
		cp coverage_${PREVIOUS_COVERAGE}.fa coverage_${COVERAGE}.fa
		for i in $(seq -f "%02g" ${PREVIOUS_COVERAGE} $(( ${COVERAGE}-1 )) ); do
            cat chunk-${i} >> coverage_${COVERAGE}.fa
		done
        rm -f coverage_${PREVIOUS_COVERAGE}.fa
    fi
    if [ ${ID1} -eq ${CHECKPOINT_INDIVIDUAL} -a ${LENGTH} -eq ${CHECKPOINT_LENGTH} -a ${COVERAGE} -le ${CHECKPOINT_COVERAGE} ]; then
        PREVIOUS_COVERAGE=${COVERAGE}
        continue
    fi
    INFIX="i${ID1}_i${ID2}_l${LENGTH}_c${COVERAGE}"
    bash ${DOCKER_DIR}/reads2svs_impl.sh ${INFIX} coverage_${COVERAGE}.bam coverage_${COVERAGE}.fa ${COVERAGE} ${SAMPLE_ID} ${N_THREADS} ${REFERENCE_FA} ${REFERENCE_FAI} ${REFERENCE_TANDEM_REPEATS} ${BUCKET_DIR} ${USE_PBSV} ${USE_SNIFFLES1} ${USE_SNIFFLES2} ${USE_HIFIASM} ${USE_PAV} ${USE_PAFTOOLS} ${KEEP_ASSEMBLIES} ${WORK_DIR} ${DOCKER_DIR} 0
    
    # Next iteration
    echo "${SAMPLE_ID} ${LENGTH} ${COVERAGE}" >> ${CHECKPOINT_FILE}
    gsutil cp ${CHECKPOINT_FILE} ${BUCKET_DIR}/checkpoints/ || echo "Cannot upload checkpoint file to <${BUCKET_DIR}/checkpoints/>"
    PREVIOUS_COVERAGE=${COVERAGE}
    tree -L 2 || echo ""
    df -h
done
rm -f coverage_* chunk-*

# Cleaning the bucket
gsutil rm -f ${BUCKET_DIR}/reads/${READS_FILE} || echo "Error removing file ${BUCKET_DIR}/reads/${READS_FILE}"
gsutil -m rm -f "${BUCKET_DIR}/alignments/${ALIGNMENTS_PREFIX}_chunk-*.bam" || echo "Error removing some of the files ${BUCKET_DIR}/alignments/${ALIGNMENTS_PREFIX}_chunk-*.bam"
if [ ${KEEP_ASSEMBLIES} -ne 1 ]; then
    gsutil -m rm -f "${BUCKET_DIR}/assemblies/assembly_i${ID1}_i${ID2}_l${LENGTH}_*" || echo "Error removing some of the files ${BUCKET_DIR}/assemblies/assembly_i${ID1}_i${ID2}_l${LENGTH}_*"
fi
