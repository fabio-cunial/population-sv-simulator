#!/bin/bash
#
# Splits the reads of a diploid individual into 1x chunks and aligns each chunk
# to the reference. Then, for each coverage, merges a prefix of the chunks and
# runs the callers on the merged BAM.
#
# Remark: for chr1 at 30x coverage, the procedure adds approx. 22 GB to the
# current disk usage. Peak RAM with one thread and no hifiasm: 2 GB (from
# minimap2).
#
READS_FILE=$1
SAMPLE_ID=$2  # SM field in the .sam file (needed later for joint calling)
LENGTH=$3
MIN_COVERAGE=$4
MAX_COVERAGE=$5
COVERAGES=$6  # String, separated by "-".
REFERENCE_FA=$7
REFERENCE_MMI=$8
REFERENCE_TANDEM_REPEATS=$9
CHECKPOINT_FILE=${10}
BUCKET_DIR=${11}  # Root dir of the simulation in the bucket
USE_PBSV=${12}
USE_SNIFFLES1=${13}
USE_SNIFFLES2=${14}
USE_HIFIASM=${15}

GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
TIME_COMMAND="/usr/bin/time --verbose"
COVERAGES=$(echo ${COVERAGES} | tr '-' ' ')
N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
MINIMAP_COMMAND="minimap2 -t ${N_THREADS} -aYx map-hifi --eqx"
READ_GROUP="@RG\tID:movie\tSM:${SAMPLE_ID}"
ID1=${SAMPLE_ID}
ID2=$((${SAMPLE_ID} + 1))

set -euxo pipefail
echo "Running <reads2svs.sh> on ${N_THREADS} cores on the following node:"
lscpu
cat /proc/meminfo

# Splitting the reads into chunks equal to 1x of a diploid individual, and 
# aligning each chunk to the reference in isolation.
N_ROWS=$(wc -l < ${READS_FILE})
N_ROWS_1X=$(( ( ${N_ROWS} / (4*${MAX_COVERAGE}) ) * 2 ))
rm -f chunk-*
split -d -l ${N_ROWS_1X} ${READS_FILE} chunk-
rm -f ${READS_FILE}
ALIGNMENTS_PREFIX="alignments_i${ID1}_i${ID2}_l${LENGTH}_c${MAX_COVERAGE}"
for CHUNK in $( find . -maxdepth 1 -name 'chunk-*' ); do
    TEST=$(gsutil -q stat ${BUCKET_DIR}/alignments/${ALIGNMENTS_PREFIX}_${CHUNK}.bam || echo 1)
    if [ ${TEST} != 1 ]; then
        gsutil cp ${BUCKET_DIR}/alignments/${ALIGNMENTS_PREFIX}_${CHUNK}.bam ${CHUNK}.bam
    else
        mv ${CHUNK} ${CHUNK}.fa
    	${TIME_COMMAND} ${MINIMAP_COMMAND} -R ${READ_GROUP} ${REFERENCE_MMI} ${CHUNK}.fa > ${CHUNK}.sam
        mv ${CHUNK}.fa ${CHUNK}
    	${TIME_COMMAND} samtools calmd -@ ${N_THREADS} -b ${CHUNK}.sam ${REFERENCE_FA} > ${CHUNK}.1.bam
    	rm -f ${CHUNK}.sam
    	${TIME_COMMAND} samtools sort -@ ${N_THREADS} ${CHUNK}.1.bam > ${CHUNK}.bam
    	rm -f ${CHUNK}.1.bam
        gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp ${CHUNK}.bam ${BUCKET_DIR}/alignments/${ALIGNMENTS_PREFIX}_${CHUNK}.bam
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
    if [ ${ID1} -eq ${CHECKPOINT_INDIVIDUAL} -a ${LENGTH} -eq ${CHECKPOINT_LENGTH} -a ${COVERAGE} -le ${CHECKPOINT_COVERAGE} ]; then
        continue
    fi
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
    
	# PBSV
    if [ ${USE_PBSV} -eq 1 ]; then
        PREFIX="pbsv_i${ID1}_i${ID2}_l${LENGTH}_c${COVERAGE}"
        TEST=$(gsutil -q stat ${BUCKET_DIR}/signatures/${PREFIX}.svsig.gz || echo 1)
        if [ ${TEST} != 1 ]; then
            ${TIME_COMMAND} gsutil cp ${BUCKET_DIR}/signatures/${PREFIX}.svsig.gz .
        else 
            # <discover> is sequential
        	${TIME_COMMAND} pbsv discover --tandem-repeats ${REFERENCE_TANDEM_REPEATS} coverage_${COVERAGE}.bam ${PREFIX}.svsig.gz
            gsutil cp ${PREFIX}.svsig.gz ${BUCKET_DIR}/signatures/
        fi
        TEST=$(gsutil -q stat ${BUCKET_DIR}/vcfs/${PREFIX}.vcf || echo 1)
        if [ ${TEST} = 1 ]; then
        	${TIME_COMMAND} pbsv call -j ${N_THREADS} --ccs ${REFERENCE_FA} ${PREFIX}.svsig.gz ${PREFIX}.vcf
            gsutil cp ${PREFIX}.vcf ${BUCKET_DIR}/vcfs/
        fi
        rm -f ${PREFIX}.*
    fi
	
	# SNIFFLES 1
    if [ ${USE_SNIFFLES1} -eq 1 ]; then 
        PREFIX="sniffles1_i${ID1}_i${ID2}_l${LENGTH}_c${COVERAGE}"
        TEST=$(gsutil -q stat ${BUCKET_DIR}/vcfs/${PREFIX}.vcf || echo 1)
        if [ ${TEST} = 1 ]; then
            ${TIME_COMMAND} sniffles1 -t ${N_THREADS} -m coverage_${COVERAGE}.bam -v ${PREFIX}.vcf
            gsutil cp ${PREFIX}.vcf ${BUCKET_DIR}/vcfs/
            rm -f ${PREFIX}.*
        fi
    fi
    
    # SNIFFLES 2
    if [ ${USE_SNIFFLES2} -eq 1 ]; then 
        PREFIX="sniffles2_i${ID1}_i${ID2}_l${LENGTH}_c${COVERAGE}"
        TEST=$(gsutil -q stat ${BUCKET_DIR}/signatures/${PREFIX}.vcf || echo 1)
        if [ ${TEST} = 1 ]; then
    	    ${TIME_COMMAND} sniffles --threads ${N_THREADS} --tandem-repeats ${REFERENCE_TANDEM_REPEATS} --reference ${REFERENCE_FA} --sample-id ${SAMPLE_ID} --input coverage_${COVERAGE}.bam --vcf ${PREFIX}.vcf --snf ${PREFIX}.snf
            gsutil cp ${PREFIX}.snf ${BUCKET_DIR}/signatures/
            gsutil cp ${PREFIX}.vcf ${BUCKET_DIR}/vcfs/
            rm -f ${PREFIX}.*
        fi
    fi
    
    # HIFIASM
    if [ ${USE_HIFIASM} -eq 1 ]; then 
        PREFIX="assembly_i${ID1}_i${ID2}_l${LENGTH}_c${COVERAGE}"
        TEST1=$(gsutil -q stat ${BUCKET_DIR}/assemblies/${PREFIX}_h1.fa || echo 1)
        TEST2=$(gsutil -q stat ${BUCKET_DIR}/assemblies/${PREFIX}_h2.fa || echo 1)
        if [ ${TEST1} = 1 -o ${TEST2} = 1 ]; then
            ${TIME_COMMAND} hifiasm --hom-cov $(( ${COVERAGE}*2 )) -o tmpasm -t ${N_THREADS} coverage_${COVERAGE}.fa
            awk '/^S/{print ">"$2; print $3}' tmpasm.bp.hap1.p_ctg.gfa > ${PREFIX}_h1.fa
            awk '/^S/{print ">"$2; print $3}' tmpasm.bp.hap2.p_ctg.gfa > ${PREFIX}_h2.fa
            rm -rf tmpasm*            
            ${TIME_COMMAND} gsutil cp ${PREFIX}_h1.fa ${BUCKET_DIR}/assemblies/
            ${TIME_COMMAND} gsutil cp ${PREFIX}_h2.fa ${BUCKET_DIR}/assemblies/
            rm -f ${PREFIX}*
        fi
    fi
    
    # Next iteration
    echo "${SAMPLE_ID} ${LENGTH} ${COVERAGE}" >> ${CHECKPOINT_FILE}
    gsutil cp ${CHECKPOINT_FILE} ${BUCKET_DIR}/checkpoints/
    PREVIOUS_COVERAGE=${COVERAGE}
done
rm -f coverage_* chunk-*

# Cleaning the bucket
gsutil rm -f ${BUCKET_DIR}/reads/${READS_FILE}
gsutil rm -f "${BUCKET_DIR}/alignments/${ALIGNMENTS_PREFIX}_chunk-*.bam"
