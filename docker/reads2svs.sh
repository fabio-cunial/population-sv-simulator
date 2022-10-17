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
MIN_COVERAGE=$3
MAX_COVERAGE=$4
COVERAGES=$5  # String, separated by "-".
REFERENCE_FA=$6
REFERENCE_MMI=$7
REFERENCE_TANDEM_REPEATS=$8
CHECKPOINT_FILE=$9
BUCKET_ADDRESS=${10}  # Root dir of the simulation in the bucket
USE_PBSV=${11}
USE_SNIFFLES1=${12}
USE_SNIFFLES2=${13}
USE_HIFIASM=${14}

TIME_COMMAND="/usr/bin/time --verbose"
COVERAGES=$(echo ${COVERAGES} | tr '-' ' ')
N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
MINIMAP_COMMAND="minimap2 -t ${N_THREADS} -aYx map-hifi --eqx"
READ_GROUP="@RG\tID:movie\tSM:${SAMPLE_ID}"

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
for CHUNK in $(find . -maxdepth 1 "chunk-*"); do
	mv ${CHUNK} ${CHUNK}.fa
done
ALIGNMENTS_FILE="alignments_i${ID1}_i${ID2}_l${LENGTH}_c${MAX_COVERAGE}.tar"
gsutil -q stat ${BUCKET_ADDRESS}/alignments/${ALIGNMENTS_FILE}
if [ $? -eq 0 ]; then
    ${TIME_COMMAND} gsutil cp ${BUCKET_ADDRESS}/alignments/${ALIGNMENTS_FILE} .
    tar -xf ${ALIGNMENTS_FILE}
    rm -f ${ALIGNMENTS_FILE}
else
    for CHUNK in $(find . -maxdepth 1 "chunk-*"); do
    	${TIME_COMMAND} ${MINIMAP_COMMAND} -R ${READ_GROUP} ${REFERENCE_MMI} ${CHUNK}.fa > ${CHUNK}.sam
    	${TIME_COMMAND} samtools calmd -@ ${N_THREADS} -b ${CHUNK}.sam ${REFERENCE_FA} > ${CHUNK}.1.bam
    	rm -f ${CHUNK}.sam
    	${TIME_COMMAND} samtools sort -@ ${N_THREADS} ${CHUNK}.1.bam > ${CHUNK}.bam
    	rm -f ${CHUNK}.1.bam
    done
    tar -cf ${ALIGNMENTS_FILE} chunk-*.bam
    ${TIME_COMMAND} gsutil cp ${ALIGNMENTS_FILE} ${BUCKET_ADDRESS}/alignments/
    rm -f ${ALIGNMENTS_FILE}
done

# Building the BAM and reads file of the smallest coverage
echo "Starting coverage ${MIN_COVERAGE}..."
rm -f coverage_${MIN_COVERAGE}.bam coverage_${MIN_COVERAGE}.fa
IDS=""
for i in $(seq -f "%02g" 0 $(( ${MIN_COVERAGE}-1 )) ); do
	IDS="${IDS} chunk-${i}.bam"
    cat chunk-${i}.fa >> coverage_${MIN_COVERAGE}.fa
done
${TIME_COMMAND} samtools merge -@ ${N_THREADS} -o coverage_${MIN_COVERAGE}.bam ${IDS}
samtools index -@ ${N_THREADS} coverage_${MIN_COVERAGE}.bam

# Iterating over coverages
PREVIOUS_COVERAGE="0"
for COVERAGE in ${COVERAGES}; do
    echo "Starting coverage ${COVERAGE}..."
    CHECKPOINT_COVERAGE=$(tail -n1 ${CHECKPOINT_FILE} | awk '{ print $3 }')
    if [ ${COVERAGE} -lt ${CHECKPOINT_COVERAGE} ]; then
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
            cat chunk-${i}.fa >> coverage_${COVERAGE}.fa
		done
        rm -f coverage_${PREVIOUS_COVERAGE}.fa
    fi
    
	# PBSV
    if [ ${USE_PBSV} -eq 1 ]; then
        PREFIX="pbsv_i${SAMPLE_ID}_i$((${SAMPLE_ID} + 1))_l${LENGTH}_c${COVERAGE}"
        gsutil -q stat ${BUCKET_ADDRESS}/signatures/${PREFIX}.svsig.gz
        if [ $? -eq 0 ]; then
            ${TIME_COMMAND} gsutil cp ${BUCKET_ADDRESS}/signatures/${PREFIX}.svsig.gz .
        else 
            # <discover> is sequential
        	${TIME_COMMAND} pbsv discover --tandem-repeats ${REFERENCE_TANDEM_REPEATS} coverage_${COVERAGE}.bam ${PREFIX}.svsig.gz
            gsutil cp ${PREFIX}.svsig.gz ${BUCKET_ADDRESS}/signatures/
        fi
        gsutil -q stat ${BUCKET_ADDRESS}/vcfs/${PREFIX}.vcf
        if [ ! $? -eq 0 ]; then
        	${TIME_COMMAND} pbsv call -j ${N_THREADS} --ccs ${REFERENCE_FA} ${PREFIX}.svsig.gz ${PREFIX}.vcf
            gsutil cp ${PREFIX}.vcf ${BUCKET_ADDRESS}/vcfs/
        fi
        rm -f ${PREFIX}.*
    fi
	
	# SNIFFLES 1
    if [ ${USE_SNIFFLES1} -eq 1 ]; then 
        PREFIX="sniffles1_i${SAMPLE_ID}_i$((${SAMPLE_ID} + 1))_l${LENGTH}_c${COVERAGE}"
        gsutil -q stat ${BUCKET_ADDRESS}/vcfs/${PREFIX}.vcf
        if [ ! $? -eq 0 ]; then
            ${TIME_COMMAND} sniffles1 -t ${N_THREADS} -m coverage_${COVERAGE}.bam -v ${PREFIX}.vcf
            gsutil cp ${PREFIX}.vcf ${BUCKET_ADDRESS}/vcfs/
            rm -f ${PREFIX}.*
        fi
    fi
    
    # SNIFFLES 2
    if [ ${USE_SNIFFLES2} -eq 1 ]; then 
        PREFIX="sniffles2_i${SAMPLE_ID}_i$((${SAMPLE_ID} + 1))_l${LENGTH}_c${COVERAGE}"
        gsutil -q stat ${BUCKET_ADDRESS}/signatures/${PREFIX}.vcf
        if [ ! $? -eq 0 ]; then
    	    ${TIME_COMMAND} sniffles --threads ${N_THREADS} --tandem-repeats ${REFERENCE_TANDEM_REPEATS} --reference ${REFERENCE_FA} --sample-id ${SAMPLE_ID} --input coverage_${COVERAGE}.bam --vcf ${PREFIX}.vcf --snf ${PREFIX}.snf
            gsutil cp ${PREFIX}.snf ${BUCKET_ADDRESS}/signatures/
            gsutil cp ${PREFIX}.vcf ${BUCKET_ADDRESS}/vcfs/
            rm -f ${PREFIX}.*
        fi
    fi
    
    # HIFIASM
    if [ ${USE_HIFIASM} -eq 1 ]; then 
        PREFIX="assembly_i${SAMPLE_ID}_i$((${SAMPLE_ID} + 1))_l${LENGTH}_c${COVERAGE}"
        gsutil -q stat ${BUCKET_ADDRESS}/assemblies/${PREFIX}.tar
        if [ ! $? -eq 0 ]; then
            ${TIME_COMMAND} hifiasm --hom-cov $(( ${COVERAGE}*2 )) -o tmpasm -t ${N_THREADS} coverage_${COVERAGE}.fa
            awk '/^S/{print ">"$2; print $3}' tmpasm.bp.hap1.p_ctg.gfa > ${PREFIX}_h1.fa
            awk '/^S/{print ">"$2; print $3}' tmpasm.bp.hap2.p_ctg.gfa > ${PREFIX}_h2.fa
            rm -rf tmpasm*
            tar -cf ${PREFIX}.tar ${PREFIX}_*.fa
            rm -rf ${PREFIX}_*.fa
            ${TIME_COMMAND} gsutil cp ${PREFIX}.tar ${BUCKET_ADDRESS}/assemblies/
            rm -f ${PREFIX}.*
        fi
    fi
    
    # Next iteration
    echo "${SAMPLE_ID} ${LENGTH} ${COVERAGE}" >> ${CHECKPOINT_FILE}
    gsutil cp ${CHECKPOINT_FILE} ${BUCKET_ADDRESS}/checkpoints/
    PREVIOUS_COVERAGE=${COVERAGE}
done
rm -f coverage_* chunk-*
gsutil rm -f ~{bucket_address}/alignments/${ALIGNMENTS_FILE}
