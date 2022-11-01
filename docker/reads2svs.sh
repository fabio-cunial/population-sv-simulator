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
# Time analysis on a 2-core node:
# minimap2 of one chunk: ~2 mins and 1.8 GB of RAM
# samtools sort of one chunk: ~3 s and 400 MB of RAM
# samtools merge: ~20 s and 10 MB of RAM
# pbsv discover: 
# pbsv call: 
# sniffles1: 
# sniffles2: 
# hifiasm:
# PAV: Max requirements from the original WDL: 8 threads, 32 GB of RAM, 1 GB of disk.
#
READS_FILE=$1
SAMPLE_ID=$2  # SM field in the .sam file (needed later for joint calling)
LENGTH=$3
MIN_COVERAGE=$4
MAX_COVERAGE=$5
COVERAGES=$6  # String, separated by "-".
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
WORK_DIR=${18}
DOCKER_DIR=${19}
KEEP_ASSEMBLIES=${20}

GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
TIME_COMMAND="/usr/bin/time --verbose"
COVERAGES=$(echo ${COVERAGES} | tr '-' ' ')
N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
MINIMAP_COMMAND="minimap2 -t ${N_THREADS} -aYx map-hifi --eqx"
READ_GROUP="@RG\tID:movie\tSM:${SAMPLE_ID}"
REFERENCE_LENGTH=$(wc -c < ${REFERENCE_FA})
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
N_ROWS_1X=$(( ( ${N_ROWS} / (4*${MAX_COVERAGE}) ) * 2 ))
rm -f chunk-*
split -d -l ${N_ROWS_1X} ${READS_FILE} chunk-
rm -f ${READS_FILE}
ALIGNMENTS_PREFIX="alignments_i${ID1}_i${ID2}_l${LENGTH}_c${MAX_COVERAGE}"
for CHUNK in $( find . -maxdepth 1 -name 'chunk-*' ); do
    FILE_NAME="${BUCKET_DIR}/alignments/${ALIGNMENTS_PREFIX}_$(basename ${CHUNK}).bam"
    TEST=$(gsutil -q stat ${FILE_NAME} && echo 0 || echo 1)
    if [ ${TEST} -eq 0 ]; then
        gsutil cp ${FILE_NAME} ${CHUNK}.bam
    else
        mv ${CHUNK} ${CHUNK}.fa
    	${TIME_COMMAND} ${MINIMAP_COMMAND} -R ${READ_GROUP} ${REFERENCE_MMI} ${CHUNK}.fa > ${CHUNK}.sam
        mv ${CHUNK}.fa ${CHUNK}
    	${TIME_COMMAND} samtools calmd -@ ${N_THREADS} -b ${CHUNK}.sam ${REFERENCE_FA} > ${CHUNK}.1.bam
    	rm -f ${CHUNK}.sam
    	${TIME_COMMAND} samtools sort -@ ${N_THREADS} ${CHUNK}.1.bam > ${CHUNK}.bam
    	rm -f ${CHUNK}.1.bam
        gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp ${CHUNK}.bam ${FILE_NAME}
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
    INFIX="i${ID1}_i${ID2}_l${LENGTH}_c${COVERAGE}"
    if [ ${USE_PBSV} -eq 1 ]; then
        PREFIX="pbsv_${INFIX}"
        TEST=$(gsutil -q stat ${BUCKET_DIR}/signatures/${PREFIX}.svsig.gz && echo 0 || echo 1)
        if [ ${TEST} -eq 0 ]; then
            ${TIME_COMMAND} gsutil cp ${BUCKET_DIR}/signatures/${PREFIX}.svsig.gz .
        else 
            # <discover> is sequential
        	${TIME_COMMAND} pbsv discover --tandem-repeats ${REFERENCE_TANDEM_REPEATS} coverage_${COVERAGE}.bam ${PREFIX}.svsig.gz
            gsutil cp ${PREFIX}.svsig.gz ${BUCKET_DIR}/signatures/
        fi
        TEST=$(gsutil -q stat ${BUCKET_DIR}/vcfs/${PREFIX}.vcf && echo 0 || echo 1)
        if [ ${TEST} -eq 1 ]; then
        	${TIME_COMMAND} pbsv call -j ${N_THREADS} --ccs ${REFERENCE_FA} ${PREFIX}.svsig.gz ${PREFIX}.vcf
            gsutil cp ${PREFIX}.vcf ${BUCKET_DIR}/vcfs/
        fi
        rm -f ${PREFIX}.*
    fi
	
	# SNIFFLES 1
    if [ ${USE_SNIFFLES1} -eq 1 ]; then 
        PREFIX="sniffles1_${INFIX}"
        TEST=$(gsutil -q stat ${BUCKET_DIR}/vcfs/${PREFIX}.vcf && echo 0 || echo 1)
        if [ ${TEST} -eq 1 ]; then
            ${TIME_COMMAND} sniffles1 -t ${N_THREADS} -m coverage_${COVERAGE}.bam -v ${PREFIX}.vcf
            gsutil cp ${PREFIX}.vcf ${BUCKET_DIR}/vcfs/
            rm -f ${PREFIX}.*
        fi
    fi
    
    # SNIFFLES 2
    if [ ${USE_SNIFFLES2} -eq 1 ]; then 
        PREFIX="sniffles2_${INFIX}"
        TEST=$(gsutil -q stat ${BUCKET_DIR}/signatures/${PREFIX}.vcf && echo 0 || echo 1)
        if [ ${TEST} -eq 1 ]; then
    	    ${TIME_COMMAND} sniffles --threads ${N_THREADS} --tandem-repeats ${REFERENCE_TANDEM_REPEATS} --reference ${REFERENCE_FA} --sample-id ${SAMPLE_ID} --input coverage_${COVERAGE}.bam --vcf ${PREFIX}.vcf --snf ${PREFIX}.snf
            gsutil cp ${PREFIX}.snf ${BUCKET_DIR}/signatures/
            gsutil cp ${PREFIX}.vcf ${BUCKET_DIR}/vcfs/
            rm -f ${PREFIX}.*
        fi
    fi
    
    # HIFIASM + QUAST
    PREFIX="assembly_${INFIX}"
    if [ ${USE_HIFIASM} -eq 1 -o ${USE_PAV} -eq 1 ]; then 
        TEST1=$(gsutil -q stat ${BUCKET_DIR}/assemblies/${PREFIX}_h1.fa && echo 0 || echo 1)
        TEST2=$(gsutil -q stat ${BUCKET_DIR}/assemblies/${PREFIX}_h2.fa && echo 0 || echo 1)
        if [ ${TEST1} -eq 1 -o ${TEST2} -eq 1 ]; then
            ${TIME_COMMAND} hifiasm -t ${N_THREADS} --hom-cov $(( ${COVERAGE}*2 )) -o tmpasm coverage_${COVERAGE}.fa
            awk '/^S/{print ">"$2; print $3}' tmpasm.bp.hap1.p_ctg.gfa > ${PREFIX}_h1.fa
            awk '/^S/{print ">"$2; print $3}' tmpasm.bp.hap2.p_ctg.gfa > ${PREFIX}_h2.fa
            ${TIME_COMMAND} gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp ${PREFIX}_h1.fa ${BUCKET_DIR}/assemblies/
            ${TIME_COMMAND} gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp ${PREFIX}_h2.fa ${BUCKET_DIR}/assemblies/
            awk '/^S/{print ">"$2; print $3}' tmpasm.bp.p_ctg.gfa > ${PREFIX}_primaryContigs.fa
            quast --threads ${N_THREADS} --eukaryote --large --est-ref-size ${REFERENCE_LENGTH} --no-gc ${PREFIX}_primaryContigs.fa
            rm -f ${PREFIX}_primaryContigs.fa
            tar -czf ${PREFIX}_quast.tar.gz quast_results/ 
            rm -rf quast_results/
            gsutil cp ${PREFIX}_quast.tar.gz ${BUCKET_DIR}/assemblies/
            rm -f ${PREFIX}_quast.tar.gz
            rm -rf tmpasm*
        elif [ ${USE_PAV} -eq 1 ]; then
            TEST=$(gsutil -q stat ${BUCKET_DIR}/vcfs/pav_i${ID1}_i${ID2}_l${LENGTH}_c${COVERAGE}.vcf && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                ${TIME_COMMAND} gsutil cp ${BUCKET_DIR}/assemblies/${PREFIX}_h1.fa ${PREFIX}_h1.fa
                ${TIME_COMMAND} gsutil cp ${BUCKET_DIR}/assemblies/${PREFIX}_h2.fa ${PREFIX}_h2.fa
            fi
        fi
    fi
    
    # PAV
    if [ ${USE_PAV} -eq 1 ]; then
        ASSEMBLY_PREFIX=${PREFIX}
        PREFIX="pav_${INFIX}"
        TEST=$(gsutil -q stat ${BUCKET_DIR}/vcfs/${PREFIX}.vcf && echo 0 || echo 1)
        if [ ${TEST} -eq 1 ]; then
            if [ ! -e config.json ]; then
                echo "{" >> config.json
                echo "\"reference\": \"asm/ref.fa\"," >> config.json
                echo "\"asm_pattern\": \"asm/{asm_name}/{hap}.fa.gz\"," >> config.json
                echo "\"aligner\": \"minimap2\"," >> config.json
                echo "\"map_threads\": ${N_THREADS}," >> config.json
                echo "\"merge_threads\": ${N_THREADS}," >> config.json
                echo "\"inv_threads\": ${N_THREADS}," >> config.json
                echo "\"inv_threads_lg\": ${N_THREADS}," >> config.json
                echo "\"lg_batch_count\": 1," >> config.json
                echo "\"inv_sig_batch_count\": 1" >> config.json
                echo "}" >> config.json
            fi
            mkdir -p asm/sample
            cp ${REFERENCE_FA} asm/ref.fa
            cp ${REFERENCE_FAI} asm/ref.fa.fai
            bgzip -@ ${N_THREADS} --compress-level 0 -c ${ASSEMBLY_PREFIX}_h1.fa > asm/sample/h1.fa.gz
            bgzip -@ ${N_THREADS} --compress-level 0 -c ${ASSEMBLY_PREFIX}_h2.fa > asm/sample/h2.fa.gz
            mkdir -p temp/sample/align/
            cp asm/sample/h1.fa.gz temp/sample/align/contigs_h1.fa.gz
            samtools faidx temp/sample/align/contigs_h1.fa.gz
            cp asm/sample/h2.fa.gz temp/sample/align/contigs_h2.fa.gz
            samtools faidx temp/sample/align/contigs_h2.fa.gz
            bash ${DOCKER_DIR}/pav_restoreCheckpoint.sh ${BUCKET_DIR} ${INFIX} ${WORK_DIR} ${DOCKER_DIR}
            bash ${DOCKER_DIR}/pav_checkpointDaemon.sh ${BUCKET_DIR} ${INFIX} ${WORK_DIR} ${DOCKER_DIR} &
            DAEMON_ID=$!
            source activate lr-pav
            ${TIME_COMMAND} snakemake -s ${DOCKER_DIR}/pav/Snakefile --cores ${N_THREADS} pav_sample.vcf.gz
            conda deactivate
            kill -KILL ${DAEMON_ID} || echo "The process ID <${DAEMON_ID}> of <pav_checkpointDaemon.sh> cannot be found."
            bcftools filter --exclude 'SVTYPE="SNV" || (SVLEN>-40 && SVLEN<40)' pav_sample.vcf.gz > ${PREFIX}.vcf
            rm -f pav_sample.vcf.gz*
            gsutil cp ${PREFIX}.vcf ${BUCKET_DIR}/vcfs/
            rm -f ${PREFIX}.vcf
            rm -rf asm/ data/ temp/ results/ log/
            gsutil -m rm -rf "${BUCKET_DIR}/pav/${INFIX}" || echo "Cannot remove directory ${BUCKET_DIR}/pav/${INFIX}"
            # Removing local assemblies as well
            rm -f ${ASSEMBLY_PREFIX}*
        fi
    fi
    
    tree -L 2
    df -h
    
    # Next iteration
    echo "${SAMPLE_ID} ${LENGTH} ${COVERAGE}" >> ${CHECKPOINT_FILE}
    gsutil cp ${CHECKPOINT_FILE} ${BUCKET_DIR}/checkpoints/
    PREVIOUS_COVERAGE=${COVERAGE}
done
rm -f coverage_* chunk-*

# Cleaning the bucket
gsutil rm -f ${BUCKET_DIR}/reads/${READS_FILE} || echo "Cannot remove file ${BUCKET_DIR}/reads/${READS_FILE}"
gsutil -m rm -f "${BUCKET_DIR}/alignments/${ALIGNMENTS_PREFIX}_chunk-*.bam" || echo "Cannot remove some of the files ${BUCKET_DIR}/alignments/${ALIGNMENTS_PREFIX}_chunk-*.bam"
if [ ${KEEP_ASSEMBLIES} -ne 1 ]; then
    gsutil -m rm -f "${BUCKET_DIR}/assemblies/assembly_i${ID1}_i${ID2}_l${LENGTH}_*" || echo "Cannot remove some of the files ${BUCKET_DIR}/assemblies/assembly_i${ID1}_i${ID2}_l${LENGTH}_*"
fi
