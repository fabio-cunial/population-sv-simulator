#!/bin/bash
#
# Splits the reads of a diploid individual into 1x chunks and aligns each chunk
# to the reference. Then, for each coverage, merges a prefix of the chunks and
# runs the callers on the merged BAM.
#
# Resource analysis for 20x coverage of one haplotype. Intel Xeon, 2.30GHz, 8
# physical cores.
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
GSUTIL_DELAY_S="600"
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
    
	# PBSV
    PREFIX_PBSV="pbsv_${INFIX}"
    if [ ${USE_PBSV} -eq 1 ]; then
        TEST=$(gsutil -q stat ${BUCKET_DIR}/signatures/${PREFIX_PBSV}.svsig.gz && echo 0 || echo 1)
        if [ ${TEST} -eq 0 ]; then
            while : ; do
                TEST=$(gsutil cp ${BUCKET_DIR}/signatures/${PREFIX_PBSV}.svsig.gz . && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading file <${BUCKET_DIR}/signatures/${PREFIX_PBSV}.svsig.gz>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        else
            # <discover> is sequential
        	${TIME_COMMAND} pbsv discover --tandem-repeats ${REFERENCE_TANDEM_REPEATS} coverage_${COVERAGE}.bam ${PREFIX_PBSV}.svsig.gz
            while : ; do
                TEST=$(gsutil cp ${PREFIX_PBSV}.svsig.gz ${BUCKET_DIR}/signatures/ && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error uploading file <${PREFIX_PBSV}.svsig.gz>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            fi
        fi
        TEST=$(gsutil -q stat ${BUCKET_DIR}/vcfs/${PREFIX_PBSV}.vcf && echo 0 || echo 1)
        if [ ${TEST} -eq 1 ]; then
        	${TIME_COMMAND} pbsv call -j ${N_THREADS} --ccs ${REFERENCE_FA} ${PREFIX_PBSV}.svsig.gz ${PREFIX_PBSV}.vcf
            while : ; do
                TEST=$(gsutil cp ${PREFIX_PBSV}.vcf ${BUCKET_DIR}/vcfs/ && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error uploading file <${PREFIX_PBSV}.vcf>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        fi
    fi
    rm -f ${PREFIX_PBSV}*
	
	# SNIFFLES 1
    PREFIX_SNIFFLES1="sniffles1_${INFIX}"
    if [ ${USE_SNIFFLES1} -eq 1 ]; then 
        TEST=$(gsutil -q stat ${BUCKET_DIR}/vcfs/${PREFIX_SNIFFLES1}.vcf && echo 0 || echo 1)
        if [ ${TEST} -eq 1 ]; then
            ${TIME_COMMAND} sniffles1 -t ${N_THREADS} -m coverage_${COVERAGE}.bam -v ${PREFIX_SNIFFLES1}.vcf
            while : ; do
                TEST=$(gsutil cp ${PREFIX_SNIFFLES1}.vcf ${BUCKET_DIR}/vcfs/ && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error uploading file <${PREFIX_SNIFFLES1}.vcf>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        fi
    fi
    rm -f ${PREFIX_SNIFFLES1}*
    
    # SNIFFLES 2
    PREFIX_SNIFFLES2="sniffles2_${INFIX}"
    if [ ${USE_SNIFFLES2} -eq 1 ]; then 
        TEST=$(gsutil -q stat ${BUCKET_DIR}/signatures/${PREFIX_SNIFFLES2}.vcf && echo 0 || echo 1)
        if [ ${TEST} -eq 1 ]; then
    	    ${TIME_COMMAND} sniffles --threads ${N_THREADS} --tandem-repeats ${REFERENCE_TANDEM_REPEATS} --reference ${REFERENCE_FA} --sample-id ${SAMPLE_ID} --input coverage_${COVERAGE}.bam --vcf ${PREFIX_SNIFFLES2}.vcf --snf ${PREFIX_SNIFFLES2}.snf
            while : ; do
                TEST=$(gsutil cp ${PREFIX_SNIFFLES2}.snf ${BUCKET_DIR}/signatures/ && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error uploading file <${PREFIX_SNIFFLES2}.snf>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else 
                    break
                fi
            done
            while : ; do
                TEST=$(gsutil cp ${PREFIX_SNIFFLES2}.vcf ${BUCKET_DIR}/vcfs/ && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error uploading file <${PREFIX_SNIFFLES2}.vcf>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        fi
    fi
    rm -f ${PREFIX_SNIFFLES2}*
    
    # HIFIASM + QUAST
    PREFIX_ASSEMBLY="assembly_${INFIX}"
    PREFIX_PAV="pav_${INFIX}"
    PAV_VCF_PRESENT=$(gsutil -q stat ${BUCKET_DIR}/vcfs/${PREFIX_PAV}.vcf && echo 0 || echo 1)
    if [ ${USE_HIFIASM} -eq 1 -o ${USE_PAV} -eq 1 ]; then 
        TEST1=$(gsutil -q stat ${BUCKET_DIR}/assemblies/${PREFIX_ASSEMBLY}_h1.fa && echo 0 || echo 1)
        TEST2=$(gsutil -q stat ${BUCKET_DIR}/assemblies/${PREFIX_ASSEMBLY}_h2.fa && echo 0 || echo 1)
        if [ ${TEST1} -eq 0 -a ${TEST2} -eq 0 ]; then
            if [ ${USE_PAV} -eq 1 -a ${PAV_VCF_PRESENT} -eq 1 ]; then
                while : ; do
                    TEST=$(gsutil cp ${BUCKET_DIR}/assemblies/${PREFIX_ASSEMBLY}_h1.fa ${PREFIX_ASSEMBLY}_h1.fa && echo 0 || echo 1)
                    if [ ${TEST} -eq 1 ]; then
                        echo "Error downloading file <${BUCKET_DIR}/assemblies/${PREFIX_ASSEMBLY}_h1.fa>. Trying again..."
                        sleep ${GSUTIL_DELAY_S}
                    else
                        break
                    fi
                done
                while : ; do
                    TEST=$(gsutil cp ${BUCKET_DIR}/assemblies/${PREFIX_ASSEMBLY}_h2.fa ${PREFIX_ASSEMBLY}_h2.fa && echo 0 || echo 1)
                    if [ ${TEST} -eq 1 ]; then
                        echo "Error downloading file <${BUCKET_DIR}/assemblies/${PREFIX_ASSEMBLY}_h2.fa>. Trying again..."
                        sleep ${GSUTIL_DELAY_S}
                    else
                        break
                    fi
                done
            fi
        else
            ${TIME_COMMAND} hifiasm -t ${N_THREADS} --hom-cov $(( ${COVERAGE}*2 )) -o tmpasm coverage_${COVERAGE}.fa
            awk '/^S/{print ">"$2; print $3}' tmpasm.bp.hap1.p_ctg.gfa > ${PREFIX_ASSEMBLY}_h1.fa
            awk '/^S/{print ">"$2; print $3}' tmpasm.bp.hap2.p_ctg.gfa > ${PREFIX_ASSEMBLY}_h2.fa
            while : ; do
                TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp ${PREFIX_ASSEMBLY}_h1.fa ${BUCKET_DIR}/assemblies/ && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error uploading file <${PREFIX_ASSEMBLY}_h1.fa>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
            while : ; do
                TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp ${PREFIX_ASSEMBLY}_h2.fa ${BUCKET_DIR}/assemblies/ && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error uploading file <${PREFIX_ASSEMBLY}_h2.fa>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else 
                    break
                fi
            done
            awk '/^S/{print ">"$2; print $3}' tmpasm.bp.p_ctg.gfa > ${PREFIX_ASSEMBLY}_primaryContigs.fa
            quast --threads ${N_THREADS} --eukaryote --large --est-ref-size ${REFERENCE_LENGTH} --no-gc ${PREFIX_ASSEMBLY}_primaryContigs.fa
            rm -f ${PREFIX_ASSEMBLY}_primaryContigs.fa
            tar -czf ${PREFIX_ASSEMBLY}_quast.tar.gz quast_results/ 
            rm -rf quast_results/
            while : ; do
                TEST=$(gsutil cp ${PREFIX_ASSEMBLY}_quast.tar.gz ${BUCKET_DIR}/assemblies/ && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error uploading file <${PREFIX_ASSEMBLY}_quast.tar.gz>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else 
                    break
                fi
            done
            rm -rf tmpasm*
        fi
    fi
    
    # PAV
    if [ ${USE_PAV} -eq 1 -a ${PAV_VCF_PRESENT} -eq 1 ]; then
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
        bgzip -@ ${N_THREADS} --compress-level 0 -c ${PREFIX_ASSEMBLY}_h1.fa > asm/sample/h1.fa.gz
        bgzip -@ ${N_THREADS} --compress-level 0 -c ${PREFIX_ASSEMBLY}_h2.fa > asm/sample/h2.fa.gz
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
        bcftools filter --exclude 'SVTYPE="SNV" || (SVLEN>-40 && SVLEN<40)' pav_sample.vcf.gz > ${PREFIX_PAV}.vcf
        rm -f pav_sample.vcf.gz*
        while : ; do
            TEST=$(gsutil cp ${PREFIX_PAV}.vcf ${BUCKET_DIR}/vcfs/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading file <${PREFIX_PAV}.vcf>. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        rm -rf asm/ data/ temp/ results/ log/
        gsutil -m rm -rf "${BUCKET_DIR}/pav/${INFIX}" || echo "Error removing directory <${BUCKET_DIR}/pav/${INFIX}>"
    fi
    rm -f ${PREFIX_ASSEMBLY}*
    rm -f ${PREFIX_PAV}*
    
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
