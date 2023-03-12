#!/bin/bash
#
# Calls several SV discovery tools on the same set of reads.
#
DATASET_ID=$1  # A unique identifier of the dataset
READS_BAM=$2
READS_FA=$3  # FASTA or FASTQ
COVERAGE=$4  # Of each haplotype
SAMPLE_ID=$5  # SM field in the .sam file (needed later for joint calling)
N_THREADS=$6
REFERENCE_FA=$7
REFERENCE_FAI=$8
REFERENCE_TANDEM_REPEATS=$9
BUCKET_DIR=${10}  # Remote bucket dir where all output files will be stored
USE_PBSV=${11}
USE_SNIFFLES1=${12}
USE_SNIFFLES2=${13}
USE_HIFIASM=${14}
USE_PAV=${15}
USE_PAFTOOLS=${16}
KEEP_ASSEMBLIES=${17}
WORK_DIR=${18}
DOCKER_DIR=${19}
PBSV_PARALLELIZE_BY_CHROMOSOME=${20}

GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
GSUTIL_DELAY_S="600"
TIME_COMMAND="/usr/bin/time --verbose"
REFERENCE_LENGTH=$(cut -f 2 ${REFERENCE_FAI} | awk '{s+=$1} END {print s}')
REFERENCE_LENGTH=$(printf '%d' ${REFERENCE_LENGTH})
MINIMAP2_ASM_FLAG="asm5"  # asm5/asm10/asm20 for ~0.1/1/5% sequence divergence
PAFTOOLS_MIN_SV_LENGTH="40"

set -euxo pipefail
cd ${WORK_DIR}

# PBSV
PREFIX_PBSV="pbsv_${DATASET_ID}"
if [ ${USE_PBSV} -eq 1 ]; then
    samtools view -H ${READS_BAM} | grep '^@SQ' | cut -f2 | cut -d':' -f2 > contigs.txt
    if [ ${PBSV_PARALLELIZE_BY_CHROMOSOME} -eq 1 ]; then
        contig=$(head -n 1 contigs.txt)
        TEST=$(gsutil -q stat ${BUCKET_DIR}/signatures/${PREFIX_PBSV}.${contig}.svsig.gz && echo 0 || echo 1)
    else
        TEST=$(gsutil -q stat ${BUCKET_DIR}/signatures/${PREFIX_PBSV}.svsig.gz && echo 0 || echo 1)
    fi
    if [ ${TEST} -eq 0 ]; then
        while : ; do
            TEST=$(gsutil -m cp "${BUCKET_DIR}/signatures/${PREFIX_PBSV}*.svsig.gz" . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading files <${BUCKET_DIR}/signatures/${PREFIX_PBSV}*.svsig.gz>. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
    else
        # <discover> is sequential
        if [ ${PBSV_PARALLELIZE_BY_CHROMOSOME} -eq 1 ]; then
            while read CONTIG; do
                ${TIME_COMMAND} pbsv discover --region ${CONTIG} --tandem-repeats ${REFERENCE_TANDEM_REPEATS} ${READS_BAM} ${PREFIX_PBSV}.${CONTIG}.svsig.gz
            done < contigs.txt
        else
        	${TIME_COMMAND} pbsv discover --tandem-repeats ${REFERENCE_TANDEM_REPEATS} ${READS_BAM} ${PREFIX_PBSV}.svsig.gz
        fi
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp "${PREFIX_PBSV}*.svsig.gz" ${BUCKET_DIR}/signatures/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading file <${PREFIX_PBSV}.svsig.gz>. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
    fi
    TEST=$(gsutil -q stat ${BUCKET_DIR}/vcfs/${PREFIX_PBSV}.vcf && echo 0 || echo 1)
    if [ ${TEST} -eq 1 ]; then
    	${TIME_COMMAND} pbsv call -j ${N_THREADS} --ccs ${REFERENCE_FA} ${PREFIX_PBSV}*.svsig.gz ${PREFIX_PBSV}.vcf
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp ${PREFIX_PBSV}.vcf ${BUCKET_DIR}/vcfs/ && echo 0 || echo 1)
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
PREFIX_SNIFFLES1="sniffles1_${DATASET_ID}"
if [ ${USE_SNIFFLES1} -eq 1 ]; then 
    TEST=$(gsutil -q stat ${BUCKET_DIR}/vcfs/${PREFIX_SNIFFLES1}.vcf && echo 0 || echo 1)
    if [ ${TEST} -eq 1 ]; then
        ${TIME_COMMAND} sniffles1 -t ${N_THREADS} -m ${READS_BAM} -v ${PREFIX_SNIFFLES1}.vcf
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp ${PREFIX_SNIFFLES1}.vcf ${BUCKET_DIR}/vcfs/ && echo 0 || echo 1)
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
PREFIX_SNIFFLES2="sniffles2_${DATASET_ID}"
if [ ${USE_SNIFFLES2} -eq 1 ]; then 
    TEST=$(gsutil -q stat ${BUCKET_DIR}/signatures/${PREFIX_SNIFFLES2}.snf && echo 0 || echo 1)
    if [ ${TEST} -eq 1 ]; then
	    ${TIME_COMMAND} sniffles --threads ${N_THREADS} --tandem-repeats ${REFERENCE_TANDEM_REPEATS} --reference ${REFERENCE_FA} --sample-id ${SAMPLE_ID} --input ${READS_BAM} --vcf ${PREFIX_SNIFFLES2}.vcf --snf ${PREFIX_SNIFFLES2}.snf
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp ${PREFIX_SNIFFLES2}.snf ${BUCKET_DIR}/signatures/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading file <${PREFIX_SNIFFLES2}.snf>. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else 
                break
            fi
        done
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp ${PREFIX_SNIFFLES2}.vcf ${BUCKET_DIR}/vcfs/ && echo 0 || echo 1)
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
PREFIX_ASSEMBLY="assembly_${DATASET_ID}"
PREFIX_PAV="pav_${DATASET_ID}"
PREFIX_PAFTOOLS="paftools_${DATASET_ID}"
PAV_VCF_PRESENT=$(gsutil -q stat ${BUCKET_DIR}/vcfs/${PREFIX_PAV}.vcf && echo 0 || echo 1)
PAFTOOLS_VCF_PRESENT=$(gsutil -q stat ${BUCKET_DIR}/vcfs/${PREFIX_PAFTOOLS}.vcf && echo 0 || echo 1)
if [ ${USE_HIFIASM} -eq 1 -o ${USE_PAV} -eq 1 -o ${USE_PAFTOOLS} -eq 1 ]; then 
    TEST1=$(gsutil -q stat ${BUCKET_DIR}/assemblies/${PREFIX_ASSEMBLY}_h1.fa && echo 0 || echo 1)
    TEST2=$(gsutil -q stat ${BUCKET_DIR}/assemblies/${PREFIX_ASSEMBLY}_h2.fa && echo 0 || echo 1)
    if [ ${TEST1} -eq 0 -a ${TEST2} -eq 0 ]; then
        if [ ${USE_PAV} -eq 1 -a ${PAV_VCF_PRESENT} -eq 1 ] || [ ${USE_PAFTOOLS} -eq 1 -a ${PAFTOOLS_VCF_PRESENT} -eq 1 ]; then
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
        ${TIME_COMMAND} hifiasm -t ${N_THREADS} --hom-cov $(echo "scale=8; ${COVERAGE}*2.0" | bc) -o tmpasm ${READS_FA}
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

# PAFTOOLS
if [ ${USE_PAFTOOLS} -eq 1 -a ${PAFTOOLS_VCF_PRESENT} -eq 1 ]; then
    TEST=$(gsutil -q stat ${BUCKET_DIR}/vcfs/${PREFIX_PAFTOOLS}.vcf && echo 0 || echo 1)
    if [ ${TEST} -eq 1 ]; then
        cat ${PREFIX_ASSEMBLY}_h1.fa ${PREFIX_ASSEMBLY}_h2.fa > ${PREFIX_ASSEMBLY}_h1_h2.fa
        ${TIME_COMMAND} minimap2 -t ${N_THREADS} -x ${MINIMAP2_ASM_FLAG} -c --cs ${REFERENCE_FA} ${PREFIX_ASSEMBLY}_h1_h2.fa -o ${PREFIX_ASSEMBLY}_h1_h2.paf
        rm -f ${PREFIX_ASSEMBLY}_h1_h2.fa
        ${TIME_COMMAND} sort --parallel=${N_THREADS} -k6,6 -k8,8n ${PREFIX_ASSEMBLY}_h1_h2.paf > ${PREFIX_ASSEMBLY}_h1_h2.sorted.paf
        rm -f ${PREFIX_ASSEMBLY}_h1_h2.paf
        ${TIME_COMMAND} paftools.js call -f ${REFERENCE_FA} ${PREFIX_ASSEMBLY}_h1_h2.sorted.paf > ${PREFIX_PAFTOOLS}_raw.vcf
        rm -f ${PREFIX_ASSEMBLY}_h1_h2.sorted.paf
        java -cp ${DOCKER_DIR} FormatPaftoolsVCF ${PREFIX_PAFTOOLS}_raw.vcf ${PREFIX_PAFTOOLS}.vcf ${PAFTOOLS_MIN_SV_LENGTH}
        rm -f ${PREFIX_PAFTOOLS}_raw.vcf
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp ${PREFIX_PAFTOOLS}.vcf ${BUCKET_DIR}/vcfs/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading file <${PREFIX_PAFTOOLS}.vcf>. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else 
                break
            fi
        done
    fi
fi
rm -f ${PREFIX_PAFTOOLS}*

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
    bash ${DOCKER_DIR}/pav_restoreCheckpoint.sh ${BUCKET_DIR} ${DATASET_ID} ${WORK_DIR} ${DOCKER_DIR}
    bash ${DOCKER_DIR}/pav_checkpointDaemon.sh ${BUCKET_DIR} ${DATASET_ID} ${WORK_DIR} ${DOCKER_DIR} &
    DAEMON_ID=$!
    source activate lr-pav
    ${TIME_COMMAND} snakemake -s ${DOCKER_DIR}/pav/Snakefile --cores ${N_THREADS} pav_sample.vcf.gz
    conda deactivate
    kill -KILL ${DAEMON_ID} || echo "The process ID <${DAEMON_ID}> of <pav_checkpointDaemon.sh> cannot be found."
    bcftools filter --exclude 'SVTYPE="SNV" || (SVLEN>-40 && SVLEN<40)' pav_sample.vcf.gz > ${PREFIX_PAV}.vcf
    rm -f pav_sample.vcf.gz*
    while : ; do
        TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp ${PREFIX_PAV}.vcf ${BUCKET_DIR}/vcfs/ && echo 0 || echo 1)
        if [ ${TEST} -eq 1 ]; then
            echo "Error uploading file <${PREFIX_PAV}.vcf>. Trying again..."
            sleep ${GSUTIL_DELAY_S}
        else
            break
        fi
    done
    rm -rf asm/ data/ temp/ results/ log/
    gsutil -m rm -rf "${BUCKET_DIR}/pav/${DATASET_ID}" || echo "Error removing directory <${BUCKET_DIR}/pav/${DATASET_ID}>"
fi
rm -f ${PREFIX_PAV}*

rm -f ${PREFIX_ASSEMBLY}*
