#!/bin/bash
#
BUCKET_DIR=$1
CHILD_ID=$2
TRUTH_VCF_PREFIX=$3  # String before <_caller>
PRECURSOR_VCF_PREFIX=$4  # Of the individuals used to build the ground truth
MEASURED_CHARACTER_CODE=$5  # Single character after <reads_>
CALLERS=$6
VALUES=$7  # Separated by -
SV_LENGTHS=$8  # Separated by -. Increasing. Must not include zero.
ONLY_PASS=$9
WORK_DIR=${10}

set -euxo pipefail
cd ${WORK_DIR}

GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
GSUTIL_DELAY_S="600"
CALLERS=$(echo ${CALLERS} | tr '-' ' ')
VALUES=$(echo ${VALUES} | tr '-' ' ')
SV_LENGTHS=$(echo ${SV_LENGTHS} | tr '-' ' ')


# Builds the histogram of one VCF
#
function histogramThread() {
    local REMOTE_VCF_ADDRESS=$1
    local LOCAL_ID=$2  # Unique among all parallel jobs
    local IS_GZ=$3
    
    # Downloading the VCF
    TEST=$(gsutil -q stat "${REMOTE_VCF_ADDRESS}" && echo 0 || echo 1)
    if [ ${TEST} -eq 1 ]; then
        return
    fi
    if [ ${IS_GZ} -eq 1 ]; then
        DOWNLOADED_FILE=${LOCAL_ID}.vcf.gz
    else
        DOWNLOADED_FILE=${LOCAL_ID}.vcf
    fi
    while : ; do
        TEST=$(gsutil cp "${REMOTE_VCF_ADDRESS}" ./${DOWNLOADED_FILE} && echo 0 || echo 1)
        if [ ${TEST} -eq 1 ]; then
            echo "Error downloading file <${REMOTE_VCF_ADDRESS}>. Trying again..."
            sleep ${GSUTIL_DELAY_S}
        else
            break
        fi
    done
    if [ ${IS_GZ} -eq 1 ]; then
        gunzip ${DOWNLOADED_FILE}
    fi
    
    # Fixing the VCF format issues of some callers
    bcftools view -h ${LOCAL_ID}.vcf > ${LOCAL_ID}.header.txt
    tail -n 1 ${LOCAL_ID}.header.txt > ${LOCAL_ID}.columns.txt
    N_LINES=$(wc -l < ${LOCAL_ID}.header.txt)
    head -n $(( ${N_LINES} - 1 )) ${LOCAL_ID}.header.txt > ${LOCAL_ID}.newFile.vcf
    echo "##FILTER=<ID=STRANDBIAS,Description=\"STRANDBIAS\">" >> ${LOCAL_ID}.newFile.vcf
    cat ${LOCAL_ID}.columns.txt >> ${LOCAL_ID}.newFile.vcf
    bcftools view -H ${LOCAL_ID}.vcf >> ${LOCAL_ID}.newFile.vcf
    rm -f ${LOCAL_ID}.vcf
    mv ${LOCAL_ID}.newFile.vcf ${LOCAL_ID}.vcf
    rm -f ${LOCAL_ID}.header.txt ${LOCAL_ID}.columns.txt
    bcftools sort --output-type z --output ${LOCAL_ID}.vcf.gz ${LOCAL_ID}.vcf
    rm -f ${LOCAL_ID}.vcf
    tabix ${LOCAL_ID}.vcf.gz
    
    # Building the histogram
    OUTPUT_FILE=${LOCAL_ID}.histogram
    rm -f ${OUTPUT_FILE}; touch ${OUTPUT_FILE}
    PREVIOUS_SV_LENGTH="0"
    for sv_length in ${SV_LENGTHS}; do
        FILTER_STRING="(SVLEN>0 && SVLEN>${PREVIOUS_SV_LENGTH} && SVLEN<=${sv_length}) || (SVLEN<0 && SVLEN<-${PREVIOUS_SV_LENGTH} && SVLEN>=-${sv_length})"
        if [ ${ONLY_PASS} -eq 1 ]; then
            FILTER_STRING="(${FILTER_STRING}) && FILTER=\"PASS\""
        fi
        N_VARIANTS=$(bcftools filter --threads 0 --include "${FILTER_STRING}" --output-type v ${LOCAL_ID}.vcf.gz | bcftools view -H | wc -l)
        echo -n "${N_VARIANTS}," >> ${OUTPUT_FILE}
        PREVIOUS_SV_LENGTH=${sv_length}
    done
    echo "" >> ${OUTPUT_FILE}
}


# Main program
#
for caller in ${CALLERS}; do
    while : ; do
        TEST=$(gsutil cp ${BUCKET_DIR}/trios_info/${CHILD_ID}.parents . && echo 0 || echo 1)
        if [ ${TEST} -eq 1 ]; then
            echo "Error downloading file <${BUCKET_DIR}/trios_info/${CHILD_ID}.parents>. Trying again..."
            sleep ${GSUTIL_DELAY_S}
        else
            break
        fi
    done
    i=0
    while read PARENT_ID; do
        histogramThread "${BUCKET_DIR}/${CHILD_ID}/${PRECURSOR_VCF_PREFIX}_${PARENT_ID}/vcfs/${caller}_${PARENT_ID}.vcf" precursor_${i} 0 &
        i=$(( $i + 1 ))
    done < ${CHILD_ID}.parents
    histogramThread "${BUCKET_DIR}/${CHILD_ID}/${TRUTH_VCF_PREFIX}_${caller}_truth.vcf.gz" truth 1 &
    for value in ${VALUES}; do
        vcfThread "${BUCKET_DIR}/${CHILD_ID}/reads_${MEASURED_CHARACTER_CODE}${value}/vcfs/${CALLER}_${CHILD_ID}.vcf" measured_${value} 0 &
    done
    wait
    rm -f ${caller}_svLengths.histogram
    cat truth.histogram precursor_*.histogram > ${caller}_svLengths.histogram
    for value in ${VALUES}; do
        cat measured_${value}.histogram > ${caller}_svLengths.histogram
    done
    while : ; do
        TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp ${caller}_svLengths.histogram ${BUCKET_DIR}/${CHILD_ID}/ && echo 0 || echo 1)
        if [ ${TEST} -eq 1 ]; then
            echo "Error uploading SV length histograms (${caller}). Trying again..."
            sleep ${GSUTIL_DELAY_S}
        else
            break
        fi
    done
done