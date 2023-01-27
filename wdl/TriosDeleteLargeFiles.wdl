version 1.0

import "TriosTruthVCFUtils.wdl"


# Removes just the largest files that were cached by the trio analysis pipeline
# to handle preemption or failure.
#
workflow TriosDeleteLargeFiles {
    input {
        File children_ids
        String bucket_dir
        Int answer_question1
        Array[Float] question1_left_weights
        Int answer_question2
        Array[Float] question2_left_coverages
    }
    parameter_meta {
        children_ids: "The trio children to be processed"
    }
    
    call TriosDeleteLargeFilesImpl { 
        input:
            children_ids = children_ids,
            bucket_dir = bucket_dir,
            answer_question1 = answer_question1,
            question1_left_weights = question1_left_weights,
            answer_question2 = answer_question2,
            question2_left_coverages = question2_left_coverages
    }
    output {
    }
}


task TriosDeleteLargeFilesImpl {
    input {
        File children_ids
        String bucket_dir
        Int answer_question1
        Array[Float] question1_left_weights
        Int answer_question2
        Array[Float] question2_left_coverages
    }
    
    command <<<
        set -ux
        
        GSUTIL_DELAY_S="600"
        
        while read CHILD_ID; do
            while : ; do
                TEST=$(gsutil cp ~{bucket_dir}/trios_info/${CHILD_ID}.parents . && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading file <~{bucket_dir}/trios_info/${CHILD_ID}.parents>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done  
            if [ ~{answer_question1} -eq 1 ]; then
                :  # To be implemented...
            fi
            if [ ~{answer_question2} -eq 1 ]; then
                VALUES=~{sep='-' question2_left_coverages}
                VALUES=$(echo ${VALUES} | tr '-' ' ')
                gsutil -m rm "~{bucket_dir}/${CHILD_ID}/bins/bin_*.bin"
                while read PARENT; do
                    DIRECTORY="~{bucket_dir}/${CHILD_ID}/long_coverage_${PARENT}"
                    gsutil -m rm ${DIRECTORY}/reads.bam ${DIRECTORY}/reads.bam.bai ${DIRECTORY}/reads.fastq
                done < ${CHILD_ID}.parents
                DIRECTORY="~{bucket_dir}/${CHILD_ID}/long_coverage_${CHILD_ID}"
                gsutil -m rm ${DIRECTORY}/reads.bam ${DIRECTORY}/reads.bam.bai ${DIRECTORY}/reads.fastq
                for VALUE in ${VALUES}; do
                    DIRECTORY="~{bucket_dir}/${CHILD_ID}/reads_c${VALUE}"
                    TEST=$(gsutil -q stat ${DIRECTORY}/coverage_${VALUE}.bam && echo 0 || echo 1)
                    if [ ${TEST} -eq 0 ]; then
                        gsutil -m rm ${DIRECTORY}/coverage_${VALUE}.bam ${DIRECTORY}/coverage_${VALUE}.bam.bai ${DIRECTORY}/coverage_${VALUE}.fastq
                    fi
                done
                DIRECTORY="~{bucket_dir}/${CHILD_ID}"
                gsutil -m rm ${DIRECTORY}/reads_maxCoverage_left.fastq ${DIRECTORY}/reads_maxCoverage_right.fastq 
                gsutil -m rm -r ${DIRECTORY}/reads_question2_chunks/
            fi
        done < ~{children_ids}
    >>>
    
    output {
    }
    runtime {
        docker: "fcunial/simulation"
        cpu: 1
        preemptible: 0
    }
}
