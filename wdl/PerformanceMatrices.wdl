version 1.0


# 
workflow PerformanceMatrices {
    input {
        String bucket_dir
        Array[String] callers
        Array[Int] coverages
        Array[Int] lengths

        File reference_fa

    }
    parameter_meta {
        bucket_dir: "a"
    }
    Array[String] svTypes = ["ANY", "DEL", "DUP", "INS", "INV"]
    Array[Int] contextTypes = [0, 1, 2, 3, 4]
    # 0: non-satellite repeat;
    # 1: satellite repeat;
    # 2: both of the above;
    # 3: none of the above, but segmental duplication;
    # 4: none of the above.
    Array[Float] repeatFractions = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    call GetChunks {
        input:
            callers = callers,
            svTypes = svTypes,
            contextTypes = contextTypes,
            repeatFractions = repeatFractions
            n_chunks = n_chunks
    }
    scatter(chunk in GetChunks.chunks) {
        ----------->
        call JointCallingImpl { 
            input:
                task_description = description,
                bucket_dir = bucket_dir,
                reference_fa = reference_fa,
                n_individuals = n_individuals,
                max_signature_file_size = max_signature_file_size,
                n_cpus = n_cpus
        }
    }
    output {
    }
}


# Collects all workpackage IDs into the rows of a list with format:
#
# <caller svType repeatFraction previousRepeatFraction contextTypeStart
# contextTypeEnd>
#
# and splits such a list into a given number of chunks.
task GetChunks {
    input {
        Array[String] callers
        Array[String] svTypes
        Array[Int] contextTypes
        Array[Float] repeatFractions
        Int n_chunks
    }
    command <<<
        set -euxo pipefail
        
        WORKPACKAGES_FILE="workpackages.txt"
        CALLERS=~{sep='-' callers}
        SV_TYPES=~{sep='-' svTypes}
        CONTEXT_TYPES=~{sep='-' contextTypes}
        REPEAT_FRACTIONS=~{sep='-' repeatFractions}
        for caller in ${CALLERS}; do
            # Specific SV type
            for svType in ${SV_TYPES}; do
                # Specific context
                for contextTypeStart in ${CONTEXT_TYPES}; do
                    for contextTypeEnd in ${CONTEXT_TYPES}; do
                        echo "${caller} ${svType} -1 -1 ${contextTypeStart} ${contextTypeEnd}" >> ${WORKPACKAGES_FILE}
                    done
                done
                previousRepeatFraction="0"
                for repeatFraction in ${REPEAT_FRACTIONS}; do
                    echo "${caller} ${svType} ${repeatFraction} ${previousRepeatFraction} -1 -1" >> ${WORKPACKAGES_FILE}
                    previousRepeatFraction=${repeatFraction}
                done
                # Any context
                echo "${caller} ${svType} -1 -1 -1 -1" >> ${WORKPACKAGES_FILE}
            done
            # Any SV type
            # Specific context
            for contextTypeStart in ${CONTEXT_TYPES}; do
                for contextTypeEnd in ${CONTEXT_TYPES}; do
                    echo "${caller} -1 -1 -1 ${contextTypeStart} ${contextTypeEnd}" >> ${WORKPACKAGES_FILE}
                done
            done
            previousRepeatFraction="0"
            for repeatFraction in ${REPEAT_FRACTIONS}; do
                echo "${caller} -1 ${repeatFraction} ${previousRepeatFraction} -1 -1" >> ${WORKPACKAGES_FILE}
                previousRepeatFraction=${repeatFraction}
            done
            # Any context
            echo "${caller} -1 -1 -1 -1 -1" >> ${WORKPACKAGES_FILE}
        done
        split -d -n ~{n_chunks} ${WORKPACKAGES_FILE} chunk-
    >>>
    output {
        Array[File] chunks = glob("chunk-*")
    }
    runtime {
        docker: "ubuntu:latest"
    }
}


# Processes every workpackage in a given chunk sequentially.
task ProcessChunk {
    input {
        File chunk_file
        Array[Int] read_lengths
        Array[Int] coverages
        File reference_fa
        Float max_vcf_size
        Int n_individuals
        String bucket_dir_measured_vcfs
        String bucket_dir_ground_truth_vcfs
        String bucket_dir_matrices
        Int n_cpus
    }
    parameter_meta {
        chunk_file: "Every row is a workpackage ID with format <caller svType repeatFraction previousRepeatFraction contextTypeStart contextTypeEnd>."
        n_individuals: "Total number of diploid individuals in the population."
        max_vcf_size: "Max size of a single VCF file, in GB."
        n_cpus: "The program proceeds through sequential iterations, but in each iteration it processes all VCF files using all the cores of the machine."
    }
    Int ram_size_gb_pre = max_vcf_size*4
    Int ram_size_gb = if 16 > ram_size_gb_pre then 16 else ram_size_gb_pre
    Int disk_size_gb = ceil(max_vcf_size*n_individuals*2) + ceil((size(reference_fa, "GB")))
    String docker_dir = "/simulation"
    String work_dir = "/cromwell_root/simulation"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        echo "Running <ProcessChunk> on ${N_THREADS} cores on the following node:"
        lscpu
        cat /proc/meminfo
        
        mkdir -p ground_truth_vcfs/
        ${TIME_COMMAND} gsutil -m cp "~{bucket_dir_ground_truth_vcfs}/groundTruth_individual_*.vcf" ground_truth_vcfs/
        while read LINE; do
            caller=$(echo ${LINE} | awk '{print $1}')
            svType=$(echo ${LINE} | awk '{print $2}')
            repeatFraction=$(echo ${LINE} | awk '{print $3}')
            previousRepeatFraction=$(echo ${LINE} | awk '{print $4}')
            contextTypeStart=$(echo ${LINE} | awk '{print $5}')
            contextTypeEnd=$(echo ${LINE} | awk '{print $6}')
            FILTER_STRING=""; PREFIX=""
            if [ ${svType} -ne -1 ]; then
                FILTER_STRING="SVTYPE=\"${svType}\""
                PREFIX="${caller}_sv${svType}"
            else
                FILTER_STRING=""
                PREFIX="${caller}"
            fi
            if [ ${repeatFraction} -eq -1 ]; then
                if [ ${contextTypeStart} -ne -1 ]; then
                    if [ ${#FILTER_STRING} -eq 0 ]; then
                        FILTER_STRING="REPEATS_START=${contextTypeStart} && REPEATS_END=${contextTypeEnd}"
                    else
                        FILTER_STRING="${FILTER_STRING} && REPEATS_START=${contextTypeStart} && REPEATS_END=${contextTypeEnd}"
                    fi
                    PREFIX="${PREFIX}_rs${contextTypeStart}_re${contextTypeEnd}"
                else
                    :  # NOP
                fi
            else
                if [ ${#FILTER_STRING} -eq 0 ]; then
                    FILTER_STRING="REPEATS_FRACTION>${previousRepeatFraction} && REPEATS_FRACTION<=${repeatFraction}"
                else
                    FILTER_STRING="${FILTER_STRING} && REPEATS_FRACTION>${previousRepeatFraction} && REPEATS_FRACTION<=${repeatFraction}"
                fi
                PREFIX="${PREFIX}_rf${repeatFraction}"
            fi
            if [ ${#FILTER_STRING} -eq 0 ]; then
                FILTER_STRING="+"
            else
                FILTER_STRING=$(echo ${FILTER_STRING} | tr ' ' '+')
            fi
            TP_MATRIX="${PREFIX}_matrix_tp.txt"
            FP_MATRIX="${PREFIX}_matrix_fp.txt"
            FN_MATRIX="${PREFIX}_matrix_fn.txt"
            PRECISION_MATRIX="${PREFIX}_matrix_precision.txt"
            RECALL_MATRIX="${PREFIX}_matrix_recall.txt"
            F1_MATRIX="${PREFIX}_matrix_f1.txt"
            TEST1=$(gsutil -q stat ${TP_MATRIX} && echo 0 || echo 1)
            TEST2=$(gsutil -q stat ${FP_MATRIX} && echo 0 || echo 1)
            TEST3=$(gsutil -q stat ${FN_MATRIX} && echo 0 || echo 1)
            TEST4=$(gsutil -q stat ${PRECISION_MATRIX} && echo 0 || echo 1)
            TEST5=$(gsutil -q stat ${RECALL_MATRIX} && echo 0 || echo 1)
            TEST6=$(gsutil -q stat ${F1_MATRIX} && echo 0 || echo 1)
            if [ ${TEST1} -eq 0 -a ${TEST2} -eq 0 -a ${TEST3} -eq 0 -a ${TEST4} -eq 0 -a ${TEST5} -eq 0 -a ${TEST6} -eq 0 ]; then
                continue
            fi
            for readLength in ~{read_lengths}; do
                for coverage in ~{coverages}; do
                    echo -n "${READ_LENGTH},${COVERAGE}," >> ${TP_MATRIX}
                    echo -n "${READ_LENGTH},${COVERAGE}," >> ${FP_MATRIX}
                    echo -n "${READ_LENGTH},${COVERAGE}," >> ${FN_MATRIX}
                    echo -n "${READ_LENGTH},${COVERAGE}," >> ${PRECISION_MATRIX}
                    echo -n "${READ_LENGTH},${COVERAGE}," >> ${RECALL_MATRIX}
                    echo -n "${READ_LENGTH},${COVERAGE}," >> ${F1_MATRIX}
                    rm -rf measured_vcfs/; mkdir -p measured_vcfs/
                    ${TIME_COMMAND} gsutil -m cp "${BUCKET_DIR_MEASURED_VCFS}/${CALLER}_*_l${READ_LENGTH}_c${COVERAGE}.vcf" measured_vcfs/
                    bash ~{docker_dir}/performance_matrices_impl.sh ${FILTER_STRING} ~{reference_fa} ${N_THREADS} ~{work_dir} ${TP_MATRIX} ${FP_MATRIX} ${FN_MATRIX} ${PRECISION_MATRIX} ${RECALL_MATRIX} ${F1_MATRIX}
                    echo "" >> ${TP_MATRIX}
                    echo "" >> ${FP_MATRIX}
                    echo "" >> ${FN_MATRIX}
                    echo "" >> ${PRECISION_MATRIX}
                    echo "" >> ${RECALL_MATRIX}
                    echo "" >> ${F1_MATRIX}
                done
            done
            ${TIME_COMMAND} gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp ${TP_MATRIX} ${FP_MATRIX} ${FN_MATRIX} ${PRECISION_MATRIX} ${RECALL_MATRIX} ${F1_MATRIX} ~{bucket_dir_matrices}
        done < ~{chunk_file}
    >>>
    
    output {
    }
    runtime {
        docker: "fcunial/simulation"
        cpu: n_cpus
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 1
    }
}
