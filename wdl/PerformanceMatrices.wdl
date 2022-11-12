version 1.0

import "AnnotateVCFs.wdl"


# The program: 
# 1. Compares each experimental VCF to its corresponding ground truth,
#    collecting per-individual performance measures. 
# 2. Merges all experimental VCFs and compares the merge to the set of distinct
#    ground-truth variants in the population. 
# 3. Compares the joint-calling files (produced by some callers) to the set of
#    distinct ground-truth variants in the population.
# To enable granular data analysis downstream, this is performed separately for
# every tuple (caller, X, repeat context), where X is an SV type, an SV length,
# or an SV frequency.
#
# Remark: this workflow can be run while the VCF files are being created by the
# simulation workflow.
#
workflow PerformanceMatrices {
    input {        
        Array[String] callers
        Array[String] svTypes
        Array[Int] svLengths
        Array[Float] svFrequencies
        Array[Int] contextTypes
        Array[Float] repeatFractions
        Array[Int] read_lengths
        Array[Int] coverages
        Int only_pass
        Int n_nodes
        Int n_cpus_per_node
        String reference_fa
        String reference_fai
        String bucket_dir_experimental_vcfs
        String bucket_dir_ground_truth_vcfs
        String bucket_dir_allIndividuals_vcfs
        String bucket_dir_matrices
        Float max_vcf_size
        Int n_individuals
        Int annotate_vcfs
        Int annotate_vcfs_min_distance
        Int annotate_vcfs_max_distance
        File annotate_vcfs_repeat_intervals
        File annotate_vcfs_segdups
    }
    parameter_meta {
        contextTypes: "0=non-satellite repeat; 1=satellite repeat; 2=both of the above; 3=none of the above, but segmental duplication; 4=none of the above."
        repeatFractions: "Fractions of one."
        only_pass: "(0/1) Use only variants with FILTER=PASS."
        n_nodes: "Number of machines over which to distribute the computation."
        reference_fa: "Address in a bucket"
        reference_fai: "Address in a bucket"
        bucket_dir_experimental_vcfs: "Input directory"
        bucket_dir_ground_truth_vcfs: "Input directory"
        bucket_dir_allIndividuals_vcfs: "Output directory: the merged VCFs over all individuals will be stored here."
        bucket_dir_matrices: "Output directory: the matrices with performance measures will be stored here."
        max_vcf_size: "Max size of a single VCF file, in GB."
        n_individuals: "Total number of diploid individuals in the population. Used just to estimate space."
    }
    if (annotate_vcfs == 1) {
        call AnnotateVCFs.AnnotateVCFs {
            input:
                bucket_dir = bucket_dir_experimental_vcfs,
                max_vcf_size = max_vcf_size,
                min_distance = annotate_vcfs_min_distance,
                max_distance = annotate_vcfs_max_distance,
                reference_fa = reference_fa,
                repeat_intervals = annotate_vcfs_repeat_intervals,
                segdups = annotate_vcfs_segdups,
                n_nodes = n_nodes
        }
    }
    call GetChunks {
        input:
            callers = callers,
            svTypes = svTypes,
            svLengths = svLengths,
            svFrequencies = svFrequencies,
            contextTypes = contextTypes,
            repeatFractions = repeatFractions,
            n_chunks = n_nodes,
            force_sequentiality = if annotate_vcfs==1 then AnnotateVCFs.force_sequentiality else [1]
    }
    scatter(chunk_file in GetChunks.chunks) {
        call ProcessChunk { 
            input:
                chunk_file = chunk_file,
                read_lengths = read_lengths,
                coverages = coverages,
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                max_vcf_size = max_vcf_size,
                n_individuals = n_individuals,
                only_pass = only_pass,
                bucket_dir_experimental_vcfs = if annotate_vcfs==1 then bucket_dir_experimental_vcfs+"/annotated" else bucket_dir_experimental_vcfs,
                bucket_dir_ground_truth_vcfs = bucket_dir_ground_truth_vcfs,
                bucket_dir_allIndividuals_vcfs = bucket_dir_allIndividuals_vcfs,
                bucket_dir_matrices = bucket_dir_matrices,
                n_cpus = n_cpus_per_node
        }
    }
    output {
    }
}


# Collects all workpackage IDs into the rows of a file with format:
#
# caller svType svLength previousSvLength svFrequency previousSvFrequency
# repeatFraction previousRepeatFraction contextTypeStart contextTypeEnd
#
# where only one of {svType,svLength,svFrequency} can be non-null. Then, splits
# such a list into a given number of chunks.
#
task GetChunks {
    input {
        Array[String] callers
        Array[String] svTypes
        Array[Int] svLengths
        Array[Float] svFrequencies
        Array[Int] contextTypes
        Array[Float] repeatFractions
        Int n_chunks
        Array[Int]? force_sequentiality
    }
    command <<<
        set -euxo pipefail
        
        CALLERS=~{sep='-' callers}
        CALLERS=$(echo ${CALLERS} | tr '-' ' ')
        SV_TYPES=~{sep='-' svTypes}
        SV_TYPES=$(echo ${SV_TYPES} | tr '-' ' ')
        SV_LENGTHS=~{sep='-' svLengths}
        SV_LENGTHS=$(echo ${SV_LENGTHS} | tr '-' ' ')
        SV_FREQUENCIES=~{sep='-' svFrequencies}
        SV_FREQUENCIES=$(echo ${SV_FREQUENCIES} | tr '-' ' ')
        CONTEXT_TYPES=~{sep='-' contextTypes}
        CONTEXT_TYPES=$(echo ${CONTEXT_TYPES} | tr '-' ' ')
        REPEAT_FRACTIONS=~{sep='-' repeatFractions}
        REPEAT_FRACTIONS=$(echo ${REPEAT_FRACTIONS} | tr '-' ' ')
        for caller in ${CALLERS}; do
            # 1. SV TYPE
            for svType in ${SV_TYPES}; do
                # 1.1 Specific context
                for contextTypeStart in ${CONTEXT_TYPES}; do
                    for contextTypeEnd in ${CONTEXT_TYPES}; do
                        echo "${caller} ${svType} -1 -1 -1 -1 -1 -1 ${contextTypeStart} ${contextTypeEnd}" >> tmp.txt
                    done
                done
                previousRepeatFraction="0"
                for repeatFraction in ${REPEAT_FRACTIONS}; do
                    echo "${caller} ${svType} -1 -1 -1 -1 ${repeatFraction} ${previousRepeatFraction} -1 -1" >> tmp.txt
                    previousRepeatFraction=${repeatFraction}
                done
                # 1.2 Any context
                echo "${caller} ${svType} -1 -1 -1 -1 -1 -1 -1 -1" >> tmp.txt
            done
            # 2. SV LENGTH
            previousSvLength="0"
            for svLength in ${SV_LENGTHS}; do
                # 2.1 Specific context
                for contextTypeStart in ${CONTEXT_TYPES}; do
                    for contextTypeEnd in ${CONTEXT_TYPES}; do
                        echo "${caller} -1 ${svLength} ${previousSvLength} -1 -1 -1 -1 ${contextTypeStart} ${contextTypeEnd}" >> tmp.txt
                    done
                done
                previousRepeatFraction="0"
                for repeatFraction in ${REPEAT_FRACTIONS}; do
                    echo "${caller} -1 ${svLength} ${previousSvLength} -1 -1 ${repeatFraction} ${previousRepeatFraction} -1 -1" >> tmp.txt
                    previousRepeatFraction=${repeatFraction}
                done
                # 2.2 Any context
                echo "${caller} -1 ${svLength} ${previousSvLength} -1 -1 -1 -1 -1 -1" >> tmp.txt
                previousSvLength=${svLength}
            done
            # 3. SV FREQUENCY
            previousSvFrequency="0"
            for svFrequency in ${SV_FREQUENCIES}; do
                # 3.1 Specific context
                for contextTypeStart in ${CONTEXT_TYPES}; do
                    for contextTypeEnd in ${CONTEXT_TYPES}; do
                        echo "${caller} -1 -1 -1 ${svFrequency} ${previousSvFrequency} -1 -1 ${contextTypeStart} ${contextTypeEnd}" >> tmp.txt
                    done
                done
                previousRepeatFraction="0"
                for repeatFraction in ${REPEAT_FRACTIONS}; do
                    echo "${caller} -1 -1 -1 ${svFrequency} ${previousSvFrequency} ${repeatFraction} ${previousRepeatFraction} -1 -1" >> tmp.txt
                    previousRepeatFraction=${repeatFraction}
                done
                # 3.2 Any context
                echo "${caller} -1 -1 -1 ${svFrequency} ${previousSvFrequency} -1 -1 -1 -1" >> tmp.txt
                previousSvFrequency=${svFrequency}
            done
            # 4. NO (TYPE,LENGTH,FREQUENCY) CONSTRAINT
            # 4.1 Specific context
            for contextTypeStart in ${CONTEXT_TYPES}; do
                for contextTypeEnd in ${CONTEXT_TYPES}; do
                    echo "${caller} -1 -1 -1 -1 -1 -1 -1 ${contextTypeStart} ${contextTypeEnd}" >> tmp.txt
                done
            done
            previousRepeatFraction="0"
            for repeatFraction in ${REPEAT_FRACTIONS}; do
                echo "${caller} -1 -1 -1 -1 -1 ${repeatFraction} ${previousRepeatFraction} -1 -1" >> tmp.txt
                previousRepeatFraction=${repeatFraction}
            done
            # 4.2 Any context
            echo "${caller} -1 -1 -1 -1 -1 -1 -1 -1 -1" >> tmp.txt
        done
        shuf tmp.txt > workpackages.txt  # For better balancing
        N_LINES_PER_CHUNK="0"
        N_LINES=$(wc -l < workpackages.txt)
        if [ ~{n_chunks} -ge ${N_LINES} ]; then
            N_LINES_PER_CHUNK="1"
        else
            N_LINES_PER_CHUNK=$(( ${N_LINES} / ~{n_chunks} ))
        fi
        split -d -l ${N_LINES_PER_CHUNK} workpackages.txt chunk-
    >>>
    output {
        Array[File] chunks = glob("chunk-*")
    }
    runtime {
        docker: "ubuntu:latest"
    }
}


# Processes every workpackage in a given chunk, sequentially.
#
# Remark: the task checkpoints by using performance matrices in 
# <bucket_dir_matrices> and merged VCFs in <bucket_dir_allIndividuals_vcfs>.
#
task ProcessChunk {
    input {
        File chunk_file
        Array[Int] read_lengths
        Array[Int] coverages
        String reference_fa
        String reference_fai
        Float max_vcf_size
        Int n_individuals
        Int only_pass
        String bucket_dir_experimental_vcfs
        String bucket_dir_ground_truth_vcfs
        String bucket_dir_allIndividuals_vcfs
        String bucket_dir_matrices
        Int n_cpus
    }
    parameter_meta {
        chunk_file: "Every row is a workpackage ID created by $GetChunks$."
        reference_fa: "Address in a bucket"
        reference_fai: "Address in a bucket"
        n_individuals: "Total number of diploid individuals in the population. Used just to estimate space."
        max_vcf_size: "Max size of a single VCF file, in GB."
        only_pass: "Use only calls with FILTER=PASS (0/1)."
        n_cpus: "The program proceeds sequentially, but in each iteration it processes the VCF files from all individuals using all the cores of the machine."
    }
    Int ram_size_gb_pre = ceil(max_vcf_size*2*n_cpus)
    Int ram_size_gb = if 16 > ram_size_gb_pre then 16 else ram_size_gb_pre
    Int disk_size_gb = ceil(max_vcf_size*n_individuals*10) + ceil((size(reference_fa, "GB")))
    String docker_dir = "/simulation"
    String work_dir = "/cromwell_root/simulation"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        echo "Running <ProcessChunk> on ${N_THREADS} cores on the following node:"
        lscpu
        cat /proc/meminfo
        df -h
        
        while : ; do
            TEST=$(gsutil cp ~{reference_fa} reference.fa && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading file <~{reference_fa}>. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        while : ; do
            TEST=$(gsutil cp ~{reference_fai} reference.fa.fai && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading file <~{reference_fai}>. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        mkdir -p ground_truth_vcfs/
        while : ; do
            TEST=$(gsutil -m cp "~{bucket_dir_ground_truth_vcfs}/groundTruth_individual_*.vcf" ground_truth_vcfs/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading ground truth files from <~{bucket_dir_ground_truth_vcfs}>. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done        
        while : ; do
            TEST=$(gsutil -m cp "~{bucket_dir_ground_truth_vcfs}/groundTruth_joint.vcf" ground_truth_vcfs/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading ground truth file <~{bucket_dir_ground_truth_vcfs}/groundTruth_joint.vcf>. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        # Adding the SAMPLE column to the joint file
        bcftools view --threads 0 -h ground_truth_vcfs/groundTruth_joint.vcf > header.txt
        N_ROWS=$(wc -l < header.txt)
        head -n $((${N_ROWS} - 1)) header.txt > new.vcf
        rm -f header.txt
        echo "##INFO=<ID=HAPLOTYPES,Type=String,Description=\"haplotypes\">" >> new.vcf
        echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE" >> new.vcf
        bcftools view --threads 0 -H ground_truth_vcfs/groundTruth_joint.vcf | awk '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tGT\t.\n", $1, $2, $3, $4, $5, $6, $7, $8)}' >> new.vcf
        rm -f ground_truth_vcfs/groundTruth_joint.vcf
        cat new.vcf | sed '/contig=<ID=0>/d' | sed '/bcftools_/d' > ground_truth_vcfs/groundTruth_joint.vcf
        rm -f new.vcf
        ${TIME_COMMAND} bcftools sort --output-type z --output ground_truth_vcfs/groundTruth_joint.vcf.gz ground_truth_vcfs/groundTruth_joint.vcf
        tabix ground_truth_vcfs/groundTruth_joint.vcf.gz
        READ_LENGTHS=~{sep='-' read_lengths}
        READ_LENGTHS=$(echo ${READ_LENGTHS} | tr '-' ' ')
        COVERAGES=~{sep='-' coverages}
        COVERAGES=$(echo ${COVERAGES} | tr '-' ' ')
        while read LINE; do
            caller=$(echo ${LINE} | awk '{print $1}')
            svType=$(echo ${LINE} | awk '{print $2}')
            svLength=$(echo ${LINE} | awk '{print $3}')
            previousSvLength=$(echo ${LINE} | awk '{print $4}')
            svFrequency=$(echo ${LINE} | awk '{print $5}')
            previousSvFrequency=$(echo ${LINE} | awk '{print $6}')
            repeatFraction=$(echo ${LINE} | awk '{print $7}')
            previousRepeatFraction=$(echo ${LINE} | awk '{print $8}')
            contextTypeStart=$(echo ${LINE} | awk '{print $9}')
            contextTypeEnd=$(echo ${LINE} | awk '{print $10}')
            FILTER_STRING=""; FILTER_STRING_FREQUENCY=""; PREFIX=""
            if [ ${svType} != "-1" ]; then
                FILTER_STRING="SVTYPE=\"${svType}\""
                PREFIX="${caller}_svt${svType}"
            else
                FILTER_STRING=""
                PREFIX="${caller}"
            fi
            if [ ${svLength} != "-1" ]; then
                if [ ${#FILTER_STRING} -eq 0 ]; then
                    FILTER_STRING="(SVLEN>${previousSvLength} && SVLEN<=${svLength}) || (SVLEN>=-${svLength} && SVLEN<-${previousSvLength})"
                else
                    FILTER_STRING="${FILTER_STRING} && ((SVLEN>${previousSvLength} && SVLEN<=${svLength}) || (SVLEN>=-${svLength} && SVLEN<-${previousSvLength}))"
                fi
                PREFIX="${PREFIX}_svl${svLength}"
            else
                :  # NOP
            fi
            if [ ${svFrequency} != "-1" ]; then
                PREFIX="${PREFIX}_svf${svFrequency}"
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
            if [ ~{only_pass} -eq 1 ]; then
                if [ ${#FILTER_STRING} -eq 0 ]; then
                    FILTER_STRING="FILTER=\"PASS\""
                else
                    FILTER_STRING="FILTER=\"PASS\" && ${FILTER_STRING}"
                fi
            fi
            if [ ${svFrequency} != "-1" ]; then
                FILTER_STRING_FREQUENCY="AF>${previousSvFrequency} && AF<=${svFrequency}"
            fi
            if [ ${#FILTER_STRING} -eq 0 ]; then
                FILTER_STRING="+"
            else
                FILTER_STRING=$(echo ${FILTER_STRING} | tr ' ' '+')
            fi
            if [ ${#FILTER_STRING_FREQUENCY} -eq 0 ]; then
                FILTER_STRING_FREQUENCY="+"
            else
                FILTER_STRING_FREQUENCY=$(echo ${FILTER_STRING_FREQUENCY} | tr ' ' '+')
            fi
            TP_MATRIX="${PREFIX}_matrix_tp.txt"
            FP_MATRIX="${PREFIX}_matrix_fp.txt"
            FN_MATRIX="${PREFIX}_matrix_fn.txt"
            PRECISION_MATRIX="${PREFIX}_matrix_precision.txt"
            RECALL_MATRIX="${PREFIX}_matrix_recall.txt"
            F1_MATRIX="${PREFIX}_matrix_f1.txt"
            TP_MATRIX_MERGE="${PREFIX}_matrix_merge_tp.txt"
            FP_MATRIX_MERGE="${PREFIX}_matrix_merge_fp.txt"
            FN_MATRIX_MERGE="${PREFIX}_matrix_merge_fn.txt"
            PRECISION_MATRIX_MERGE="${PREFIX}_matrix_merge_precision.txt"
            RECALL_MATRIX_MERGE="${PREFIX}_matrix_merge_recall.txt"
            F1_MATRIX_MERGE="${PREFIX}_matrix_merge_f1.txt"
            TP_MATRIX_JOINT="${PREFIX}_matrix_joint_tp.txt"
            FP_MATRIX_JOINT="${PREFIX}_matrix_joint_fp.txt"
            FN_MATRIX_JOINT="${PREFIX}_matrix_joint_fn.txt"
            PRECISION_MATRIX_JOINT="${PREFIX}_matrix_joint_precision.txt"
            RECALL_MATRIX_JOINT="${PREFIX}_matrix_joint_recall.txt"
            F1_MATRIX_JOINT="${PREFIX}_matrix_joint_f1.txt"
            TEST1=$(gsutil -q stat ~{bucket_dir_matrices}/${TP_MATRIX} && echo 0 || echo 1)
            TEST2=$(gsutil -q stat ~{bucket_dir_matrices}/${FP_MATRIX} && echo 0 || echo 1)
            TEST3=$(gsutil -q stat ~{bucket_dir_matrices}/${FN_MATRIX} && echo 0 || echo 1)
            TEST4=$(gsutil -q stat ~{bucket_dir_matrices}/${PRECISION_MATRIX} && echo 0 || echo 1)
            TEST5=$(gsutil -q stat ~{bucket_dir_matrices}/${RECALL_MATRIX} && echo 0 || echo 1)
            TEST6=$(gsutil -q stat ~{bucket_dir_matrices}/${F1_MATRIX} && echo 0 || echo 1)
            TEST7=$(gsutil -q stat ~{bucket_dir_matrices}/${TP_MATRIX_MERGE} && echo 0 || echo 1)
            TEST8=$(gsutil -q stat ~{bucket_dir_matrices}/${FP_MATRIX_MERGE} && echo 0 || echo 1)
            TEST9=$(gsutil -q stat ~{bucket_dir_matrices}/${FN_MATRIX_MERGE} && echo 0 || echo 1)
            TEST10=$(gsutil -q stat ~{bucket_dir_matrices}/${PRECISION_MATRIX_MERGE} && echo 0 || echo 1)
            TEST11=$(gsutil -q stat ~{bucket_dir_matrices}/${RECALL_MATRIX_MERGE} && echo 0 || echo 1)
            TEST12=$(gsutil -q stat ~{bucket_dir_matrices}/${F1_MATRIX_MERGE} && echo 0 || echo 1)            
            if [ ${TEST1} -eq 0 -a ${TEST2} -eq 0 -a ${TEST3} -eq 0 -a ${TEST4} -eq 0 -a ${TEST5} -eq 0 -a ${TEST6} -eq 0 -a ${TEST7} -eq 0 -a ${TEST8} -eq 0 -a ${TEST9} -eq 0 -a ${TEST10} -eq 0 -a ${TEST11} -eq 0 -a ${TEST12} -eq 0 ]; then
                # We do not test the existence of JOINT matrices for simplicity.
                continue
            fi
            for readLength in ${READ_LENGTHS}; do
                for coverage in ${COVERAGES}; do
                    rm -rf experimental_vcfs/; mkdir -p experimental_vcfs/
                    while : ; do
                        TEST=$(gsutil -m cp "~{bucket_dir_experimental_vcfs}/${caller}_*_l${readLength}_c${coverage}_annotated.vcf" experimental_vcfs/ && echo 0 || echo 1)
                        if [ ${TEST} -eq 1 ]; then
                            echo "Error downloading experimental VCFs from <~{bucket_dir_experimental_vcfs}>. Trying again..."
                            sleep ${GSUTIL_DELAY_S}
                        else
                            break
                        fi
                    done
                    JOINT_CALLING_FILE="null"
                    TEST=$(gsutil -q stat "~{bucket_dir_experimental_vcfs}/joint_${caller}_l${readLength}_c${coverage}_annotated.vcf" && echo 0 || echo 1)
                    if [ ${TEST} -eq 0 ]; then
                        while : ; do
                            TEST=$(gsutil -m cp "~{bucket_dir_experimental_vcfs}/joint_${caller}_l${readLength}_c${coverage}_annotated.vcf" experimental_vcfs/ && echo 0 || echo 1)
                            if [ ${TEST} -eq 1 ]; then
                                echo "Error downloading file <~{bucket_dir_experimental_vcfs}/joint_${caller}_l${readLength}_c${coverage}_annotated.vcf>. Trying again..."
                                sleep ${GSUTIL_DELAY_S}
                            else
                                break
                            fi
                        done
                        JOINT_CALLING_FILE="experimental_vcfs/joint_${caller}_l${readLength}_c${coverage}_annotated.vcf"
                    fi
                    echo -n "${readLength},${coverage}," >> ${TP_MATRIX}
                    echo -n "${readLength},${coverage}," >> ${FP_MATRIX}
                    echo -n "${readLength},${coverage}," >> ${FN_MATRIX}
                    echo -n "${readLength},${coverage}," >> ${PRECISION_MATRIX}
                    echo -n "${readLength},${coverage}," >> ${RECALL_MATRIX}
                    echo -n "${readLength},${coverage}," >> ${F1_MATRIX}
                    echo -n "${readLength},${coverage}," >> ${TP_MATRIX_MERGE}
                    echo -n "${readLength},${coverage}," >> ${FP_MATRIX_MERGE}
                    echo -n "${readLength},${coverage}," >> ${FN_MATRIX_MERGE}
                    echo -n "${readLength},${coverage}," >> ${PRECISION_MATRIX_MERGE}
                    echo -n "${readLength},${coverage}," >> ${RECALL_MATRIX_MERGE}
                    echo -n "${readLength},${coverage}," >> ${F1_MATRIX_MERGE}
                    if [ ${JOINT_CALLING_FILE} != "null" ]; then
                        echo -n "${readLength},${coverage}," >> ${TP_MATRIX_JOINT}
                        echo -n "${readLength},${coverage}," >> ${FP_MATRIX_JOINT}
                        echo -n "${readLength},${coverage}," >> ${FN_MATRIX_JOINT}
                        echo -n "${readLength},${coverage}," >> ${PRECISION_MATRIX_JOINT}
                        echo -n "${readLength},${coverage}," >> ${RECALL_MATRIX_JOINT}
                        echo -n "${readLength},${coverage}," >> ${F1_MATRIX_JOINT}
                    fi
                    bash ~{docker_dir}/performance_matrices_impl.sh ${FILTER_STRING} ${FILTER_STRING_FREQUENCY} ${JOINT_CALLING_FILE} ~{reference_fa} ~{reference_fai} ${N_THREADS} ~{work_dir} "${PREFIX}_l${readLength}_c${coverage}" ${caller} ${readLength} ${coverage} ~{bucket_dir_allIndividuals_vcfs} ${TP_MATRIX} ${FP_MATRIX} ${FN_MATRIX} ${PRECISION_MATRIX} ${RECALL_MATRIX} ${F1_MATRIX}    ${TP_MATRIX_MERGE} ${FP_MATRIX_MERGE} ${FN_MATRIX_MERGE} ${PRECISION_MATRIX_MERGE} ${RECALL_MATRIX_MERGE} ${F1_MATRIX_MERGE}    ${TP_MATRIX_JOINT} ${FP_MATRIX_JOINT} ${FN_MATRIX_JOINT} ${PRECISION_MATRIX_JOINT} ${RECALL_MATRIX_JOINT} ${F1_MATRIX_JOINT}
                    echo "" >> ${TP_MATRIX}; echo "" >> ${FP_MATRIX}; echo "" >> ${FN_MATRIX}; echo "" >> ${PRECISION_MATRIX}; echo "" >> ${RECALL_MATRIX}; echo "" >> ${F1_MATRIX}
                    echo "" >> ${TP_MATRIX_MERGE}; echo "" >> ${FP_MATRIX_MERGE}; echo "" >> ${FN_MATRIX_MERGE}; echo "" >> ${PRECISION_MATRIX_MERGE}; echo "" >> ${RECALL_MATRIX_MERGE}; echo "" >> ${F1_MATRIX_MERGE}
                    if [ ${JOINT_CALLING_FILE} != "null" ]; then
                        echo "" >> ${TP_MATRIX_JOINT}; echo "" >> ${FP_MATRIX_JOINT}; echo "" >> ${FN_MATRIX_JOINT}; echo "" >> ${PRECISION_MATRIX_JOINT}; echo "" >> ${RECALL_MATRIX_JOINT}; echo "" >> ${F1_MATRIX_JOINT}
                    fi
                done
            done
            while : ; do
                TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp ${TP_MATRIX} ${FP_MATRIX} ${FN_MATRIX} ${PRECISION_MATRIX} ${RECALL_MATRIX} ${F1_MATRIX} ~{bucket_dir_matrices} && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error uploading per-individual matrices. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
            while : ; do
                TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp ${TP_MATRIX_MERGE} ${FP_MATRIX_MERGE} ${FN_MATRIX_MERGE} ${PRECISION_MATRIX_MERGE} ${RECALL_MATRIX_MERGE} ${F1_MATRIX_MERGE} ~{bucket_dir_matrices} && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error uploading merged matrices. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
            if [ -e ${TP_MATRIX_JOINT} ]; then
                while : ; do
                    TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp ${TP_MATRIX_JOINT} ${FP_MATRIX_JOINT} ${FN_MATRIX_JOINT} ${PRECISION_MATRIX_JOINT} ${RECALL_MATRIX_JOINT} ${F1_MATRIX_JOINT} ~{bucket_dir_matrices} && echo 0 || echo 1)
                    if [ ${TEST} -eq 1 ]; then
                        echo "Error uploading joint matrices. Trying again..."
                        sleep ${GSUTIL_DELAY_S}
                    else
                        break
                    fi
                done
            fi
        done < ~{chunk_file}
    >>>
    
    output {
    }
    runtime {
        docker: "fcunial/simulation"
        cpu: n_cpus
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
