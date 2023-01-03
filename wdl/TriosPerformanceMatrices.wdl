version 1.0

import "AnnotateVCFs.wdl"


# The program compares each experimental VCF to its corresponding ground truth,
# collecting per-individual performance measures. This is performed separately
# for every tuple (caller, left weight, SV length).
#
workflow TriosPerformanceMatrices {
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
    }
    parameter_meta {
        svFrequencies: "They should be fractions $X/(2*n_individuals)$, where $X$ is an integer number of haplotypes."
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






task ChildPerformanceMatrices {
    input {
        String child_id
        String bucket_dir
        Array[String] callers
        Array[Float] left_weights
        Array[Int] sv_lengths
        Int only_pass
    }
    parameter_meta {
    }
    
    Int ram_size_gb_pre = ceil(max_vcf_size*2*n_cpus)
    Int ram_size_gb = if 16 > ram_size_gb_pre then 16 else ram_size_gb_pre
    Int disk_size_gb = ceil(max_vcf_size*n_individuals*20) + ceil((size(reference_fa, "GB")))
    String docker_dir = "/simulation"
    String work_dir = "/cromwell_root/simulation"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        CALLERS=~{sep='-' callers}
        CALLERS=$(echo ${CALLERS} | tr '-' ' ')
        WEIGHTS=~{sep='-' left_weights}
        WEIGHTS=$(echo ${WEIGHTS} | tr '-' ' ')
        SV_LENGTHS=~{sep='-' sv_lengths}
        SV_LENGTHS=$(echo ${SV_LENGTHS} | tr '-' ' ')
        for caller in ${CALLERS}; do
            while : ; do
                TEST=$(gsutil -m cp "~{bucket_dir}/~{child_id}/${caller}_truth.vcf.gz" . && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading file <~{bucket_dir}/~{child_id}/${caller}_truth.vcf.gz>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
            TP_MATRIX="${caller}_matrix_tp.txt"
            FP_MATRIX="${caller}_matrix_fp.txt"
            FN_MATRIX="${caller}_matrix_fn.txt"
            PRECISION_MATRIX="${caller}_matrix_precision.txt"
            RECALL_MATRIX="${caller}_matrix_recall.txt"
            F1_MATRIX="${caller}_matrix_f1.txt"
            touch ${TP_MATRIX} ${FP_MATRIX} ${FN_MATRIX} ${PRECISION_MATRIX} ${RECALL_MATRIX} ${F1_MATRIX}
            for weight in ${WEIGHTS}; do
                TEST=$(gsutil -q stat "~{bucket_dir}/~{child_id}/reads_${weight}/${caller}_~{child_id}.vcf" . && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    continue
                fi
                while : ; do
                    TEST=$(gsutil -m cp "~{bucket_dir}/~{child_id}/reads_${weight}/${caller}_~{child_id}.vcf" . && echo 0 || echo 1)
                    if [ ${TEST} -eq 1 ]; then
                        echo "Error downloading file <~{bucket_dir}/~{child_id}/reads_${weight}/${caller}_~{child_id}.vcf>. Trying again..."
                        sleep ${GSUTIL_DELAY_S}
                    else
                        break
                    fi
                done
                for sv_length in ${SV_LENGTHS}; do
                    echo -n "${weight},${sv_length}," >> ${TP_MATRIX}
                    echo -n "${weight},${sv_length}," >> ${FP_MATRIX}
                    echo -n "${weight},${sv_length}," >> ${FN_MATRIX}
                    echo -n "${weight},${sv_length}," >> ${PRECISION_MATRIX}
                    echo -n "${weight},${sv_length}," >> ${RECALL_MATRIX}
                    echo -n "${weight},${sv_length}," >> ${F1_MATRIX}
                    
                    
                    ----------------->
                    
                    
                    
                    echo "" >> ${TP_MATRIX}; echo "" >> ${FP_MATRIX}; echo "" >> ${FN_MATRIX}; echo "" >> ${PRECISION_MATRIX}; echo "" >> ${RECALL_MATRIX}; echo "" >> ${F1_MATRIX}
                done
            done
            while : ; do
                TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp ${TP_MATRIX} ${FP_MATRIX} ${FN_MATRIX} ${PRECISION_MATRIX} ${RECALL_MATRIX} ${F1_MATRIX} ~{bucket_dir}/~{child_id}/ && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error uploading per-individual matrices. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        done
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
