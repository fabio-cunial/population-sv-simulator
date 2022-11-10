version 1.0


# Adds repeat annotations to every VCF file in a bucket.
#
# Remark: this workflow can be run while the VCF files are being created by the
# simulation workflow.
#
workflow AnnotateVCFs {
    input {
        String bucket_dir
        Float max_vcf_size
        Int min_distance
        Int max_distance
        File reference_fa
        File repeat_intervals
        File segdups
        Int n_nodes
    }
    parameter_meta {
        max_vcf_size: "Max size of a single VCF file, in GB."
        min_distance: "Defined in <AnnotateVCF.java>."
        max_distance: "Defined in <AnnotateVCF.java>."
        n_nodes: "Number of machines over which to distribute the computation (which is sequential on each machine)."
    }
    call GetChunks {
        input:
            bucket_dir = bucket_dir,
            n_chunks = n_nodes,
            min_files_per_chunk = 2000  # Fast enough in practice
    }
    scatter(chunk_file in GetChunks.chunks) {
        call ProcessChunk {
            input:
                chunk_file = chunk_file,
                n_files = GetChunks.n_files_per_chunk,
                max_vcf_size = max_vcf_size,
                min_distance = min_distance,
                max_distance = max_distance,
                reference_fa = reference_fa,
                repeat_intervals = repeat_intervals,
                segdups = segdups,
                bucket_dir = bucket_dir
        }
    }
    output {
        Array[Int] force_sequentiality = ProcessChunk.force_sequentiality
    }
}


task GetChunks {
    input {
        String bucket_dir
        Int n_chunks
        Int min_files_per_chunk
    }
    parameter_meta {
        min_files_per_chunk: "Guarantees at least this number of files in a chunk (if there are enough files)."
    }
    command <<<
        set -euxo pipefail
        gsutil ls ~{bucket_dir} > tmp.txt
        shuf tmp.txt > workpackages.txt  # For better balancing
        N_FILES=$(wc -l < workpackages.txt)
        if [ ${N_FILES} -le ~{min_files_per_chunk} ]; then
            cp workpackages.txt chunk-00
            echo ${N_FILES} > n_files_per_chunk.txt
        else
            N_FILES_PER_CHUNK=$(( ${N_FILES} / ~{n_chunks} ))
            if [ ${N_FILES_PER_CHUNK} -lt ~{min_files_per_chunk} ]; then
                N_FILES_PER_CHUNK=~{min_files_per_chunk}
            fi
            split -d -l ${N_FILES_PER_CHUNK} workpackages.txt chunk-
            echo ${N_FILES_PER_CHUNK} > n_files_per_chunk.txt
        fi
    >>>
    output {
        Array[File] chunks = glob("chunk-*")
        Int n_files_per_chunk = read_int("n_files_per_chunk.txt")
    }
    runtime {
        docker: "fcunial/simulation"
    }
}


task ProcessChunk {
    input {
        File chunk_file
        Int n_files
        Float max_vcf_size
        Int min_distance
        Int max_distance
        File reference_fa
        File repeat_intervals
        File segdups
        String bucket_dir
    }
    parameter_meta {
        max_vcf_size: "Max size of a single VCF file, in GB."
        bucket_dir: "Annotated VCFs are stored in subdirectory <annotated> of this directory."
    }
    Int disk_size_gb = ceil(n_files*max_vcf_size*2) + ceil((size(reference_fa, "GB")))
    String docker_dir = "/simulation"
    String work_dir = "/cromwell_root/simulation"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        TIME_COMMAND="/usr/bin/time --verbose"
        echo "Running <ProcessChunk> on the following node:"
        lscpu
        cat /proc/meminfo
        
        MEMORY=$(grep 'MemAvailable:' /proc/meminfo | awk '{print $2}' )
        MEMORY=$(( ${MEMORY} / 1000000 - 1 ))  # We leave 1 GB to the OS
        cat ~{chunk_file} | gsutil -m cp -I .
        find . -maxdepth 1 -name "*.vcf" > list.txt
        ${TIME_COMMAND} java -Xmx${MEMORY}G -cp ~{docker_dir} AnnotateVCF list.txt ~{min_distance} ~{max_distance} ~{reference_fa} ~{repeat_intervals} ~{segdups}
        find . -maxdepth 1 -name "*_annotated.vcf" > list.txt
        gsutil cp list.txt ~{bucket_dir}/annotated/list.txt
        cat list.txt | gsutil ${GSUTIL_UPLOAD_THRESHOLD} -m cp -I ~{bucket_dir}/annotated/
        echo "1" > force_sequentiality.txt
    >>>
    
    output {
        Int force_sequentiality = read_int("~{work_dir}/force_sequentiality.txt")
    }
    runtime {
        docker: "fcunial/simulation"
        cpu: 1  # The program is sequential
        memory: "24 GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
