version 1.0


# Runs PBSV and Sniffles 2 in joint mode, on each (caller, length, coverage)
# tuple in parallel.
#
# Resource analysis for 20x coverage of one haplotype over 300 individuals.
# Intel Xeon, 2.30GHz, 8 physical cores, 64 GB of RAM.
#
# TASK            | TIME   | RAM    | CORES | COMMENT
# sniffles2 joint | 2 m    | 1 GB   |   1   | 
# pbsv joint      | 14-23h | 18 GB  |  3.5  | coverage 4
#                 | 43 h   | 28 GB  |  3.7  | coverage 8
#                 | 50 h   | 33 GB  |  3.6  | coverage 12
#                 | ?? h   | ?? GB  |  ???  | coverage 16
#                 | ?? h   | ?? GB  |  ???  | coverage 20
#
# Remark: PBSV processes 1 Mbp chunk after another, in sequence rather than in
# parallel.
#
#
#
# Remark: for chr1 at 30x, for ~90 individuals, 4 threads: joint calling with
# PBSV takes ~40 hours and 18 GB of RAM (and it uses <=3 cores rather than
# the 4 assigned), and joint calling with Sniffles 2 takes just seconds.
#
workflow JointCalling {
    input {
        String bucket_dir
        Array[Int] coverages
        Array[Int] lengths
        Int n_individuals
        File reference_fa
        Float max_signature_file_size_pbsv = 0.020
        Float max_signature_file_size_sniffles2 = 0.003
        Int use_pbsv
        Int use_sniffles2
        Array[Int]? force_sequentiality
    }
    parameter_meta {
        bucket_dir: "URL of the root directory of the simulation."
        n_individuals: "Total number of diploid individuals in the population."
        max_signature_file_size: "In GB."
        use_pbsv: "1=yes, 0=no."
        force_sequentiality: "Fake input, just to make sure this workflow is executed after some previous step."
    }
    call CreateWorkpackages {
        input:
            coverages = coverages,
            lengths = lengths
    }
    scatter(description in CreateWorkpackages.tasks) {
        if (use_pbsv == 1) {
            call JointCallingPbsv { 
                input:
                    task_description = description,
                    bucket_dir = bucket_dir,
                    reference_fa = reference_fa,
                    n_individuals = n_individuals,
                    max_signature_file_size = max_signature_file_size_pbsv
            }
        }
        if (use_sniffles2 == 1) {
            call JointCallingSniffles2 { 
                input:
                    task_description = description,
                    bucket_dir = bucket_dir,
                    n_individuals = n_individuals,
                    max_signature_file_size = max_signature_file_size_sniffles2
            }
        }
    }
    output {
    }
}


task CreateWorkpackages {
    input {
        Array[Int] coverages
        Array[Int] lengths
    }
    command <<<
        set -euxo pipefail
        LENGTHS=~{sep='-' lengths}
        LENGTHS=$(echo ${LENGTHS} | tr '-' ' ')
        COVERAGES=~{sep='-' coverages}
        COVERAGES=$(echo ${COVERAGES} | tr '-' ' ')
        for LENGTH in ${LENGTHS}; do
            for COVERAGE in ${COVERAGES}; do
                echo "${LENGTH} ${COVERAGE}"
            done
        done
    >>>
    output {
        Array[String] tasks = read_lines(stdout())
    }
    runtime {
        docker: "ubuntu:latest"
    }
}


task JointCallingPbsv {
    input {
        String task_description
        String bucket_dir
        File reference_fa
        Int n_individuals
        Float max_signature_file_size
    }
    parameter_meta {
        task_description: "Format: <LENGTH COVERAGE>."
        n_individuals: "Total number of diploid individuals in the population."
        max_signature_file_size: "In GB."
    }
    # PBSV loads all input files in RAM
    Int ram_size_gb = ceil(n_individuals*max_signature_file_size*10) + ceil(size(reference_fa, "GB"))
    Int disk_size_gb = ram_size_gb
    String docker_dir = "/simulation"
    String work_dir = "/cromwell_root/simulation"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        GSUTIL_DELAY_S="600"
        
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        echo "Running <JointCallingPbsv> on ${N_THREADS} cores on the following node:"
        lscpu
        cat /proc/meminfo
        
        LENGTH=$(echo ~{task_description} | awk '{ print $1 }')
        COVERAGE=$(echo ~{task_description} | awk '{ print $2 }')
        while : ; do
            TEST=$(gsutil -m cp "~{bucket_dir}/signatures/pbsv_*_*_l${LENGTH}_c${COVERAGE}.svsig.gz" . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading files <~{bucket_dir}/signatures/pbsv_*_*_l${LENGTH}_c${COVERAGE}.svsig.gz>. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        find . -maxdepth 1 -name "*.svsig.gz" > filenames.fofn
        VCF_FILE="joint_pbsv_l${LENGTH}_c${COVERAGE}.vcf"
        ${TIME_COMMAND} pbsv call --log-level INFO -j ${N_THREADS} --ccs ~{reference_fa} filenames.fofn ${VCF_FILE}
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp ${VCF_FILE} ~{bucket_dir}/vcfs/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading file <${VCF_FILE}>. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
    >>>
    
    output {
    }
    runtime {
        docker: "fcunial/simulation"
        cpu: 6
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}


task JointCallingSniffles2 {
    input {
        String task_description
        String bucket_dir
        Int n_individuals
        Float max_signature_file_size
    }
    parameter_meta {
        task_description: "Format: <LENGTH COVERAGE>."
        n_individuals: "Total number of diploid individuals in the population."
        max_signature_file_size: "In GB."
    }
    Int disk_size_gb = ceil(n_individuals*max_signature_file_size) + 10
    String docker_dir = "/simulation"
    String work_dir = "/cromwell_root/simulation"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        GSUTIL_DELAY_S="600"
        
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        echo "Running <JointCallingSniffles2> on the following node:"
        lscpu
        cat /proc/meminfo
        
        LENGTH=$(echo ~{task_description} | awk '{ print $1 }')
        COVERAGE=$(echo ~{task_description} | awk '{ print $2 }')
        while : ; do
            TEST=$(gsutil -m cp "~{bucket_dir}/signatures/sniffles2_*_*_l${LENGTH}_c${COVERAGE}.snf" . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading files <~{bucket_dir}/signatures/sniffles2_*_*_l${LENGTH}_c${COVERAGE}.snf>. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        find . -maxdepth 1 -name "*.snf" > filenames.tsv
        VCF_FILE="joint_sniffles2_l${LENGTH}_c${COVERAGE}.vcf"
        ${TIME_COMMAND} sniffles --threads ${N_THREADS} --input filenames.tsv --vcf ${VCF_FILE}
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp ${VCF_FILE} ~{bucket_dir}/vcfs/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading file <${VCF_FILE}>. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
    >>>
    
    output {
    }
    runtime {
        docker: "fcunial/simulation"
        cpu: 1
        memory: "8 GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 3
    }
}