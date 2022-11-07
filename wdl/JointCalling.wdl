version 1.0


# Runs PBSV and Sniffles 2 in joint mode, on each (caller, length, coverage)
# tuple in parallel.
#
workflow JointCalling {
    input {
        String bucket_dir
        Array[Int] coverages
        Array[Int] lengths
        Int n_individuals
        File reference_fa
        Float max_signature_file_size = 0.020
        Int n_cpus
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
            lengths = lengths,
            use_pbsv = use_pbsv,
            use_sniffles2 = use_sniffles2
    }
    scatter(description in CreateWorkpackages.tasks) {
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


task CreateWorkpackages {
    input {
        Array[Int] coverages
        Array[Int] lengths
        Int use_pbsv
        Int use_sniffles2
    }
    command <<<
        set -euxo pipefail
        LENGTHS=~{sep='-' lengths}
        LENGTHS=$(echo ${LENGTHS} | tr '-' ' ')
        COVERAGES=~{sep='-' coverages}
        COVERAGES=$(echo ${COVERAGES} | tr '-' ' ')
        for LENGTH in ${LENGTHS}; do
            for COVERAGE in ${COVERAGES}; do
                if [ ~{use_pbsv} -eq 1 ]; then
                    echo "1 ${LENGTH} ${COVERAGE}"
                fi
                if [ ~{use_sniffles2} -eq 1 ]; then
                    echo "2 ${LENGTH} ${COVERAGE}"
                fi
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


# Remark: for chr1 at 30x, for ~90 individuals (each <.svsig.gz> file takes ~15
# MB, and each <.snf> file takes ~2.5 MB) and for 4 threads: joint calling with
# PBSV takes ~40 hours and 18 GB of RAM (and it uses <=3 cores rather than
# the 4 assigned), and joint calling with Sniffles 2 takes just seconds.
#
#
# Resource analysis for 20x coverage of one haplotype. Intel ???, ???GHz, ???
# physical cores.
#
# TASK                    | TIME   | RAM    | CORES | COMMENT
# pbsv joint              | ? s    | ??? MB |   ?   |
# sniffles2 joint         | ?? s   | ? GB   |   ?   |
#
task JointCallingImpl {
    input {
        String task_description
        String bucket_dir
        File reference_fa
        Int n_individuals
        Float max_signature_file_size
        Int n_cpus
    }
    parameter_meta {
        task_description: "Format: <CALLER LENGTH COVERAGE>, where CALLER is 1 (PBSV) or 2 (Sniffles2)."
        n_individuals: "Total number of diploid individuals in the population."
        max_signature_file_size: "In GB."
    }
    # Both PBSV and Sniffles 2 seem to load all input files in RAM
    Int ram_size_gb = ceil(n_individuals*max_signature_file_size)*2 + ceil(size(reference_fa, "GB"))
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
        echo "Running <JointCallingImpl> on ${N_THREADS} cores on the following node:"
        lscpu
        cat /proc/meminfo
        
        CALLER=$(echo ~{task_description} | awk '{ print $1 }')
        LENGTH=$(echo ~{task_description} | awk '{ print $2 }')
        COVERAGE=$(echo ~{task_description} | awk '{ print $3 }')
        if [ ${CALLER} -eq 1 ]; then
            # PBSV
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
            ${TIME_COMMAND} pbsv call -j ${N_THREADS} --ccs ~{reference_fa} filenames.fofn ${VCF_FILE}
            while : ; do
                TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp ${VCF_FILE} ~{bucket_dir}/vcfs/ && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error uploading file <${VCF_FILE}>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        else
            # Sniffles 2
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
        fi
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
