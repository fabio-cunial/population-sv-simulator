version 1.0


# Calls SVs (on the entire genome) from reads coming from children of trios,
# subsampled at several values of coverage.
#
# Remark: the workflow assumes that the number of children is low and it
# processes each child on a distinct machine.
#
workflow CallRealTrios {
    input {
        Int n_children
        File child_fa_addresses
        File reference_fa
        File reference_fai
        File reference_mmi
        File reference_tandem_repeats
        Array[Int] coverages
        String bucket_dir
        Int n_cpus_single_sample
        Int use_pbsv
        Int use_sniffles1
        Int use_sniffles2
        Int use_hifiasm
        Int use_pav
        Int use_paftools
        Int delete_bucket_dir
        Int keep_assemblies
    }
    parameter_meta {
        child_fa_addresses: "Text file with one line per address of a child reads fasta file."
        reference_mmi: "The minimap2 index of <reference_fa>."
        reference_tandem_repeats: "Used just for masking the reference by some callers."
        coverages: "Of one haplotype, sorted in increasing order."
        bucket_dir: "URL of a directory in a bucket. Stores both temporary files and output files."
        use_pbsv: "1=yes, 0=no."
        delete_bucket_dir: "Erases <bucket_dir> before running the simulation (1=yes, 0=no)."
        keep_assemblies: "Stores two assemblies per individual and per configuration in <bucket_dir>. 1=yes, 0=no."
    }
    if (delete_bucket_dir == 1) {
        call DeleteBucketDir {
            input:
                bucket_dir = bucket_dir
        }
    }
    Array[String] addresses = read_lines(child_fa_addresses)
    scatter(i in range(length(addresses))) {
        call ProcessChild {
            input:
                child_id = 2*i,
                child_reads = addresses[i],
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                reference_mmi = reference_mmi,
                reference_tandem_repeats = reference_tandem_repeats,
                coverages = coverages,
                bucket_dir = bucket_dir,
                n_cpus = n_cpus_single_sample,
                use_pbsv = use_pbsv,
                use_sniffles1 = use_sniffles1,
                use_sniffles2 = use_sniffles2,
                use_hifiasm = use_hifiasm,
                use_pav = use_pav,
                use_paftools = use_paftools,
                keep_assemblies = keep_assemblies,
                force_sequentiality = DeleteBucketDir.force_sequentiality
        }
    }
    output {
    }
}


task DeleteBucketDir {
    input { 
        String bucket_dir
    }    
    command <<<
        set -euxo pipefail
        GSUTIL_DELAY_S="60"
        
        TEST=$(gsutil -q stat ~{bucket_dir} && echo 0 || echo 1)
        if [ ${TEST} -eq 0 ]; then
            while ; do
                TEST=$(gsutil -m rm -rf ~{bucket_dir} && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]
                    echo "Error deleting <~{bucket_dir}>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        fi
        echo "1" > force_sequentiality.txt
    >>>
    output {
        Int force_sequentiality = read_int("force_sequentiality.txt")
    }
    runtime {
        docker: "fcunial/simulation"
    }
}


# Calls variants from the reads of a child, subsampled at several values of
# coverage.
#
# Remark: the task checkpoints by storing files in <bucket_dir>: (1) a
# checkpoint file for each individual, that contains lines <individual length
# coverage>, used to skip through the coverage loop; (2) the files generated
# by each step, used to skip that step if they already exist.
#
# File size analysis for ????x coverage of one haplotype:
# FILE     | PROGRAM   | SIZE 
# svsig.gz | pbsv      | 
# .snf     | sniffles2 | 
# .vcf     | any       | 
#
task ProcessChild {
    input {
        Int child_id
        File child_reads
        File reference_fa
        File reference_fai
        File reference_mmi
        File reference_tandem_repeats
        Array[Int] coverages
        String bucket_dir
        Int n_cpus
        Int use_pbsv
        Int use_sniffles1
        Int use_sniffles2
        Int use_hifiasm
        Int use_pav
        Int use_paftools
        Int keep_assemblies
        Int? force_sequentiality
    }
    parameter_meta {
        child_id: "Even number, starting from zero."
        reference_mmi: "The minimap2 index of <reference_fa>."
        reference_tandem_repeats: "Used just for masking the reference by some callers."
        coverages: "Of one haplotype, sorted in increasing order."
        n_cpus: "Used only for setting the runtime."
        use_pbsv: "1=yes, 0=no."
        keep_assemblies: "Stores two assemblies per individual and per configuration in the remote bucket. 1=yes, 0=no."
        force_sequentiality: "An artificial dependency used just to enforce sequential execution."
    }
    
    Int min_coverage = coverages[0]
    Int last = length(coverages) - 1
    Int max_coverage = coverages[last]
    Int ram_size_gb = 4 + ceil(size(child_reads, "GB") * 3)
    # *3: from loosely rounding up PAV; +4: to leave some free RAM to the OS.
    Int disk_size_gb = ram_size_gb*2 + ceil( size(reference_fa, "GB")*3 + size(reference_mmi, "GB") )
    String docker_dir = "/simulation"
    String work_dir = "/cromwell_root/simulation"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        GSUTIL_DELAY_S="600"
        FAKE_LENGTH="0"
        
        CHECKPOINT_FILE="checkpoint_c~{child_id}.txt"
        TEST=$(gsutil -q stat ~{bucket_dir}/checkpoints/${CHECKPOINT_FILE} && echo 0 || echo 1)
        if [ ${TEST} -eq 0 ]; then
            while : ; do
                TEST=$(gsutil cp ~{bucket_dir}/checkpoints/${CHECKPOINT_FILE} . && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading file <~{bucket_dir}/checkpoints/${CHECKPOINT_FILE}>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        else
            echo "0 0 0" > ${CHECKPOINT_FILE}
        fi
        COVERAGES=~{sep='-' coverages}
        bash ~{docker_dir}/reads2svs.sh ~{child_reads} ~{child_id} ${FAKE_LENGTH} ~{min_coverage} ~{max_coverage} ${COVERAGES} ~{reference_fa} ~{reference_fai} ~{reference_mmi} ~{reference_tandem_repeats} ${CHECKPOINT_FILE} ~{bucket_dir} ~{use_pbsv} ~{use_sniffles1} ~{use_sniffles2} ~{use_hifiasm} ~{use_pav} ~{use_paftools} ~{work_dir} ~{docker_dir} ~{keep_assemblies}
    >>>
    
    output {
    }
    runtime {
        docker: "fcunial/simulation"
        cpu: n_cpus
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 3
    }
}
