version 1.0


# The program compares the VCF file for each tuple (trio child, value, caller)
# to its corresponding "truth" file created by a <TrioChildCreateTruthVCFs*>
# script. SVs of different lengths are analyzed separately. The results are
# collected in a matrix for every child and caller, whose lines follow the
# format <coverage, svLength, metric>.
#
workflow TriosCoverageEffectPerformanceMatrices {
    input {
        File children_ids
        String bucket_dir
        Array[Float] right_coverages
        Array[String] callers
        Array[Int] sv_lengths
        Int only_pass
        File reference_fa
    }
    parameter_meta {
        children_ids: "The trio children to be processed"
        only_pass: "(0/1) Use only variants with FILTER=PASS."
    }    

    scatter(child_id in read_lines(children_ids)) {
        call ChildPerformanceMatrices { 
            input:
                child_id = child_id,
                bucket_dir = bucket_dir,
                right_coverages = right_coverages,
                callers = callers,
                sv_lengths = sv_lengths,
                only_pass = only_pass,
                reference_fa = reference_fa
        }
    }
    output {
    }
}


task ChildPerformanceMatrices {
    input {
        String child_id
        String bucket_dir
        Array[Float] right_coverages
        Array[String] callers
        Array[Int] sv_lengths
        Int only_pass
        File reference_fa
    }
    parameter_meta {
        only_pass: "Use only calls with FILTER=PASS (0/1)."
    }
    
    String docker_dir = "/simulation"
    String work_dir = "/cromwell_root/simulation"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        CALLERS=~{sep='-' callers}
        SV_LENGTHS=~{sep='-' sv_lengths}
        VALUES=~{sep='-' right_coverages}
        bash ~{docker_dir}/triosPerformanceMatrices_impl.sh ~{bucket_dir}/~{child_id}/coverage_effect c ~{bucket_dir}/~{child_id} long_coverage ~{child_id} ${CALLERS} ${VALUES} ${SV_LENGTHS} ~{only_pass} ~{reference_fa} ~{work_dir}
    >>>
    
    output {
    }
    runtime {
        docker: "fcunial/simulation"
        cpu: 10  # Arbitrary
        memory: "16GB"  # Arbitrary
        disks: "local-disk 50 HDD"  # Arbitrary
        preemptible: 0  # The double loop to create a single matrix can take hours and there is no checkpointing in the middle
    }
}
