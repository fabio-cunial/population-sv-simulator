version 1.0


# For every trio child and caller, the workflow creates a matrix whose columns
# are SV length bins, and whose rows are: 
# - the ground truth VCF; 
# - each one of the 3 VCFs used to compute the ground truth VCF;
# - each VCF obtained by creating a simulated read length distribution.
# A cell of the matrix stores the total number of SVs in a given VCF and length
# bin (SVs are not decomposed by type).
#
workflow TriosCreateSVLengthHistograms {
    input {
        String bucket_dir
        File children_ids
        Array[String] callers
        Array[Int] sv_lengths
        Int answer_question1
        Array[Float] question1_left_weights
        Int answer_question2
        Array[Float] question2_left_coverages
        Int only_pass
        Int n_cpus
    }
    parameter_meta {
        n_cpus: "In a single compute node"
    }
    
    scatter(child in read_lines(children_ids)) {
        call TrioChildCreateSVLengthHistograms {
            input:
                bucket_dir = bucket_dir,
                child_id = child,
                callers = callers,
                sv_lengths = sv_lengths,
                answer_question1 = answer_question1,
                question1_left_weights = question1_left_weights,
                answer_question2 = answer_question2,
                question2_left_coverages = question2_left_coverages,
                only_pass = only_pass,
                n_cpus = n_cpus
        }
    }        
    
    output {
    }
}


task TrioChildCreateSVLengthHistograms {
    input {
        String bucket_dir
        String child_id
        Array[String] callers
        Array[Int] sv_lengths
        Int answer_question1
        Array[Float] question1_left_weights
        Int answer_question2
        Array[Float] question2_left_coverages
        Int only_pass
        Int n_cpus
    }
    
    String docker_dir = "/simulation"
    String work_dir = "/cromwell_root/simulation"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        CALLERS=~{sep='-' callers}
        SV_LENGTHS=~{sep='-' sv_lengths}
        if [ ~{answer_question1} -eq 1 ]; then
            TRUTH_VCF_PREFIX="full_coverage"
            PRECURSOR_VCF_PREFIX="full_coverage"
            MEASURED_CHARACTER_CODE="w"
            VALUES=~{sep='-' question1_left_weights}
            bash ~{docker_dir}/triosGetSVLengthHistograms_impl.sh ~{bucket_dir} ~{child_id} ${TRUTH_VCF_PREFIX} ${PRECURSOR_VCF_PREFIX} ${MEASURED_CHARACTER_CODE} ${CALLERS} ${VALUES} ${SV_LENGTHS} ~{only_pass} ~{work_dir}
        fi
        if [ ~{answer_question2} -eq 1 ]; then
            TRUTH_VCF_PREFIX="long_coverage"
            PRECURSOR_VCF_PREFIX="long_coverage"
            MEASURED_CHARACTER_CODE="c"
            VALUES=~{sep='-' question2_left_coverages}
            bash ~{docker_dir}/triosGetSVLengthHistograms_impl.sh ~{bucket_dir} ~{child_id} ${TRUTH_VCF_PREFIX} ${PRECURSOR_VCF_PREFIX} ${MEASURED_CHARACTER_CODE} ${CALLERS} ${VALUES} ${SV_LENGTHS} ~{only_pass} ~{work_dir}
        fi
    >>>
    
    output {
    }
    runtime {
        docker: "fcunial/simulation"
        cpu: n_cpus
        memory: "16GB"  # Arbitrary
        disks: "local-disk 32 HDD"  # Arbitrary
        bootDiskSizeGb: 15
        preemptible: 0
    }
}