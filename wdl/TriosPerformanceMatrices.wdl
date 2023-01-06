version 1.0

import "AnnotateVCFs.wdl"


# The program compares the VCF file for each tuple (trio child, value, caller)
# to its corresponding "truth" file created by a <TrioChildCreateTruthVCFs*>
# script. SVs of different lengths are analyzed separately. The results are
# collected in a matrix for every child and caller, whose lines follow the
# format <value, svLength, metric>. By <value> we mean a left weight in
# Question 1 and a left coverage in Question 2.
#
workflow TriosPerformanceMatrices {
    input {
        File children_ids
        String bucket_dir
        Int answer_question1
        Array[Float] question1_left_weights
        Int answer_question2
        Array[Int] question2_left_coverages
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
                answer_question1 = answer_question1,
                question1_left_weights = question1_left_weights,
                answer_question2 = answer_question2,
                question2_left_coverages = question2_left_coverages,
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
        Int answer_question1
        Array[Float] question1_left_weights
        Int answer_question2
        Array[Int] question2_left_coverages
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
        
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        TRUVARI_BENCH_FLAGS=" "  # Using the default settings for now
        
        CALLERS=~{sep='-' callers}
        CALLERS=$(echo ${CALLERS} | tr '-' ' ')
        SV_LENGTHS=~{sep='-' sv_lengths}
        SV_LENGTHS=$(echo ${SV_LENGTHS} | tr '-' ' ')
        
        if [ ~{answer_question1} -eq 1 ]; then
            VALUES=~{sep='-' question1_left_weights}
            VALUES=$(echo ${VALUES} | tr '-' ' ')
            bash ~{docker_dir}/triosPerformanceMatrices_impl.sh ~{bucket_dir} ~{child_id} ${CALLERS} ${VALUES} ${SV_LENGTHS} full_coverage w ~{only_pass} ~{reference_fa} ~{work_dir}
        fi
        if [ ~{answer_question2} -eq 1 ]; then
            VALUES=~{sep='-' question2_left_coverages}
            VALUES=$(echo ${VALUES} | tr '-' ' ')
            bash ~{docker_dir}/triosPerformanceMatrices_impl.sh ~{bucket_dir} ~{child_id} ${CALLERS} ${VALUES} ${SV_LENGTHS} long_coverage c ~{only_pass} ~{reference_fa} ~{work_dir}
        fi
    >>>
    
    output {
    }
    runtime {
        docker: "fcunial/simulation"
        cpu: 1
        memory: "8GB"  # Arbitrary
        disks: "local-disk 50 HDD"  # Arbitrary
        preemptible: 3
    }
}
