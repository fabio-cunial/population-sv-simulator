version 1.0

import "TrioChildCreateTruthVCFs2.wdl"


# Runs several callers on every child of a trio and on its parents, using, for
# each individual, a readset whose length histogram is identical to the full
# readset, but all the modes except the one centered on the rightmost peak are
# removed. Finally, keeps only the variants of a child that occur also in some
# parent.
#
workflow TriosCreateTruthVCFs2 {
    input {
        String bucket_dir
        File children_ids
        Int bin_length
        Int max_read_length
        File reference_fa
        File reference_fai
        File reference_tandem_repeats
        Int n_cpus
        Int use_pbsv
        Int use_sniffles1
        Int use_sniffles2
        Int use_hifiasm
        Int use_pav
        Int use_paftools
        Int keep_assemblies
    }
    parameter_meta {
        children_ids: "The trio children to be processed"
        n_cpus: "For the bottleneck part of creating VCFs from full-coverage read sets"
    }
    
    scatter(child_id in read_lines(children_ids)) {
        call TrioChildCreateTruthVCFs2.TrioChildCreateTruthVCFs2 {
            input:
                child_id = child_id,
                bucket_dir = bucket_dir,
                bin_length = bin_length,
                max_read_length = max_read_length,
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                reference_tandem_repeats = reference_tandem_repeats,
                n_cpus = n_cpus,
                use_pbsv = use_pbsv,
                use_sniffles1 = use_sniffles1,
                use_sniffles2 = use_sniffles2,
                use_hifiasm = use_hifiasm,
                use_pav = use_pav,
                use_paftools = use_paftools,
                keep_assemblies = keep_assemblies
        }
    }
    
    output {
    }
}
