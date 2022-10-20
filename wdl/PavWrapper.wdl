version 1.0

import "../pav-wdl/pav.wdl"


# Calls the existing PAV WDL for every triplet (diploid individual, length,
# coverage) in parallel.
workflow PavWrapper {
    input {
        String bucket_dir
        Int n_haplotypes
        Array[Int] coverages
        Array[Int] lengths
        File reference_fa
        File reference_fai
        Array[Int] force_sequentiality
    }
    parameter_meta {
        n_haplotypes: "Total number of haplotypes in the simulated population."
        force_sequentiality: "Fake input, just to make sure that this workflow is executed after some previous step has completed."
    }
    call CreateDescriptions {
        input:
            n_haplotypes = n_haplotypes,
            coverages = coverages,
            lengths = lengths
    }
    scatter(description in CreateDescriptions.descriptions) {
        call ReadDescription {
            input:
                description = description,
                reference_fa = reference_fa
        }
        call pav.pav {
            input:
                ref = reference_fa,
                refFai = reference_fai,
                hapOne = bucket_dir + "/assemblies/assembly_i" + ReadDescription.id1 + "_i" + ReadDescription.id2 + "_l" + ReadDescription.length + "_c" + ReadDescription.coverage + "_h1.fa",
                hapTwo = bucket_dir + "/assemblies/assembly_i" + ReadDescription.id1 + "_i" + ReadDescription.id2 + "_l" + ReadDescription.length + "_c" + ReadDescription.coverage + "_h2.fa",
                sample = "assembly_i" + ReadDescription.id1 + "_i" + ReadDescription.id2 + "_l" + ReadDescription.length + "_c" + ReadDescription.coverage,
                config = ReadDescription.config
        }
        call StoreVCF {
            input:
                vcf = pav.vcf,
                bucket_dir = bucket_dir,
                id1 = ReadDescription.id1,
                id2 = ReadDescription.id2,
                length = ReadDescription.length,
                coverage = ReadDescription.coverage
        }
    }
    output {
    }
}


task CreateDescriptions {
    input {
        Int n_haplotypes
        Array[Int] coverages
        Array[Int] lengths
    }
    command <<<
        set -euxo pipefail
        COVERAGES=~{sep='-' coverages}
        COVERAGES=$(echo ${COVERAGES} | tr '-' ' ')
        LENGTHS=~{sep='-' lengths}
        LENGTHS=$(echo ${LENGTHS} | tr '-' ' ')
        for ID1 in $(seq 0 2 $(( ~{n_haplotypes} - 1 ))); do
            ID2=$(( ${ID1} + 1 ))
            for LENGTH in ${LENGTHS}; do
                for COVERAGE in ${COVERAGES}; do
                    echo "${ID1}-${ID2}-${LENGTH}-${COVERAGE}"
                done
            done
        done
    >>>
    output {
        Array[String] descriptions = read_lines(stdout())
    }
    runtime {
        docker: "ubuntu:latest"
    }
}


task ReadDescription {
    input {
        File description
        String reference_fa
    }
    parameter_meta {
        reference_fa: "Path of the reference in the simulation bucket"
    }
    command <<<
        set -euxo pipefail
        IFS='-' read -ra TOKENS <<< ~{description}
        echo ${TOKENS[0]} > id1.txt
        echo ${TOKENS[1]} > id2.txt
        echo ${TOKENS[2]} > length.txt
        echo ${TOKENS[3]} > coverage.txt
        echo "{" >> config.json
        echo "\"reference\": \"~{reference_fa}\"," >> config.json
        echo "\"asm_pattern\": \"\{asm_name\}_\{hap\}.fa\"," >> config.json
        echo "\"aligner\": \"minimap2\"" >> config.json
        echo "}" >> config.json
    >>>
    output {
        Int id1 = read_int("id1.txt")
        Int id2 = read_int("id2.txt")
        Int length = read_int("length.txt")
        Int coverage = read_int("coverage.txt")
        File config = "config.json"
    }
    runtime {
        docker: "ubuntu:latest"
    }
}


task StoreVCF {
    input {
        File vcf
        String bucket_dir
        Int id1
        Int id2
        Int length
        Int coverage
    }
    command <<<
        set -euxo pipefail
        FILENAME="pav_i~{id1}_i~{id2}_l~{length}_c~{coverage}.vcf"
        cp ~{vcf} ${FILENAME}
        gsutil cp ${FILENAME} ~{bucket_dir}/vcfs/
    >>>
    output {
    }
    runtime {
        docker: "fcunial/simulation"
        disks: "local-disk " + (ceil(size(vcf, "GB"))*2) + " HDD"
    }
}
