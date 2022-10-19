version 1.0CreateDescriptions

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
    }
    parameter_meta {
        n_haplotypes: "Total number of haplotypes in the simulated population."
    }
    call CreateDescriptions {
        input:
            n_haplotypes = n_haplotypes,
            coverages = coverages,
            lengths = lengths
    }
    scatter(description in CreateDescriptions.tasks) {
        call RetrieveAssembly {
            input:
                filename = description,
                reference_fa = reference_fa
        }
        call pav.pav {
            input:
                ref = reference_fa,
                refFai = reference_fai,
                hapOne = bucket_dir + "/assemblies/assembly_" + RetrieveAssembly.id1 + "_" + RetrieveAssembly.id2 + "_" + RetrieveAssembly.length + "_" + RetrieveAssembly.coverage + "_h1.fa",
                hapTwo = bucket_dir + "/assemblies/assembly_" + RetrieveAssembly.id1 + "_" + RetrieveAssembly.id2 + "_" + RetrieveAssembly.length + "_" + RetrieveAssembly.coverage + "_h2.fa",
                sample = "assembly_" + RetrieveAssembly.id1 + "_" + RetrieveAssembly.id2 + "_" + RetrieveAssembly.length + "_" + RetrieveAssembly.coverage + ".fa",
                config = RetrieveAssembly.config
        }
        call ArchiveVCF {
            input:
                vcf = pav.vcf,
                bucket_dir = bucket_dir,
                id1 = RetrieveAssembly.id1,
                id2 = RetrieveAssembly.id2,
                length = RetrieveAssembly.length,
                coverage = RetrieveAssembly.coverage
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


# Remark: the task creates also the <config.json> file for PAV.
task ReadDescription {
    input {
        String description
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


task ArchiveVCF {
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
        FILENAME="pav_~{id1}_~{id2}_~{length}_~{coverage}.vcf"
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
