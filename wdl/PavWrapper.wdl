version 1.0

import "../pav-wdl/pav.wdl"


# Calls the existing PAV WDL for every triplet (diploid individual, length,
# coverage) in parallel.
workflow PavWrapper {
    input {
        String bucket_address
        Int n_haplotypes
        Array[Int] coverages
        Array[Int] lengths
        File reference_fa
        File reference_fai
    }
    parameter_meta {
        n_haplotypes: "Total number of haplotypes in the simulated population."
    }
    call CreateWorkpackages {
        input:
            n_haplotypes = n_haplotypes,
            coverages = coverages,
            lengths = lengths
    }
    scatter(description in CreateWorkpackages.tasks) {
        call RetrieveAssembly {
            input:
                bucket_address = bucket_address,
                filename = description,
                haplotype_size = size(reference_fa, "GB"),
                reference_fa = reference_fa
        }
        call pav.pav {
            input:
                ref = reference_fa,
                refFai = reference_fai,
                hapOne = RetrieveAssembly.h1,
                hapTwo = RetrieveAssembly.h2,
                sample = "assembly_" + RetrieveAssembly.id1 + "_" + RetrieveAssembly.id2 + "_" + RetrieveAssembly.length + "_" + RetrieveAssembly.coverage + "_h1.fa",
                config = RetrieveAssembly.config
        }
        call ArchiveVCF {
            input:
                vcf = pav.vcf,
                bucket_address = bucket_address,
                id1 = RetrieveAssembly.id1,
                id2 = RetrieveAssembly.id2,
                length = RetrieveAssembly.length,
                coverage = RetrieveAssembly.coverage
        }
    }
    output {
    }
}


task CreateWorkpackages {
    input {
        Int n_haplotypes
        Array[Int] coverages
        Array[Int] lengths
    }
    command <<<
        set -euxo pipefail
        for ID1 in $(seq 0 2 $(( ~{n_haplotypes} - 1 ))); do
            ID2=$(( ${ID1} + 1 ))
            for LENGTH in ~{sep=' ' lengths}; do
                for COVERAGE in ~{sep=' ' coverages}; do
                    echo "assembly_i${ID1}_i${ID2}_l${LENGTH}_c${COVERAGE}.tar"
                done
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


# Remark: the task creates also the <config.json> file for PAV.
task RetrieveAssembly {
    input {
        String bucket_address
        String filename
        Float haplotype_size
        String reference_fa
    }
    parameter_meta {
        filename: "TAR archive containing the two haplotypes. Includes the <.tar> suffix."
        haplotype_size: "Of one haplotype. In GB."
    }
    command <<<
        set -euxo pipefail
        gsutil cp ~{bucket_address}/~{filename} .
        tar -xf ~{filename}
        rm -f ~{filename}
        FILENAME=$(basename -s .tar ~{filename})
        IFS='_' read -ra TOKENS <<< ${FILENAME}
        echo ${TOKENS[1]} > id1.txt
        echo ${TOKENS[2]} > id2.txt
        echo ${TOKENS[3]} > length.txt
        echo ${TOKENS[4]} > coverage.txt
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
        File h1 = "assembly_" + id1 + "_" + id2 + "_" + length + "_" + coverage + "_h1.fa"
        File h2 = "assembly_" + id1 + "_" + id2 + "_" + length + "_" + coverage + "_h2.fa"
        File config = "config.json"
    }
    runtime {
        docker: "fcunial/simulation"
        disks: "local-disk " + (ceil(haplotype_size*4) + 4) + " HDD"
        preemptible: 3
    }
}


task ArchiveVCF {
    input {
        File vcf
        String bucket_address
        Int id1
        Int id2
        Int length
        Int coverage
    }
    command <<<
        set -euxo pipefail
        FILENAME="pav_~{id1}_~{id2}_~{length}_~{coverage}.vcf"
        cp ~{vcf} ${FILENAME}
        gsutil cp ${FILENAME} ~{bucket_address}/vcfs/
    >>>
    output {
    }
    runtime {
        docker: "fcunial/simulation"
        disks: "local-disk " + (ceil(size(vcf, "GB"))*2) + " HDD"
    }
}
