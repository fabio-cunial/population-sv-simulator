version 1.0

import "JointCalling.wdl"
import "PavWrapper.wdl"


# Given the config files created by the simulator and describing a set of
# simulated haplotypes, this workflow generates the haplotype sequences,
# samples CCS reads from them for several values of coverage and mean length,
# aligns the reads to the reference, and calls SVs.
workflow SimulatePopulation {
    input {
        File reference_fa
        File reference_fai
        File reference_mmi
        File reference_tandem_repeats
        Int n_haplotypes
        File haplotype2variants_file
        File variants_file
        Array[Int] coverages
        Int length_min
        Int length_max
        Int length_stdev
        Array[Int] length_means
        String bucket_address
        Int n_nodes
        Int n_cpus_single_sample
        Int n_cpus_joint
        Int use_pbsv
        Int use_sniffles1
        Int use_sniffles2
        Int use_pav
    }
    parameter_meta {
        reference_fa: "A version of GRCh37 that contains just chr1 and whose header is '>chr1'."
        reference_mmi: "The minimap2 index of <reference_fa>."
        reference_tandem_repeats: "A version that contains just chr1 of <hg19.trf.bed> (downloadable from <https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.trf.bed.gz>)."
        n_haplotypes: "Total number of haplotypes in the simulated population"
        haplotype2variants_file: "Created by the population simulator"
        variants_file: "Created by the population simulator"
        coverages: "Of one haplotype, sorted in increasing order."
        length_means: "Mean read length values, in no particular order."
        n_nodes: "Number of compute nodes to process the population in parallel. Each node processes a chunk of diploid individuals."
        bucket_address: "The full address of the directory in a bucket to be used to store temporary files and output files. The directory is deleted if it already exists."
        use_pbsv: "1=yes, 0=no."
    }
    call DeleteBucketDir {
        input:
            bucket_address = bucket_address
    }
    call GetHaplotypeChunks {
        input:
            n_haplotypes = n_haplotypes,
            n_chunks = n_nodes
    }
    scatter(i in GetHaplotypeChunks.chunks) {
        call ProcessChunkOfHaplotypes {
            input:
                 id_from = i,
                 chunk_size = GetHaplotypeChunks.chunk_size,
                 reference_fa = reference_fa,
                 reference_mmi = reference_mmi,
                 reference_tandem_repeats = reference_tandem_repeats,
                 haplotype2variants_file = haplotype2variants_file,
                 variants_file = variants_file,
                 coverages = coverages,
                 length_min = length_min,
                 length_max = length_max,
                 length_stdev = length_stdev,
                 lengths = length_means,
                 bucket_address = bucket_address,
                 n_cpus = n_cpus_single_sample,
                 use_pbsv = use_pbsv,
                 use_sniffles1 = use_sniffles1,
                 use_sniffles2 = use_sniffles2,
                 use_hifiasm = use_pav
        }
    }
    call JointCalling.JointCalling {
        input:
            bucket_address = bucket_address,
            coverages = coverages,
            lengths = length_means,
            n_individuals = n_haplotypes/2,
            reference_fa = reference_fa,
            n_cpus = n_cpus_joint,
            use_pbsv = use_pbsv,
            use_sniffles2 = use_sniffles2
    }
    if (use_pav == 1) {
        call PavWrapper.PavWrapper {
            input:
                bucket_address = bucket_address,
                n_haplotypes = n_haplotypes,
                coverages = coverages,
                lengths = length_means,
                reference_fa = reference_fa,
                reference_fai = reference_fai
        }
    }
    output {
        Array[File] log_stdout_ProcessChunkOfHaplotypes = ProcessChunkOfHaplotypes.log_stdout
        Array[File] log_stderr_ProcessChunkOfHaplotypes = ProcessChunkOfHaplotypes.log_stderr
        Array[File] log_stdout_JointCalling = JointCalling.logs_stdout
        Array[File] log_stderr_JointCalling = JointCalling.logs_stderr
    }
}


task DeleteBucketDir {
    input { 
        String bucket_address
    }    
    command <<<
        gsutil rm -rf ~{bucket_address}
    >>>
    output {
    }
    runtime {
        docker: "ubuntu:latest"
    }
}


# Returns an array with the ID of the first haplotype to be processed by each
# chunk (zero-based).
task GetHaplotypeChunks {
    input { 
        Int n_haplotypes
        Int n_chunks
    }    
    Int size = ( (n_haplotypes/2 + n_chunks - 1) / n_chunks ) * 2
    command <<<
        seq 0 ~{size} $(( ~{n_haplotypes} - 1 )) > out.txt
    >>>
    output {
        Array[String] chunks = read_lines("out.txt")
        Int chunk_size = size
    }
    runtime {
        docker: "ubuntu:latest"
    }
}


# Given a chunk of simulated haplotype IDs, the task uses every pair of
# consecutive haplotype IDs in the chunk to build a diploid individual and to
# call variants on its read alignments, where reads are simulated at several
# settings of coverage and mean read length.
#
# Remark: the task checkpoints by storing files in <bucket_address>: (1) a
# checkpoint file for each chunk, that contains lines <individual length
# coverage>, is used to skip through the three nested loops of the simulation;
# (2) the files generated by each step are used to skip that step if they
# already exist.
#
# Remark about disk space. For chr1:
# - PBSV's <.svsig.gz> files take from 2MB to 15MB for single-haplotype
#   coverages that range from 4x to 30x.
# - Sniffles2's <.snf> files take 2.5MB regardless of coverage.
# - A single <.vcf> file takes at most 4MB over all coverages.
# Every individual is processed sequentially.
task ProcessChunkOfHaplotypes {
    input {
        Int id_from
        Int chunk_size
        File reference_fa
        File reference_mmi
        File reference_tandem_repeats
        File haplotype2variants_file
        File variants_file
        Array[Int] coverages
        Int length_min
        Int length_max
        Int length_stdev
        Array[Int] lengths
        String bucket_address
        Int n_cpus
        Int use_pbsv
        Int use_sniffles1
        Int use_sniffles2
        Int use_hifiasm
    }
    parameter_meta {
        id_from: "First simulated haplotype to process (zero-based, inclusive)."
        reference_fa: "A version of GRCh37 that contains just chr1 and whose header is '>chr1'."
        reference_mmi: "The minimap2 index of <reference_fa>."
        reference_tandem_repeats: "A version that contains just chr1 of <hg19.trf.bed> (downloadable from <https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.trf.bed.gz>)."
        haplotype2variants_file: "Created by the population simulator"
        variants_file: "Created by the population simulator"
        coverages: "Of one haplotype, sorted in increasing order."
        lengths: "Mean read length values, in no particular order."
        n_cpus: "Used only for setting the runtime."
        use_pbsv: "1=yes, 0=no."
    }
    
    Int min_coverage = coverages[0]
    Int last = length(coverages) - 1
    Int max_coverage = coverages[last]
    Int ram_size_gb = ceil((size(reference_fa, "GB")*max_coverage*2) * 3)  # *3 because of hifiasm
    Int disk_size_image = 20
    Int disk_size_tools = 5
    Int disk_size_gb = disk_size_image + ram_size_gb*2 + ceil( size(reference_fa, "GB") + size(reference_mmi, "GB") + size(reference_tandem_repeats, "GB") + size(haplotype2variants_file, "GB") + size(variants_file, "GB") ) + disk_size_tools
    String work_dir = "/simulation"
    
    command <<<
        set -euxo pipefail
        ID_TO=$(( ~{id_from} + ~{chunk_size} - 1 ))
        
        CHECKPOINT_FILE = "checkpoint_i~{id_from}_i${ID_TO}.txt"
        gsutil -q stat ~{bucket_address}/checkpoints/${CHECKPOINT_FILE}
        if [ $? -eq 0 ]; then
            gsutil cp ~{bucket_address}/checkpoints/${CHECKPOINT_FILE} .
        else
            echo "0 0 0" > ${CHECKPOINT_FILE}
        fi
        COVERAGES=~{sep='-' coverages}
        LENGTHS=~{sep=' ' lengths}
        for ID1 in $(seq ~{id_from} 2 ${ID_TO}); do
            CHECKPOINT_INDIVIDUAL=$(tail -n1 ${CHECKPOINT_FILE} | awk '{ print $1 }')
            if [ ${ID1} -lt ${CHECKPOINT_INDIVIDUAL} ]; then
                continue
            fi
            ID2=$(( ${ID1} + 1 ))
            java -cp ~{work_dir} -Xmx10g PrintChromosomes ${ID1} ${ID2} ~{reference_fa} ~{haplotype2variants_file} ~{variants_file} .
            for LENGTH in ${LENGTHS}; do
                CHECKPOINT_LENGTH=$(tail -n1 ${CHECKPOINT_FILE} | awk '{ print $2 }')
                if [ ${LENGTH} -lt ${CHECKPOINT_LENGTH} ]; then
                    continue
                fi
                LOG_FILE="haplotype2reads_i${ID1}_i${ID2}_l${LENGTH}.log"
                haplotype2reads.sh ${ID1} ${ID2} ~{length_min} ~{length_max} ${LENGTH} ~{length_stdev} ~{max_coverage} ~{bucket_address} &> ${LOG_FILE}
                gsutil cp ${LOG_FILE} ~{bucket_address}/logs/
                LOG_FILE="reads2svs_i${ID1}_i${ID2}_l${LENGTH}.log"
                READS_FILE="reads_i${ID1}_i${ID2}_l${LENGTH}_c~{max_coverage}.fa"
                reads2svs.sh ${READS_FILE} ${ID1} ~{min_coverage} ~{max_coverage} ${COVERAGES} ~{reference_fa} ~{reference_mmi} ~{reference_tandem_repeats} ${CHECKPOINT_FILE} ~{bucket_address} ~{use_pbsv} ~{use_sniffles1} ~{use_sniffles2} ~{use_hifiasm} &> ${LOG_FILE}
                gsutil cp ${LOG_FILE} ~{bucket_address}/logs/
                gsutil rm -f ~{bucket_address}/reads/${READS_FILE}
            done
        done
    >>>
    
    output {
        File log_stdout = stdout()
        File log_stderr = stderr()
    }
    runtime {
        docker: "fcunial/simulation"
        cpu: n_cpus
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 10
    }
}
