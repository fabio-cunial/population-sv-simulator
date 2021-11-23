#!/bin/env bash

module load miniconda/4.8.3
module load samtools/1.10
module load minimap2/2.21
module load ucsc/202003
module load lra/1.3.1
module load htslib/1.9

snakemake -s ${SNAKEFILE} -j ${JOB_COUNT} --nt --ri -k \
    --jobname "{rulename}.{jobid}" \
    -w 60 "$@"
