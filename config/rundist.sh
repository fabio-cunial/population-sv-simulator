#!/bin/env bash

module load miniconda/4.8.3
module load samtools/1.12
module load minimap2/2.21
module load ucsc/202003
module load lra/1.3.1
module load htslib/1.9
module load seqtk/1.3

mkdir -p log

snakemake -s ${SNAKEFILE} -j ${JOB_COUNT} --nt --ri -k \
    --jobname "{rulename}.{jobid}" \
    --drmaa " -V -cwd -j y -o ./log -l centos=7 -pe serial {cluster.cpu} -l mfree={cluster.mem} -l disk_free={cluster.disk} -l h_rt={cluster.rt} {cluster.params} -w n -S /bin/bash" \
    -w 60 -u ${SITE_CONFIG_DIR}/sge.json "$@"
