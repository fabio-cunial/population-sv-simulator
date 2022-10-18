#!/bin/bash
#
# Creates a diploid individual by taking two haplotypes with consecutive IDs.
# Runs the read simulator at max coverage on each. Randomly permutes its
# output, to make sure the two haplotypes are mixed when splitting the file
# later.
#
# Remark: for chr1 at 30x coverage, the procedure adds approx. 26 GB to the
# current disk usage, the output file is ~13 GB, and the peak RAM is ~13 GB
# (from <shuf>, which loads the whole file in memory; <pbsim> takes just 240 MB
# of RAM).
#
ID1=$1
ID2=$2
LENGTH_MIN=$3
LENGTH_MAX=$4
LENGTH_MEAN=$5
LENGTH_STDEV=$6
MAX_COVERAGE=$7
BUCKET_ADDRESS=$8  # Root dir of the simulation in the bucket

TIME_COMMAND="/usr/bin/time --verbose"
PBSIM_MODEL="/simulation/pbsim_data/model_qc_ccs"
PBSIM_COMMAND="pbsim --data-type CLR --length-min ${LENGTH_MIN} --length-max ${LENGTH_MAX} --accuracy-min 0.99 --accuracy-max 1.0 --difference-ratio 6:21:73 --model_qc ${PBSIM_MODEL} --length-mean ${LENGTH_MEAN} --length-sd ${LENGTH_STDEV} --accuracy-mean 0.995 --accuracy-sd 0.002"

set -euxo pipefail
echo "Running <haplotype2reads.sh> on the following node:"
lscpu
cat /proc/meminfo

OUTPUT_FILE="reads_i${ID1}_i${ID2}_l${LENGTH_MEAN}_c${MAX_COVERAGE}.fa"
if [ $(gsutil -q stat ${BUCKET_ADDRESS}/reads/${OUTPUT_FILE}) -eq 0 ]; then
    ${TIME_COMMAND} gsutil cp ${BUCKET_ADDRESS}/reads/${OUTPUT_FILE} .
else
    ${TIME_COMMAND} ${PBSIM_COMMAND} --depth ${MAX_COVERAGE} --prefix reads_${ID1} haplotype_${ID1}.fa
    mv reads_${ID1}_0001.fastq reads_${ID1}.fa
    rm -f reads_${ID1}_0001.* 
    rm -f haplotype_${ID1}.fa
    ${TIME_COMMAND} ${PBSIM_COMMAND} --depth ${MAX_COVERAGE} --prefix reads_${ID2} haplotype_${ID2}.fa
    mv reads_${ID2}_0001.fastq reads_${ID2}.fa
    rm -f reads_${ID2}_0001.* 
    rm -f haplotype_${ID2}.fa
    cat reads_${ID2}.fa >> reads_${ID1}.fa
    rm -f reads_${ID2}.fa
    ${TIME_COMMAND} shuf reads_${ID1}.fa > reads_${ID1}_${ID2}.1
    rm -f reads_${ID1}.fa
    tr 'Z' '\n' < reads_${ID1}_${ID2}.1 > ${OUTPUT_FILE}
    rm -f reads_${ID1}_${ID2}.1
    ${TIME_COMMAND} gsutil cp ${OUTPUT_FILE} ${BUCKET_ADDRESS}/reads/
fi
