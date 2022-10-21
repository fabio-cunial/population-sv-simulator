#!/bin/bash
#
# The workflow in <https://github.com/broadinstitute/pav-wdl> translated into a
# single Bash script with checkpointing at every step.
#
# Max requirements from the original WDL: 8 threads, 32 GB of RAM, 1 GB of disk.
#
ID1=$1
ID2=$2
LENGTH=$3
COVERAGE=$4
REFERENCE_FA=$5
REFERENCE_FAI=$6
HAPLOTYPE1_FA=$7
HAPLOTYPE2_FA=$8
BUCKET_DIR=$9
N_THREADS=${10}

GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
TIME_COMMAND="/usr/bin/time --verbose"
SNAKEMAKE_COMMAND="snakemake -s pav/Snakefile --cores ${N_THREADS}"
SAMPLE_ID="i${ID1}_i${ID2}_l${LENGTH}_c${COVERAGE}"
CHROMOSOME="chr1"

set -euxo pipefail
source activate lr-pav

# Given the string of space-separated local files in input to the Snakemake
# script, the procedure runs the script only if at least one of those files
# does not already exist in the remote bucket (otherwise it downloads all
# files).
function createFiles() {
    local LOCAL_FILES=$1
    MISSING="0"
    for FILE in ${LOCAL_FILES}; do
        TEST=$(gsutil -q stat "${BUCKET_DIR}/pav/${SAMPLE_ID}/${FILE}" || echo 1)
        if [ ${TEST} = 1 ]; then
            MISSING="1"
            break
        fi
    done
    if [ ${MISSING} != 1 ]; then
        for FILE in ${LOCAL_FILES}; do
            mkdir -p $(dirname ${FILE})
            gsutil cp "${BUCKET_DIR}/pav/${SAMPLE_ID}/${FILE}" ${FILE}
        done
    else
        ${TIME_COMMAND} ${SNAKEMAKE_COMMAND} ${LOCAL_FILES}
        for FILE in ${LOCAL_FILES}; do
            gsutil ${GSUTIL_UPLOAD_THRESHOLD} cp ${FILE} "${BUCKET_DIR}/pav/${SAMPLE_ID}/${FILE}"
        done
    fi
}

# Formatting the input
rm -f config.json
echo "{" >> config.json
echo "\"reference\": \"asm/ref.fa\"," >> config.json
echo "\"asm_pattern\": \"asm/{asm_name}/{hap}.fa.gz\"" >> config.json
echo "}" >> config.json
mkdir -p asm/${SAMPLE_ID}
cp ${REFERENCE_FA} asm/ref.fa
cp ${REFERENCE_FAI} asm/ref.fa.fai
gzip -q -c --fast ${HAPLOTYPE1_FA} > asm/${SAMPLE_ID}/h1.fa.gz
gzip -q -c --fast ${HAPLOTYPE2_FA} > asm/${SAMPLE_ID}/h2.fa.gz

# Main pipeline
${SNAKEMAKE_COMMAND} data/ref/ref.fa.gz data/ref/ref.fa.gz.fai
${SNAKEMAKE_COMMAND} temp/${SAMPLE_ID}/align/contigs_h1.fa.gz temp/${SAMPLE_ID}/align/contigs_h1.fa.gz.fai
${SNAKEMAKE_COMMAND} temp/${SAMPLE_ID}/align/contigs_h2.fa.gz temp/${SAMPLE_ID}/align/contigs_h2.fa.gz.fai
createFiles "data/ref/n_gap.bed.gz"
createFiles "temp/${SAMPLE_ID}/align/pre-cut/aligned_tig_h1.sam.gz"
createFiles "temp/${SAMPLE_ID}/align/pre-cut/aligned_tig_h2.sam.gz"
createFiles "results/${SAMPLE_ID}/align/pre-cut/aligned_tig_h1.bed.gz"
createFiles "results/${SAMPLE_ID}/align/pre-cut/aligned_tig_h2.bed.gz"
createFiles "results/${SAMPLE_ID}/align/aligned_tig_h1.bed.gz"
createFiles "results/${SAMPLE_ID}/align/aligned_tig_h2.bed.gz"
for i in $(seq 1 10); do
    createFiles temp/${SAMPLE_ID}/cigar/batched/insdel_h1_${i}.bed.gz temp/${SAMPLE_ID}/cigar/batched/snv.bed_h1_${i}.gz
    createFiles temp/${SAMPLE_ID}/cigar/batched/insdel_h2_${i}.bed.gz temp/${SAMPLE_ID}/cigar/batched/snv.bed_h2_${i}.gz
done
createFiles temp/${SAMPLE_ID}/lg_sv/batch_h1.tsv.gz
createFiles temp/${SAMPLE_ID}/lg_sv/batch_h2.tsv.gz
for i in $(seq 1 10); do
    createFiles temp/${SAMPLE_ID}/lg_sv/batch/sv_ins_h1_${i}.bed.gz temp/${SAMPLE_ID}/lg_sv/batch/sv_del_h1_${i}.bed.gz temp/${SAMPLE_ID}/lg_sv/batch/sv_inv_h1_${i}.bed.gz
    createFiles temp/${SAMPLE_ID}/lg_sv/batch/sv_ins_h2_${i}.bed.gz temp/${SAMPLE_ID}/lg_sv/batch/sv_del_h2_${i}.bed.gz temp/${SAMPLE_ID}/lg_sv/batch/sv_inv_h2_${i}.bed.gz
done
createFiles temp/${SAMPLE_ID}/cigar/pre_inv/svindel_insdel_h1.bed.gz temp/${SAMPLE_ID}/cigar/pre_inv/snv_snv_h1.bed.gz
createFiles temp/${SAMPLE_ID}/cigar/pre_inv/svindel_insdel_h2.bed.gz temp/${SAMPLE_ID}/cigar/pre_inv/snv_snv_h2.bed.gz
createFiles temp/${SAMPLE_ID}/lg_sv/sv_del_h1.bed.gz
createFiles temp/${SAMPLE_ID}/lg_sv/sv_ins_h1.bed.gz
createFiles temp/${SAMPLE_ID}/lg_sv/sv_inv_h1.bed.gz
createFiles temp/${SAMPLE_ID}/lg_sv/sv_del_h2.bed.gz
createFiles temp/${SAMPLE_ID}/lg_sv/sv_ins_h2.bed.gz
createFiles temp/${SAMPLE_ID}/lg_sv/sv_inv_h2.bed.gz
createFiles temp/${SAMPLE_ID}/inv_caller/flag/cluster_indel_h1.bed.gz
createFiles temp/${SAMPLE_ID}/inv_caller/flag/cluster_snv_h1.bed.gz
createFiles temp/${SAMPLE_ID}/inv_caller/flag/cluster_indel_h2.bed.gz
createFiles temp/${SAMPLE_ID}/inv_caller/flag/cluster_snv_h2.bed.gz
createFiles temp/${SAMPLE_ID}/inv_caller/flag/insdel_indel_h1.bed.gz
createFiles temp/${SAMPLE_ID}/inv_caller/flag/insdel_indel_h2.bed.gz
createFiles temp/${SAMPLE_ID}/inv_caller/flag/insdel_sv_h1.bed.gz
createFiles temp/${SAMPLE_ID}/inv_caller/flag/insdel_sv_h2.bed.gz
createFiles results/${SAMPLE_ID}/callable/callable_regions_h1_500.bed.gz
createFiles results/${SAMPLE_ID}/callable/callable_regions_h2_500.bed.gz
createFiles results/${SAMPLE_ID}/inv_caller/flagged_regions_h1.bed.gz
createFiles results/${SAMPLE_ID}/inv_caller/flagged_regions_h2.bed.gz
for i in $(seq 1 60); do
    createFiles temp/${SAMPLE_ID}/inv_caller/batch/h1/inv_call_${i}.bed.gz
    createFiles temp/${SAMPLE_ID}/inv_caller/batch/h2/inv_call_${i}.bed.gz
done
createFiles temp/${SAMPLE_ID}/inv_caller/sv_inv_h1.bed.gz
createFiles temp/${SAMPLE_ID}/inv_caller/sv_inv_h2.bed.gz
createFiles temp/${SAMPLE_ID}/bed/integrated/h1/svindel_ins.bed.gz temp/${SAMPLE_ID}/bed/integrated/h1/svindel_del.bed.gz temp/${SAMPLE_ID}/bed/integrated/h1/snv_snv.bed.gz temp/${SAMPLE_ID}/bed/integrated/h1/sv_inv.bed.gz
createFiles temp/${SAMPLE_ID}/bed/integrated/h2/svindel_ins.bed.gz temp/${SAMPLE_ID}/bed/integrated/h2/svindel_del.bed.gz temp/${SAMPLE_ID}/bed/integrated/h2/snv_snv.bed.gz temp/${SAMPLE_ID}/bed/integrated/h2/sv_inv.bed.gz
createFiles temp/${SAMPLE_ID}/bed/bychrom/svindel_ins/${CHROMOSOME}.bed.gz
createFiles temp/${SAMPLE_ID}/bed/bychrom/svindel_del/${CHROMOSOME}.bed.gz
createFiles temp/${SAMPLE_ID}/bed/bychrom/sv_inv/${CHROMOSOME}.bed.gz
createFiles temp/${SAMPLE_ID}/bed/bychrom/snv_snv/${CHROMOSOME}.bed.gz
createFiles temp/${SAMPLE_ID}/bed/merged/snv_snv.bed.gz
createFiles temp/${SAMPLE_ID}/bed/merged/sv_inv.bed.gz
createFiles temp/${SAMPLE_ID}/bed/merged/svindel_ins.bed.gz
createFiles temp/${SAMPLE_ID}/bed/merged/svindel_del.bed.gz
createFiles results/${SAMPLE_ID}/bed/snv_snv.bed.gz results/${SAMPLE_ID}/bed/indel_ins.bed.gz results/${SAMPLE_ID}/bed/indel_del.bed.gz results/${SAMPLE_ID}/bed/sv_ins.bed.gz results/${SAMPLE_ID}/bed/sv_del.bed.gz results/${SAMPLE_ID}/bed/sv_inv.bed.gz results/${SAMPLE_ID}/bed/fa/indel_ins.fa.gz results/${SAMPLE_ID}/bed/fa/indel_del.fa.gz results/${SAMPLE_ID}/bed/fa/sv_ins.fa.gz results/${SAMPLE_ID}/bed/fa/sv_del.fa.gz results/${SAMPLE_ID}/bed/fa/sv_inv.fa.gz
createFiles data/ref/contig_info.tsv.gz
createFiles pav_${SAMPLE_ID}.vcf.gz

# Output
gunzip pav_${SAMPLE_ID}.vcf.gz
gsutil cp pav_${SAMPLE_ID}.vcf ${BUCKET_DIR}/vcfs/
rm -rf asm/ data/ temp/ results/
gsutil -m rm -f "${BUCKET_DIR}/pav/${SAMPLE_ID}/"
conda deactivate
