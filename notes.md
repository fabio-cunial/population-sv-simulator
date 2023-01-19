KIRAN:
joint calling of SVs	
simulation				1000 samples simulation is enough for now, we should do more only later.

---------------------------------------------

Space for the chromosomes:
273 MB per sample * 10000 individuals * 2 chromosomes per individual = 5,460 GB 




badread simulate --reference ${CHROMOSOME} --quantity ${COVERAGE}x --length 15000,2000 --identity 99.9,100,0.05 --error_model pacbio2016 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 --start_adapter_seq "" --end_adapter_seq "" | gzip > reads.fastq.gz


https://github.com/rrwick/Badread



# PBSIM: downloaded from https://github.com/pfaucon/PBSIM-PacBio-Simulator

# Command for the Badread simulator:
# badread simulate --reference ${POPULATION_DIR}/${CHROMOSOME_PREFIX}_${ID1}.fa --quantity ${MAX_COVERAGE}x --length ${LENGTH_MEAN},${LENGTH_STDEV} --identity 99.9,100,0.05 --error_model pacbio2016 --junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 --start_adapter_seq "" --end_adapter_seq "" > ${PREFIX}_${MAX_COVERAGE}_${ID1}.fastq



Sniffles2 Error: A tandem repeat annotations file was provided, but no matching annotations were found for any contig in the sample input file. Please check if the contig naming scheme in the tandem repeat annotations matches with the one in the input sample file. (Fatal error, exiting.)



pbsv discover ERROR: Hard clips are not supported.


sniffles1: Please select only vcf OR bedpe output format!





# function reads2svs() {
# 	local READS_FILE_PREFIX=$1  # Without .fa
# 	local SAMPLE_ID=$2  # SM field in the .sam file (needed for joint calling)
#
# 	local READ_GROUP="@RG\tID:movie\tSM:${SAMPLE_ID}"
# 	${MINIMAP_COMMAND} -R ${READ_GROUP} ${REFERENCE_FA} ${READS_FILE_PREFIX}.fa > ${READS_FILE_PREFIX}.sam
# 	rm -f ${READS_FILE_PREFIX}.fa
# 	samtools view -b ${READS_FILE_PREFIX}.sam > ${READS_FILE_PREFIX}.bam
# 	rm -f ${READS_FILE_PREFIX}.sam
# 	samtools sort ${READS_FILE_PREFIX}.bam > ${READS_FILE_PREFIX}.sorted.bam
# 	rm -f ${READS_FILE_PREFIX}.bam
# 	mv ${READS_FILE_PREFIX}.sorted.bam ${READS_FILE_PREFIX}.bam
# 	samtools index ${READS_FILE_PREFIX}.bam
#
# 	local SIGNATURE_FILE="${SIGNATURE_DIR}/$(basename ${READS_FILE_PREFIX}).svsig.gz"
# 	${PBSV_COMMAND} discover --tandem-repeats ${TANDEM_REPEATS_FILE} ${READS_FILE_PREFIX}.bam ${SIGNATURE_FILE}
# 	local VCF_FILE="${VCF_DIR}/$(basename ${READS_FILE_PREFIX})_pbsv.vcf"
# 	${PBSV_COMMAND} call --ccs ${REFERENCE_FA} ${SIGNATURE_FILE} ${VCF_FILE}
#
# 	VCF_FILE="${VCF_DIR}/$(basename ${READS_FILE_PREFIX})_sniffles1.vcf"
# 	${SNIFFLES1_COMMAND} -t 1 -m ${READS_FILE_PREFIX}.bam -v ${VCF_FILE}
# 	SIGNATURE_FILE="${SIGNATURE_DIR}/$(basename ${READS_FILE_PREFIX}).snf"
# 	VCF_FILE="${VCF_DIR}/$(basename ${READS_FILE_PREFIX})_sniffles2.vcf"
# 	${SNIFFLES2_COMMAND} --threads 1 --tandem-repeats ${TANDEM_REPEATS_FILE} --reference ${REFERENCE_FA} --input ${READS_FILE_PREFIX}.bam --vcf ${VCF_FILE} --snf ${SIGNATURE_FILE}
#
# 	rm -f ${READS_FILE_PREFIX}.bam
# }


# function processChunkOfChromosomes() {
# 	local ID_FROM=$1  # Inclusive
# 	local ID_TO=$2  # Inclusive
# 	local PREFIX="${READS_DIR}/${CHROMOSOME_PREFIX}"
#
# 	for ID1 in $(seq ${ID_FROM} 2 ${ID_TO}); do
# 		ID2=$(( ${ID1} + 1 ))
# 		chromosome2reads ${ID1} ${ID2} &> ${LOGS_DIR}/chromosome2reads_${ID1}_${ID2}.log
# 		N_ROWS=$(wc -l < ${PREFIX}_${MAX_COVERAGE}_${ID1}_${ID2}.fa)
# 		for COVERAGE in ${COVERAGES}; do
# 			N_ROWS_PRIME=$(( ( (${N_ROWS}*${COVERAGE}) / (${MAX_COVERAGE}*2) ) * 2 ))
# 			head -n ${N_ROWS_PRIME} ${PREFIX}_${MAX_COVERAGE}_${ID1}_${ID2}.fa > ${PREFIX}_${COVERAGE}_${ID1}_${ID2}.fa
# 			reads2svs ${PREFIX}_${COVERAGE}_${ID1}_${ID2} ${ID1} &> ${LOGS_DIR}/reads2svs_${COVERAGE}_${ID1}_${ID2}.log
# 		done
# 		reads2svs ${PREFIX}_${MAX_COVERAGE}_${ID1}_${ID2} ${ID1} &> ${LOGS_DIR}/reads2svs_${MAX_COVERAGE}_${ID1}_${ID2}.log
# 	done
# }



buildChromosomes
\_ referenceLength, reference
\_ chromosome2variants
\_ variants






		
	/**
	 * Variants
	 */
	private static Variant[] variants;
	private static int lastVariant, nVariants_frequent, nVariants_rare;
	
	/**
	 * For every chromosome (rows), the sorted list of variants it contains.
	 */
	private static int N_CHROMOSOMES;
	
	
	
[DONE] try running without PAV first, in the cloud.

[DONE] reengineer PAV wdl into a bash script

[DONE] above: checkpoint every file

[DONE] checkpointing files: create a function that takes care of this also in the other scripts. NO, IN THE OTHER SCRIPTS THERE IS A DIFFERENT COMMAND EACH TIME.

[NO, TOO COMPLICATED] pav: organize by chunks as well, by creating a chunk-workflow that loops and invokes the pav wdl many times. however this is probably not guaranteed to be on the same machine, it's just a way to limit excessive parallelism.

[DONE] force sequentiality in pav

[DONE] force sequentiality in SimulateHaplotypes

[DONE] reads2svs.sh: the bottleneck of the whole process is minimap2. so we should store in the bucket every single aligned chunk for robustness, otherwise we lose hours of computation if a failure occurs in that loop.

[DONE] above: remove the remote reads and alignments file only at the very end of reads2svs, not in SImulateHaplotypes. since that might be performed even when reads2svs did not complete.

[DONE] above: all the logs, even in the main script, should be archived transparently to the bucket. WRONG: leave all the logs to std so that they can be seen live in terra.

[DONE] save assembly statistics to the bucket. the actual assemblies are not necessary, but give a flag that can save them, if we need to inspect them later.  => QUAST

[DONE] use pbsim2, just to say we are using the latest, but resort to the same CLR trick we are currently using for pbsim1. they also have better length distribution and the new quality distribution affects the errors inserted in the reads. although the error is probably not important in our CCS case.

[DONE] above: give also the option of no error, i.e. create a version of pbsim that does not introduce errors but keeps the length distribution. IN THEORY WE COULD JUST SET ERROR RATE TO ZERO ON THE COMMAND LINE... THIS IS JUST TO BE COMPLETELY SURE THERE IS NO ERROR.

[DONE] above: there is also PBSIM3, which simulates CLR and feeds them to the CCS algorithm. but it does not seem to be an accepted paper yet. add it as soon as it becomes accepted.

[DONE] in the final simulation, replace chr1 with the T2T chr1. i.e. allow the haplotypes to be generated on a different, more realistic chromosome.

[DONE] make sure we run the pav script (very expensive) only if the pav vcf is missing from the remote bucket.

[DONE] above: make sure we download the assemblies from the remote bucket only if we need to run pav.

[DONE] pav: get rid of all the batches.

[DONE] pav: set all thread numbers.

[DONE] create MMI file for T2T ref

[DONE] make sure we are saving quast output

[DONE] make sure that PAV outputs insertion sequences in the VCF (otherwise we have to save some auxiliary file as well).

[DONE] PAV VCFs are huge: remove the short variants.

[DONE] implement checkpointing of PAV: function that loops over fixed set of output files before calling PAV. but we also need a function that every few minutes checks is those files exist, if they have not been touched for some time, and if they were absent from the bucket at the beginning (keep array of booleans in order not to query the bucket every time).

[DONE] fix coverage file construction bug in main loop of read2svs.sh when resuming from preemption.

[DONE] no space left on device after PAV and concatenating BAM chunks in the following iteration?!?! do not run under /simulation (small HD of the image), run under /cromwell_root (data directory in the allocated HD).
Filesystem Size Used Avail Use% Mounted on
overlay 17G 16G 1.4G 92% /
tmpfs 64M 0 64M 0% /dev
shm 32G 0 32G 0% /dev/shm
/dev/disk/by-id/google-local-disk 89G 938M 88G 2% /cromwell_root
/dev/sda1 17G 16G 1.4G 92% /google
tmpfs 16G 0 16G 0% /proc/acpi
tmpfs 16G 0 16G 0% /proc/scsi
tmpfs 16G 0 16G 0% /sys/firmware

[DONE] above: the requested HD should not count the size of the image.

- above: say in the documentation that PAV can be checkpointed, but there is a single long step that does not produce any file and that lasts 1h. so if preemption happens there, we lose 1h for every preemption.

[DONE] update documentation of WDL by removing explicit references to GRCh37

[DONE] reduce preemption attempts to 3, since they significantly slow down PAV because of the above.

- above: preemptibility should be chosen by the caller. set high preemtability if the longest atomic task is not very long (that is the speed price we pay after preemption).

- PRs with all the fixes to PAV.

[DONE] replace our TRF for T2T with the one provided by pbsv at https://github.com/PacificBiosciences/pbsv/tree/master/annotations

[DONE] do detailed time and memory analysis now that we have the actual times of everything. How much disk/RAM does hifiasm really take? workflow failed but we can still use it for timing, I think.
workflow id: dbb2fcf2-118d-42b8-884e-0d8b932c3632

[DONE] in the final simulation, take the read length mean/stdev from the technical pilot? since those are in the terra tables. ACTUALLY TAKEN FROM GELB AND OTHER SIMILAR DATASETS WITH MORE SAMPLES THAN THE TECH PILOT.

[DONE] Start with Truvari first, just to make sure people trust the results.

[DONE] above: use negative SVLEN for deletions also in the ground truth VCFs.

[DONE] end is smaller than pos in the ground truth VCFs?!

[DONE] make sure we output in the ground truth simple INS and DEL calls conforming to the VCF spec.

[DONE] make sure we also compare VCFs without filtering by SVTYPE or context.

[DONE] parallelize over the 1000 files of each individual processed in each sequential iteration, using the many cores of the node.

[WRONG] indexing the ground truth VCFs must be done at the very beginning, just one time, in a separate WDL task. NO, BECAUSE WE FILTER THE VCF FILES EVERY TIME BY THE CRITERIA IN THE MEASURED FILE.

[DONE] make indexing of measured files fault tolerant. IT'S ALREADY FAULT TOLERANT AT THE LEVEL OF COMPUTING THE MATRICES FOR EACH CONFIGURATION. MOST STEPS IN THE BASH SCRIPT ARE PROBABLY FAST. THE ONLY SLOW ONE IS MERGING ALL VCFS IN THE POPULATION, WHICH IS NOW FAUT TOLERANT.

[DONE] above: make sure that all key steps are timed.

[DONE] try to merge all single-sample calls from all individuals into a single union, and compare it to the VCF that contains the distinct variants in the population.

[DONE] above: make sure truvari collapse and bcftools merge flags are all correct.

[DONE] above: handle also joint calls by pbsv and sniffles2

[DONE] above: can the workflow run on incomplete output? can it work without joint calling files? YES TO BOTH.

[DONE] add truvari et al. to the docker. add new workflows to the dockstore descriptor

[DONE] document the checkpointing strategy in PerformanceMatrices WDL.

[DONE] create program that annotates with repeat context and surface every experimental VCF, like the ground truth VCF.

[DONE] above: this should be a separate task that is run in parallel just once.

[DONE] above: make sure that we are 0-based but the VCF is 1-based. NOT SO IMPORTANT, ANNOTATIONS HAVE TOLERANCES.

[DONE] above: handle case where length and SVEND do not match.

[DONE] above: handle insertions!

[DONE] above: loading the data structures is the time bottleneck. make the program loop over many VCFs.

[DONE] above: test it on ALL files we have already generated.

[DONE] above: randomize files in chunks, otherwise if there is an order to gsutil ls, the work might be imbalanced.

[DONE] in PerformanceMatrices: filter also by SV length and SV frequency in the population (for the latter we can only measure recall). These should be at the same level as SV type in the current process.

[ALREADY IMPLEMENTED] above: make use of joint-calling matrices optional, since only some of them have been computed.

[DONE] above: randomize order of those workpackages as well...

[NOT NEEDED] above: actually, randomize every time we create chunks throughout the code. OTHER CHUNKS IN THE CODE ARE ALREADY BALANCED.

- reduce RAM of PBSV joint if it is smaller than allocated even at max coverage. if RAM depends on coverage, one should allocate it better rather than one size for each coverage...


[DONE] fix WDL bug:
Required file output '/cromwell_root/force_sequentiality.txt' does not exist.

[DONE] handle gsutil cp errors throughout the simulation code:
+ gsutil cp pbsv_i330_i331_l22500_c20.svsig.gz gs://cunial-bucket-1/simulation-root/run/signatures/
Copying file://pbsv_i330_i331_l22500_c20.svsig.gz [Content-Type=application/octet-stream]...
/ [0 files][ 0.0 B/ 9.9 MiB] ResumableUploadException: 503 Server Error
2022/11/04 20:42:33 Starting delocalization.

[DONE] above: handle upload/download errors also in PerformanceMatrices

[DONE] call AnnotateVCFs workflow from PerformanceMatrices, optionally.

[DONE] make sure we are not splitting into more chunks than there are lines in the file

[DONE] make JointCalling its own script that can be executed independently in case the main one fails.

[DONE] above: make sniffles joint preemptible, but pbsv joint not since it's very slow. in general, separate sniffles2 and pbsv into distinct tasks with resources optimized for each.

[DONE] should the FA of pbsv joint have a FAI? NO, NOT NECESSARY, EVERYTHING RUNS WITHOUT IT.

[DONE] say that in PeformanceMatrices we do not try to handle BNDs, even though maybe that would improve performance.

[DONE] do not split by nChunks, since it cuts in the middle of lines! always split by number of lines!

[DONE] PerformanceMatrices: add lines in the header of experimental VCFs. 
[DONE] above: add a note to the java program that annotates those VCFs.

[DONE] PerformanceMatrices: groundTruth_joint.vcf should be compressed just once, rather than every time the filter string is empty.

[DONE] above: reference.fa,fai should be copied only once, not every time the shell script sis executed. this is a very expensive network transfer!

[DONE] PerformanceMatrices: handle negative SVLEN in filters!

[DONE] PerformanceMatrices: sort VCFs before feeding them to truvari or merging them.

[DONE] above: add note about VCF sorting in comments and possibly in java code.

[DONE] we could avoid sorting ground truth VCFs.

[NOT NEEDED, IT'S SAFER THIS WAY] PerformanceMatrices : sed filters should operate only on headers.

[DONE] PerformanceMatrices: disk size should be increased, since now we create several temp files in parallel. and we also have stored the merged file and the joint file.

[DONE] PerformanceMatrices: in the bash script, copy FAI as well, otherwise it is created in parallel by every thread.

[DONE] PerformanceMatrixes: make sure parentheses in filter strings are correct!

[DONE] PerformanceMatrices: make svLengths, svFrequencies, repeatFractions cumulative filters.
FILTER_STRING="((SVLEN>${previousSvLength} && SVLEN<=${svLength}) || (SVLEN>=-${svLength} && SVLEN<-${previousSvLength}))"

[DONE] PerformanceMatrices: allow some filter arrays to be empty so the user does not have to use all filtering conditions every time.

[DONE] PerformanceMatrices: do also repeat-only analysis at the highest level of the decision tree.

[DONE] PerformanceMatrices: fix chunk number bug. this applies to other parts of the code as well.

[DONE] PAFtools could be used as an alternative to PAV. might be faster since it's part of minimap2? needs conversion from its format to VCF?

[DONE] above: write PAFtools parser that keeps just SVs and produces a decent VCF? even if truvari works with it out of the box, we still have to add SVLEN etc for our filters.


- above: PAFtools should be run on each haplotype as well, since the output of using both haplotypes concatenated is much smaller than using each haplotype in isolation.



- make a version of PerformanceMatrices that treats DUPs as INSs (e.g. paftools needs that). unclear how to handle INVs...






-----> PRIORITIZE TRIOS AOU AND READ LENGTH AND COVERAGE FOR ANSWERING GP QUESTIONS IMMEDIATELY. 




- create tools to visualize assembly stats.


[DONE] organize octave code into functions that take a path as input.

[DONE] above: find a way to make fonts larger and set the axes to the same scale.

- PerformanceMatrices: ask for deleting the bucket dir, optional.



- how preemptible should the analysis be? DEPENDS ON THE LENGTH OF THE SLOWEST STEP...

[DONE] above: do detailed performance analysis of PerformanceMatrices wdl.

[DONE] test analysis workflow ASAP to catch bugs in the simulation output and to produce figures. 

[DONE] run an experiment with few individuals to test joint calling. If joint calling does not work, only then abort the current large-population run.

[DONE] remember that SVLEN is negative in deletions.

[DONE] fix SV length bug: we should focus on long SVs rather than on short.


- samtools calmd is faster if run after having sorted the bam by position. apply this fix to the simulator. remark: calmd is needed by sniffles1.


https://github.com/PacificBiosciences/sv-benchmark



to run a docker image when the GUI does not work:
docker run -dit fcunial/assemblybased




- according to Kiran there is a joint version of sniffles2, but the authors recommend no to use it.


------------------------------- KIRAN NOV 18 -----------------------------------

- read lengths [5k, 6, 7, 8, 9]. like in the critical length paper.

OBJECTIVES OF JOINT CALLING:
- joint calling: recall.
- joint calling: square off genotyping matrix. no calls per sample.

[DONE] allele freq: actually use raw counts of individuals. do more samples in joint calling become useless if we already have enough, say 3-4? because we reach the high-coverage regime of that SV.

- mine presentations on haplotype caller joint calling, completely subsume them.
- motivate my joint calling method with strong arguments
- specific examples that are missed by joint calling

- study defaults of pbsv and sniffles2

- use alternative to gnomAD-SV because short reads
- icelander genomes

[DONE] no need for other callers for now

- add results for PAV.

[DONE] share github with kiran
- merge it with long reads pipelines

- complex events: see which ones are handled by sniffles2 and pbsv

[DONE] FIX PLOTS BUG: STARS ARE NOT ALIGNED!!!!!!!!

- compare to critical length in long-read resequencing paper. completely subsume it.

- to answer Steve's questions: plot more info about the distinct simulated variants, e.g. how much they overlap, to inform truvari bench parameters.

[DONE] above: plot distribution of SVs along the chr1 positions. to see how much they occur in centromeres/telomeres.

- above: learn also the position bias from gnomAD-SV, since it was there in the original paper. replace it in the step where we sample a position at random.

- above: make position on chr1 also a filtering option, to study what happens in centromere/telomeres.

(FABIO'S IDEA: could we contact pacbio for low performance joint?)

- real All of Us trios at high coverage. we can downsample them and do the "simulation" on them. (50 trios for now) do this before Ryan comes, and before the revio device arrives. so in 2-3 months.

- evaluate SV calls in the sense of QC of SV calls: there is a lot of software already made by the SV team to do that. extensive evaluation.

- outline plan for Ryan for joint calling

--------------------------------------------------------------------------------


