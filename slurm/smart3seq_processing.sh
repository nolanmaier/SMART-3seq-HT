#! /bin/bash

#################################################
# SMART3SEQ PROCESSING SCRIPT
# this script processes input fastq.gz files 
# will demultiplex, trim, align, deduplicate, and peform some QC
#################################################


#################################################
# SETUP
#################################################

# define a function to display help information for the script
Help()
{
  echo 'Syntax: smart3seq_processing.sh -i [input directory] -o [output directory] -t [TSO barcode fasta file] -s [STAR genome index] -g [gtf reference file] -b [bed reference file]'
}

# get the input options from the command line arguments
while getopts ":i:o:t:s:g:b:h" flag; do
  case "${flag}" in
  	i) INDIR=${OPTARG} ;;
    o) OUTDIR=${OPTARG} ;;
		t) TSO_FILE=${OPTARG} ;;
		s) GENOMEDIR=${OPTARG} ;;
    g) GTFFILE=${OPTARG} ;;
		b) BEDFILE=${OPTARG} ;;
    h) Help; exit ;;
    :) echo "$0: Must supply an argument to -${OPTARG}"; exit 1;;
    \?) echo "Invalid option: -${OPTARG}"; Help; exit 1;;
  esac
done
shift "$(($OPTIND -1))"

# "SECONDS" is a special bash variable that tracks the time
SECONDS=0
echo "smart3seq processing script started at $(date)"
# print out the command submitted to initialize this script
echo 'Command Submitted:'
echo $(scontrol show job "${SLURM_JOBID}" | awk -F= '/Command=/{print $2;exit;}')
echo

# load biogrids module to access analysis software
# tools used: fqtools, cutadapt, fastqc, STAR, samtools, umi_tools, rseqc, qualimap
module load biogrids/latest
# harcode the location of the RSeQC scripts from the biogrids install
RSEQC_SCRIPTS='/programs/x86_64-linux/rseqc/4.0.0/bin.capsules'
echo

# print versions of software used for documentation
echo 'Programs used:'
fqtools -v
echo "Cutadapt version: $(cutadapt --version)"
fastqc -v
echo "STAR version: $(STAR --version)"
samtools --version | head -n 1
umi_tools -v
echo "RSEQC version: $(${RSEQC_SCRIPTS}/bam_stat.py --version) $(${RSEQC_SCRIPTS}/infer_experiment.py --version) $(${RSEQC_SCRIPTS}/read_distribution.py --version)"
qualimap --help 2> '/dev/null' | sed -n '/^QualiMap/p'
echo

# each job of the array will run one fastq file found in the input directory
RAW_FASTQ=$(find "${INDIR}" -maxdepth 1 -type f -size +0 -name '*.fastq.gz' | sed -n "${SLURM_ARRAY_TASK_ID}"p)
# define the filename prefix
RAW_PREFIX="$(basename ${RAW_FASTQ} .fastq.gz)"
echo 'This script will process:'
echo ${RAW_FASTQ}
echo

#################################################
# CUTADAPT demultiplexing using barcodes file
#################################################
SECTION_START="${SECONDS}"

# before starting check the input fastq file using fqtools
echo "Validating: $(basename ${RAW_FASTQ})"
fqtools validate "${RAW_FASTQ}"
echo "File contains $(fqtools count ${RAW_FASTQ}) reads"
echo

# make output folder
mkdir -p "${OUTDIR}/demultiplexed_fastq"
# find the outfile name prefix
PREFIX="${OUTDIR}/demultiplexed_fastq/$(basename ${RAW_FASTQ} .fastq.gz)"

# run cutadapt demultiplexing
echo "Demultiplexing: $(basename ${RAW_FASTQ})"
cutadapt --cores="${SLURM_CPUS_PER_TASK}" -g ^file:"${TSO_FILE}" -e 1 --no-indels \
	--untrimmed-output "${OUTDIR}/demultiplexed_fastq/${RAW_PREFIX}_TSONaN_OFF0.fastq.gz" \
	-o "${OUTDIR}/demultiplexed_fastq/${RAW_PREFIX}_{name}.fastq.gz" \
	"${RAW_FASTQ}" \
	> "${OUTDIR}/demultiplexed_fastq/${RAW_PREFIX}.cutadapt_demulti.log"
echo "Demultiplexed into" $(ls -1q ${OUTDIR}/demultiplexed_fastq/${RAW_PREFIX}_*.fastq.gz | wc -l) "files"
echo 'DONE'
echo

echo "Section finished in "$((${SECONDS}-${SECTION_START}))" seconds"
echo


#################################################
# FASTQC on demultiplexed files
#################################################
SECTION_START="${SECONDS}"

# cutadapt outputs an empty file for barcodes with zero reads, delete these before downstream analysis
# otherwise validate output file with fqtools
# uses xargs -P to run multiple parallel processes and reduce runtime
find "${OUTDIR}/demultiplexed_fastq" -maxdepth 1 -type f -name "${RAW_PREFIX}_*.fastq.gz" -print0 | \
	xargs -0 -I '{}' -P "${SLURM_CPUS_PER_TASK}" bash -c \
	'if test "$(fqtools count $1)" -eq 0; \
	then echo ""$(basename $1)" --> EMPTY AND REMOVED"; rm "$1"; \
	else echo ""$(basename $1)" --> VALIDATION: "$(fqtools validate $1)""; \
	fi' \
	-- '{}'
echo 'DONE'
echo

# make output folder
mkdir -p "${OUTDIR}/fastqc"

# run fastqc
echo 'Running FastQC on demultiplexed fastq files...'
fastqc "${OUTDIR}/demultiplexed_fastq/${RAW_PREFIX}_"*".fastq.gz" \
	--outdir="${OUTDIR}/fastqc" --noextract --nogroup --quiet --threads "${SLURM_CPUS_PER_TASK}"
# rename output files
echo "Renaming FastQC logs..."
for FASTQC_HTML in "${OUTDIR}/fastqc/${RAW_PREFIX}_"*"_fastqc.html"; do
	mv "${FASTQC_HTML}" $(echo "${FASTQC_HTML}" | sed 's/_fastqc.html/.raw_fastqc.html/'); done
for FASTQC_ZIP in "${OUTDIR}/fastqc/${RAW_PREFIX}_"*"_fastqc.zip"; do
	mv "${FASTQC_ZIP}" $(echo "${FASTQC_ZIP}" | sed 's/_fastqc.zip/.raw_fastqc.zip/'); done
echo 'DONE'
echo

echo "Section finished in "$((${SECONDS}-${SECTION_START}))" seconds"
echo

#################################################
# CUTADAPT trimming and UMI extraction
#################################################
SECTION_START="${SECONDS}"

# make output folder
mkdir -p "${OUTDIR}/trimmed_fastq"

# iterate through the raw demultiplexed fastq files
for DEMULTI_FASTQ in "${OUTDIR}/demultiplexed_fastq/${RAW_PREFIX}_"*".fastq.gz"; do
	echo "Using Cutadapt to trim: $(basename ${DEMULTI_FASTQ})"
	# find offset using filename
	OFFSET="${DEMULTI_FASTQ##*OFF}"
	OFFSET="${OFFSET%%.fastq.gz}"
	# cutadapt runs three times on the same file:
	# 1. trim offset (no logs or intermediate files are kept)
	# 2. trim UMI and append UMI sequence to read name (no logs or intermediate files are kept)
	# 3. a) trim 13 bases from the 5' end (corresponding to the TSO constant sequence)
	# 3. b) trim low quality bases from the 3' end (assuming 2-color nextseq chemistry)
	# 3. c) trim Illumina adapters and polyA sequences from the 3' end
	# 3. d) filter out any reads shorter than a minimum length
	# 3. e) output a new fastq file and a log file
	cutadapt --cores="${SLURM_CPUS_PER_TASK}" --quiet --cut "${OFFSET}" "${DEMULTI_FASTQ}" | \
	cutadapt --cores="${SLURM_CPUS_PER_TASK}" --quiet --cut 5 --rename='{id}_{cut_prefix} {comment}' - | \
	cutadapt --cores="${SLURM_CPUS_PER_TASK}" --cut 13 --nextseq-trim=15 --minimum-length 20 \
	-a 'polyA=A{20}' -a 'illumina=AGATCGGAAGAGC' --times 2 \
	-o "${OUTDIR}/trimmed_fastq/$(basename ${DEMULTI_FASTQ} .fastq.gz).trimmed.fastq.gz" - \
	> "${OUTDIR}/trimmed_fastq/$(basename ${DEMULTI_FASTQ} .fastq.gz).cutadapt.log"
done
echo 'DONE'
echo

echo "Section finished in "$((${SECONDS}-${SECTION_START}))" seconds"
echo

#################################################
# FASTQC on trimmed files
#################################################
SECTION_START="${SECONDS}"

# check each output file and delete it if empty
# otherwise validate output file with fqtools
# uses xargs -P to run multiple parallel processes and reduce runtime
find "${OUTDIR}/trimmed_fastq" -maxdepth 1 -type f -name "${RAW_PREFIX}_*.trimmed.fastq.gz" -print0 | \
	xargs -0 -I '{}' -P "${SLURM_CPUS_PER_TASK}" bash -c \
	'if test "$(fqtools count $1)" -eq 0; \
	then echo ""$(basename $1)" --> EMPTY AND REMOVED"; rm "$1"; \
	else echo ""$(basename $1)" --> VALIDATION: "$(fqtools validate $1)""; \
	fi' \
	-- '{}'
echo 'DONE'
echo

# make output folder
mkdir -p "${OUTDIR}/fastqc"

# run fastqc
echo 'Running FastQC on processed fastq files...'
fastqc "${OUTDIR}/trimmed_fastq/${RAW_PREFIX}_"*".trimmed.fastq.gz" \
	 --outdir="${OUTDIR}/fastqc" --noextract --nogroup --quiet --threads "${SLURM_CPUS_PER_TASK}"
echo 'DONE'
echo

echo "Section finished in "$((${SECONDS}-${SECTION_START}))" seconds"
echo

#################################################
# STAR alignment to genome index
#################################################
SECTION_START="${SECONDS}"

# make output folder
mkdir -p "${OUTDIR}/aligned_bam"

# load the STAR genome index
STAR --genomeLoad LoadAndExit --genomeDir "${GENOMEDIR}" \
	--outFileNamePrefix "${OUTDIR}/aligned_bam/LoadGenome_"
echo
#
 loop through the trimmed fastq files
for PROCESSED_FASTQ in "${OUTDIR}/trimmed_fastq/${RAW_PREFIX}_"*".trimmed.fastq.gz"; do
	echo "Aligning using STAR: $(basename ${PROCESSED_FASTQ})"
	# run STAR aligner
	# the outFilterMismatch parameters are changed from defaults to enforce more strict alignments
	# will only output uniquely mapping reads
	STAR --runThreadN "${SLURM_CPUS_PER_TASK}" \
		--genomeLoad LoadAndKeep \
		--genomeDir "${GENOMEDIR}" \
		--readFilesIn "${PROCESSED_FASTQ}" \
		--readFilesCommand zcat \
		--limitBAMsortRAM 10000000000 \
		--outFileNamePrefix "${OUTDIR}/aligned_bam/$(basename ${PROCESSED_FASTQ} .trimmed.fastq.gz)_" \
		--outSAMtype BAM SortedByCoordinate \
		--outFilterMismatchNoverLmax 0.1 \
		--outFilterMismatchNmax 999 \
		--outFilterMultimapNmax 1
	echo
done

# unload the STAR genome index
STAR --genomeLoad Remove --genomeDir "${GENOMEDIR}" \
	--outFileNamePrefix "${OUTDIR}/aligned_bam/RemoveGenome_"
echo

# once finished aligning, loop through output files to cleanup
for STAR_BAM_FILE in ${OUTDIR}/aligned_bam/${RAW_PREFIX}_*_Aligned.sortedByCoord.out.bam; do
	# check last line of logfile to determine that STAR has finished succesfully
	if [[ $(tail -1 "${STAR_BAM_FILE%Aligned.sortedByCoord.out.bam}Log.out") == 'ALL DONE!' ]]; then
		echo "Alignment completed for: $(basename ${STAR_BAM_FILE})"
		# delete matching empty temporary directories recursively (STAR is not good about cleaning up)
		echo 'Cleaning up empty temporary directories...'
		find "${STAR_BAM_FILE%Aligned.sortedByCoord.out.bam}_STARtmp" -type d -empty -delete
		# check second to last line of log to check if empty (STAR will output a .bam file even if no reads align)
		if [[ $(tail -2 "${STAR_BAM_FILE%Aligned.sortedByCoord.out.bam}Log.out" | head -1) \
			== 'WARNING: nothing to sort - no output alignments' ]]; then
			# delete the empty .bam file
			echo 'bam file is empty, removing...'
			rm "${STAR_BAM_FILE}"
		fi
	fi
done
echo 'DONE'
echo

echo "Section finished in "$((${SECONDS}-${SECTION_START}))" seconds"
echo

#################################################
# SAMTOOLS indexing and statistics
#################################################
SECTION_START="${SECONDS}"

# make output folder
mkdir -p "${OUTDIR}/samtools"

# validate remaining output files with samtools
echo 'Validating bam files...'
samtools quickcheck -v "${OUTDIR}/aligned_bam/${RAW_PREFIX}_"*"_Aligned.sortedByCoord.out.bam"
echo 'DONE'
echo

# loop through bam files
for STAR_BAM_FILE in "${OUTDIR}/aligned_bam/${RAW_PREFIX}_"*"_Aligned.sortedByCoord.out.bam"; do
	echo "Processing bam file: $(basename ${STAR_BAM_FILE})"
	# index with samtools (required for deduplication)
	echo "Indexing using Samtools..."
	samtools index "${STAR_BAM_FILE}" -@ $((${SLURM_CPUS_PER_TASK}-1))
	# generate initial qc statistics with samtools stats
	echo "Generating statistics using Samtools..."
	samtools stats "${STAR_BAM_FILE}" -@ $((${SLURM_CPUS_PER_TASK}-1)) \
		> "${OUTDIR}/samtools/$(basename ${STAR_BAM_FILE} _Aligned.sortedByCoord.out.bam).stats.aligned.out"
done
echo 'DONE'
echo

echo "Section finished in "$((${SECONDS}-${SECTION_START}))" seconds"
echo

#################################################
# UMI_TOOLS deduplication of bam files
#################################################
SECTION_START="${SECONDS}"

# make output folder
mkdir -p "${OUTDIR}/deduplicated_bam"

# deduplicate UMIs using UMI-tools
# uses xargs -P to run multiple parallel processes and reduce runtime
find "${OUTDIR}/aligned_bam" -maxdepth 1 -type f -name "${RAW_PREFIX}_*_Aligned.sortedByCoord.out.bam" -print0 | \
	xargs -0 -I '{}' -P "${SLURM_CPUS_PER_TASK}" bash -c \
	'echo "Using UMI-tools to deduplicate: "$(basename $1)""; \
	umi_tools dedup --stdin="$1" \
	--log ""$2"/deduplicated_bam/"$(basename $1 _Aligned.sortedByCoord.out.bam)".dedup.log" \
	--output-stats=""$2"/deduplicated_bam/"$(basename $1 _Aligned.sortedByCoord.out.bam)"" \
	> ""$2"/deduplicated_bam/"$(basename $1 _Aligned.sortedByCoord.out.bam)".dedup.bam"' \
	-- '{}' "${OUTDIR}"
echo 'DONE'
echo

echo "Section finished in "$((${SECONDS}-${SECTION_START}))" seconds"
echo

#################################################
# SAMTOOLS indexing and statistics of deduplicated files
#################################################
SECTION_START="${SECONDS}"

# make output folder
mkdir -p "${OUTDIR}/samtools"

# validate remaining output files with samtools
echo 'Validating deduplicated bam files...'
samtools quickcheck -v "${OUTDIR}/deduplicated_bam/${RAW_PREFIX}_"*".dedup.bam"
echo 'DONE'
echo

# loop through deduplicated bam files
for DEDUP_BAM_FILE in "${OUTDIR}/deduplicated_bam/${RAW_PREFIX}_"*".dedup.bam"; do
	echo "Processing deduplicated bam file: $(basename ${DEDUP_BAM_FILE})"
	# index the deduplicated file (required for downstream QC)
	echo "Indexing using Samtools..."
	samtools index "${DEDUP_BAM_FILE}" -@ $((${SLURM_CPUS_PER_TASK}-1))
	# generate qc statistics using samtools
	echo "Generating statistics using Samtools..."
	samtools stats "${DEDUP_BAM_FILE}" -@ $((${SLURM_CPUS_PER_TASK}-1)) \
		> "${OUTDIR}/samtools/$(basename ${DEDUP_BAM_FILE} .dedup.bam).stats.dedup.out"
done
echo 'DONE'
echo

echo "Section finished in "$((${SECONDS}-${SECTION_START}))" seconds"
echo

#################################################
# RSEQC analysis of deduplicated files
#################################################
SECTION_START="${SECONDS}"

# make output folder
mkdir -p "${OUTDIR}/rseqc"

# run RSEQC scripts to generate statistics of the deduplicated files
# uses xargs -P to run multiple parallel processes and reduce runtime
find "${OUTDIR}/deduplicated_bam" -maxdepth 1 -type f -name "${RAW_PREFIX}_*.dedup.bam" -print0 | \
	xargs -0 -I '{}' -P "${SLURM_CPUS_PER_TASK}" bash -c \
	'echo "Running RSEQC scripts on: "$(basename $1)""; \
	""$2"/bam_stat.py" -i "$1" > ""$4"/rseqc/"$(basename $1 .bam)".bamStat.txt" 2> /dev/null; \
	""$2"/infer_experiment.py" -i "$1" -r "$3" > ""$4"/rseqc/"$(basename $1 .bam)".inferExperiment.txt" 2> /dev/null; \
	""$2"/read_distribution.py" -i "$1" -r "$3" > ""$4"/rseqc/"$(basename $1 .bam)".readDistribution.txt" 2> /dev/null' \
	-- '{}' "${RSEQC_SCRIPTS}" "${BEDFILE}" "${OUTDIR}"
echo 'DONE'
echo

echo "Section finished in "$((${SECONDS}-${SECTION_START}))" seconds"
echo

#################################################
# Qualimap analysis of deduplicated files
#################################################
SECTION_START="${SECONDS}"

# run Qualimap to generate statistics of the deduplicated files
# uses xargs -P to run multiple parallel processes and reduce runtime
# make sure that the (--java-mem-size) multiplied by the number of processes (-P) does not exceed available memory
find "${OUTDIR}/deduplicated_bam" -maxdepth 1 -type f -name "${RAW_PREFIX}_*.dedup.bam" -print0 | \
	xargs -0 -I '{}' -P 4 bash -c \
	'echo "Running Qualimap on: "$(basename $1)""; \
	qualimap rnaseq --java-mem-size=12G -bam "$1" -gtf "$2" -outdir ""$3"/qualimap/"$(basename $1 .dedup.bam)"_dedup" > "/dev/null" 2>&1' \
	-- '{}' "${GTFFILE}" "${OUTDIR}"
echo 'DONE'
echo

echo "Section finished in "$((${SECONDS}-${SECTION_START}))" seconds"
echo


# report the total time
SCRIPT_DURATION=$SECONDS
echo "Script completed at $(date)"
echo "Total script duration was: $(printf '%02dh:%02dm:%02ds\n' $(($SCRIPT_DURATION/3600)) $(($SCRIPT_DURATION%3600/60)) $(($SCRIPT_DURATION%60)))"


##### END OF SCRIPT #####
