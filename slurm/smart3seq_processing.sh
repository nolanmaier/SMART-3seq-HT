#! /bin/bash

#################################################
# SMART3SEQ PROCESSING SCRIPT
# this script processes input fastq.gz files 
# will demultiplex, trim, align, deduplicate, and peform some QC
#################################################


#################################################
# SETUP
#################################################

# tell bash to exit immediately if something in the code fails
set -e
set -o pipefail

# define a function to display help information for the script
Help()
{
  echo 'Syntax: smart3seq_processing.sh -i [input directory] -o [output directory] -t [TSO barcode fasta file] -s [STAR genome index] -g [gtf reference file] -b [bed reference file]'
}

# get the input options from the command line arguments
while getopts ":i:o:t:s:g:b:h" flag; do
  case "${flag}" in
  	i) indir=${OPTARG} ;;
    o) outdir=${OPTARG} ;;
		t) tso_path=${OPTARG} ;;
		s) genomedir=${OPTARG} ;;
    g) gtffile=${OPTARG} ;;
		b) bedfile=${OPTARG} ;;
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

cores_avail="${SLURM_CPUS_PER_TASK}"

# load biogrids module to access analysis software
# tools used: fqtools, cutadapt, fastqc, STAR, samtools, umi_tools, rseqc, qualimap
module load biogrids/latest
# harcode the location of the RSeQC scripts from the biogrids install
rseqc_scripts='/programs/x86_64-linux/rseqc/4.0.0/bin.capsules'
echo

# print versions of software used for documentation
echo 'Programs used:'
echo "Cutadapt version: $(cutadapt --version)"
fastqc -v
echo "STAR version: $(STAR --version)"
samtools --version | head -n 1
umi_tools -v
echo "RSEQC version: $(${rseqc_scripts}/bam_stat.py --version) $(${rseqc_scripts}/infer_experiment.py --version) $(${rseqc_scripts}/read_distribution.py --version)"
qualimap --help 2> '/dev/null' | sed -n '/^QualiMap/p'
echo

# each job of the array will run one fastq file found in the input directory, sorted to match array id with filename if possible
raw_fastq=$(find "${indir}" -maxdepth 1 -type f -size +0 -name '*.fastq.gz' -print0 | sort -z |sed -nz "${SLURM_ARRAY_TASK_ID}"p)
# define the filename prefix
raw_prefix="$(basename ${raw_fastq} .fastq.gz)"
echo 'This script will process:'
echo ${raw_fastq}
echo


#################################################
# CUTADAPT demultiplexing using barcodes file
#################################################
section_start="${SECONDS}"

# before starting check the input fastq file using fqtools
echo "Validating: $(basename ${raw_fastq})"
cutadapt --report=minimal -o /dev/null "${raw_fastq}"
echo

# make output folder
mkdir -p "${outdir}/demultiplexed_fastq"

# find the barcode file to demultiplex with
dT_regex='dT[0-9][0-9]'
[[ $raw_prefix =~ $dT_regex ]]
sample_number="${BASH_REMATCH[0]}"
tso_file="${tso_path}/${sample_number}_barcodes.fasta"
echo "Using barcode file: $(basename ${tso_file})"

# run cutadapt demultiplexing
echo "Demultiplexing: $(basename ${raw_fastq})"
cutadapt --cores="${cores_avail}" -g ^file:"${tso_file}" -e 1 --no-indels \
	--untrimmed-output "${outdir}/demultiplexed_fastq/${raw_prefix}_TSONaN_OFF0.fastq.gz" \
	-o "${outdir}/demultiplexed_fastq/${raw_prefix}_{name}.fastq.gz" \
	"${raw_fastq}" \
	> "${outdir}/demultiplexed_fastq/${raw_prefix}.cutadapt_demulti.log"
echo "Demultiplexed into" $(ls -1q ${outdir}/demultiplexed_fastq/${raw_prefix}_*.fastq.gz | wc -l) "files"
echo 'DONE'
echo

echo "Section finished in "$((${SECONDS}-${section_start}))" seconds"
echo


#################################################
# CUTADAPT trimming and UMI extraction
#################################################
section_start="${SECONDS}"

# make output folder
mkdir -p "${outdir}/trimmed_fastq"

# iterate through the raw demultiplexed fastq files
for demulti_fastq in "${outdir}/demultiplexed_fastq/${raw_prefix}_"*".fastq.gz"; do
	echo "Using Cutadapt to trim: $(basename ${demulti_fastq})"
	
	# find offset using filename
	offset_regex='OFF([0-9])'
	[[ $demulti_fastq =~ $offset_regex ]]
	offset="${BASH_REMATCH[1]}"
	# if the offset is 3, remove 3 bases from the 5' end (using +3).
	# if the offset is 0, remove 3 bases from the 3' end (using -3) to retain identical read lengths.
	if [ "${offset}" -eq 0 ]; then
		offset=-3
	fi

	# cutadapt runs multiple times on the same file:
	# 1. trim offset (no logs or intermediate files are kept)
	# 2. trim UMI and append UMI sequence to read name (no logs or intermediate files are kept)
	# 3a. trim 13 bases from the 5' end (corresponding to the TSO constant sequence)
	# 3b. trim any low qulity bases from the 3' end (assuming two-color chemistry)
	# 3c. trim any remaining N bases from either end
	# Finally output a new fastq file and a log file, only keep reads with at least 1bp
	cutadapt --cores="${cores_avail}" --quiet --cut "${offset}" "${demulti_fastq}" | \
	cutadapt --cores="${cores_avail}" --quiet --cut 5 --rename '{id}_{cut_prefix} {comment}' - | \
	cutadapt --cores="${cores_avail}" --cut 13 --nextseq-trim 15 --adapter "illumina=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=5" \
	2> "${outdir}/trimmed_fastq/$(basename ${demulti_fastq} .fastq.gz).cutadapt_illumina.log" - | \
	cutadapt --cores="${cores_avail}" --adapter "poly_a=A{100};min_overlap=5" --trim-n --minimum-length 1 \
	-o "${outdir}/trimmed_fastq/$(basename ${demulti_fastq} .fastq.gz).trimmed.fastq.gz" - \
	> "${outdir}/trimmed_fastq/$(basename ${demulti_fastq} .fastq.gz).cutadapt.log"

	# validate the output file
	echo 'Validating trimmed file:'
	cutadapt --quiet -o /dev/null "${outdir}/trimmed_fastq/$(basename ${demulti_fastq} .fastq.gz).trimmed.fastq.gz" && echo 'OK'

done
echo 'DONE'
echo

echo "Section finished in "$((${SECONDS}-${section_start}))" seconds"
echo


#################################################
# FASTQC on trimmed files
#################################################
section_start="${SECONDS}"

# make output folder
mkdir -p "${outdir}/fastqc"

# run fastqc
echo 'Running FastQC on processed fastq files...'
fastqc "${outdir}/trimmed_fastq/${raw_prefix}_"*".trimmed.fastq.gz" \
	 --outdir="${outdir}/fastqc" --quiet --threads "${cores_avail}"
echo 'DONE'
echo

echo "Section finished in "$((${SECONDS}-${section_start}))" seconds"
echo


#################################################
# STAR alignment to genome index
#################################################
section_start="${SECONDS}"

# make output folder
mkdir -p "${outdir}/aligned_bam"

# load the STAR genome index
STAR --genomeLoad LoadAndExit --genomeDir "${genomedir}" \
	--outFileNamePrefix "${outdir}/aligned_bam/LoadGenome_"
echo

# loop through the trimmed fastq files
for processed_fastq in "${outdir}/trimmed_fastq/${raw_prefix}_"*".trimmed.fastq.gz"; do
	echo "Aligning using STAR: $(basename ${processed_fastq})"
	# run STAR aligner
	# the outFilterMismatch parameters are changed from defaults to enforce more strict alignments
	# will only output uniquely mapping reads
	STAR --runThreadN "${cores_avail}" \
		--genomeLoad LoadAndKeep \
		--genomeDir "${genomedir}" \
		--readFilesIn "${processed_fastq}" \
		--readFilesCommand zcat \
		--limitBAMsortRAM 10000000000 \
		--outFileNamePrefix "${outdir}/aligned_bam/$(basename ${processed_fastq} .trimmed.fastq.gz)_" \
		--outSAMtype BAM SortedByCoordinate \
		--outFilterMismatchNoverLmax 0.1 \
		--outFilterMismatchNmax 999 \
		--outFilterMultimapNmax 1
	echo
done

# unload the STAR genome index
STAR --genomeLoad Remove --genomeDir "${genomedir}" \
	--outFileNamePrefix "${outdir}/aligned_bam/RemoveGenome_"
echo

# once finished aligning, loop through output files to cleanup
for star_bam_file in ${outdir}/aligned_bam/${raw_prefix}_*_Aligned.sortedByCoord.out.bam; do
	# check last line of logfile to determine that STAR has finished succesfully
	if [[ $(tail -1 "${star_bam_file%Aligned.sortedByCoord.out.bam}Log.out") == 'ALL DONE!' ]]; then
		echo "Alignment completed for: $(basename ${star_bam_file})"
		# delete matching empty temporary directories recursively (STAR is not good about cleaning up)
		echo 'Cleaning up empty temporary directories...'
		find "${star_bam_file%Aligned.sortedByCoord.out.bam}_STARtmp" -type d -empty -delete
		# check second to last line of log to check if empty (STAR will output a .bam file even if no reads align)
		if [[ $(tail -2 "${star_bam_file%Aligned.sortedByCoord.out.bam}Log.out" | head -1) \
			== 'WARNING: nothing to sort - no output alignments' ]]; then
			# delete the empty .bam file
			echo 'bam file is empty, removing...'
			rm "${star_bam_file}"
		fi
	fi
done
echo 'DONE'
echo

echo "Section finished in "$((${SECONDS}-${section_start}))" seconds"
echo


#################################################
# SAMTOOLS indexing
#################################################
section_start="${SECONDS}"

# validate remaining output files with samtools
echo 'Validating bam files...'
samtools quickcheck -v "${outdir}/aligned_bam/${raw_prefix}_"*"_Aligned.sortedByCoord.out.bam"
echo 'DONE'
echo

# loop through bam files
for star_bam_file in "${outdir}/aligned_bam/${raw_prefix}_"*"_Aligned.sortedByCoord.out.bam"; do
	# index with samtools (required for deduplication)
	echo "Indexing $(basename ${star_bam_file})"
	samtools index "${star_bam_file}" -@ $((${cores_avail}-1))
done
echo 'DONE'
echo

echo "Section finished in "$((${SECONDS}-${section_start}))" seconds"
echo


#################################################
# UMI_TOOLS deduplication of bam files
#################################################
section_start="${SECONDS}"

# make output folder
mkdir -p "${outdir}/deduplicated_bam"

# deduplicate UMIs using UMI-tools
# uses xargs -P to run multiple parallel processes and reduce runtime
find "${outdir}/aligned_bam" -maxdepth 1 -type f -name "${raw_prefix}_*_Aligned.sortedByCoord.out.bam" -print0 | \
	xargs -0 -I '{}' -P "${cores_avail}" bash -c \
	'echo "Using UMI-tools to deduplicate: "$(basename $1)""; \
	umi_tools dedup --stdin="$1" \
	--log ""$2"/deduplicated_bam/"$(basename $1 _Aligned.sortedByCoord.out.bam)".umitools.log" \
	--output-stats=""$2"/deduplicated_bam/"$(basename $1 _Aligned.sortedByCoord.out.bam)"" \
	> ""$2"/deduplicated_bam/"$(basename $1 _Aligned.sortedByCoord.out.bam)".dedup.bam"' \
	-- '{}' "${outdir}"
echo 'DONE'
echo

echo "Section finished in "$((${SECONDS}-${section_start}))" seconds"
echo


#################################################
# SAMTOOLS indexing of deduplicated files
#################################################
section_start="${SECONDS}"

# validate remaining output files with samtools
echo 'Validating deduplicated bam files...'
samtools quickcheck -v "${outdir}/deduplicated_bam/${raw_prefix}_"*".dedup.bam"
echo 'DONE'
echo

# loop through deduplicated bam files
for dedup_bam_file in "${outdir}/deduplicated_bam/${raw_prefix}_"*".dedup.bam"; do
	# index the deduplicated file (required for downstream QC)
	echo "Indexing $(basename ${dedup_bam_file})"
	samtools index "${dedup_bam_file}" -@ $((${cores_avail}-1))
done
echo 'DONE'
echo

echo "Section finished in "$((${SECONDS}-${section_start}))" seconds"
echo


#################################################
# RSEQC analysis of deduplicated files
#################################################
section_start="${SECONDS}"

# make output folder
mkdir -p "${outdir}/rseqc"

# run RSEQC scripts to generate statistics of the deduplicated files
# uses xargs -P to run multiple parallel processes and reduce runtime
find "${outdir}/deduplicated_bam" -maxdepth 1 -type f -name "${raw_prefix}_*.dedup.bam" -print0 | \
	xargs -0 -I '{}' -P "${cores_avail}" bash -c \
	'echo "Running RSEQC scripts on: "$(basename $1)""; \
	""$2"/bam_stat.py" -i "$1" > ""$4"/rseqc/"$(basename $1 .bam)".bamStat.txt" 2> /dev/null; \
	""$2"/infer_experiment.py" -i "$1" -r "$3" > ""$4"/rseqc/"$(basename $1 .bam)".inferExperiment.txt" 2> /dev/null; \
	""$2"/read_distribution.py" -i "$1" -r "$3" > ""$4"/rseqc/"$(basename $1 .bam)".readDistribution.txt" 2> /dev/null' \
	-- '{}' "${rseqc_scripts}" "${bedfile}" "${outdir}"
echo 'DONE'
echo

echo "Section finished in "$((${SECONDS}-${section_start}))" seconds"
echo


#################################################
# Qualimap analysis of deduplicated files
#################################################
section_start="${SECONDS}"

# run Qualimap to generate statistics of the deduplicated files
# uses xargs -P to run multiple parallel processes and reduce runtime
# make sure that the (--java-mem-size) multiplied by the number of processes (-P) does not exceed available memory
find "${outdir}/deduplicated_bam" -maxdepth 1 -type f -name "${raw_prefix}_*.dedup.bam" -print0 | \
	xargs -0 -I '{}' -P 4 bash -c \
	'echo "Running Qualimap on: "$(basename $1)""; \
	qualimap rnaseq --java-mem-size=12G -bam "$1" -gtf "$2" -outdir ""$3"/qualimap/"$(basename $1 .dedup.bam)"_dedup" > "/dev/null" 2>&1' \
	-- '{}' "${gtffile}" "${outdir}"
echo 'DONE'
echo

echo "Section finished in "$((${SECONDS}-${section_start}))" seconds"
echo


# report the total time
script_duration=$SECONDS
echo "Script completed at $(date)"
echo "Total script duration was: $(printf '%02dh:%02dm:%02ds\n' $(($script_duration/3600)) $(($script_duration%3600/60)) $(($script_duration%60)))"

##### END OF SCRIPT #####