#! /bin/bash

#SBATCH -c 1 							# request 1 cores for job
#SBATCH -t 0-00:05						# job will run for max 5 minutes
#SBATCH -p short						# use the short partition
#SBATCH --mem=100M						# request 100 Mib memory for job
#SBATCH -o smart3seq_wrapper_%j.out		# name the slurm log like this
#SBATCH --mail-type=FAIL				# send an email on job failure

#################################################
# SMART3SEQ WRAPPER SCRIPT
# this script captures input from the user and submits additional scripts to process smart3seq data
#################################################


#################################################
# SETUP
#################################################

# define a function to display help information for the script
Help()
{
	echo 'Processes smart3seq fastq files'
	echo 'Performs demultiplexing, trimming, alignment, UMI-deduplication, and some quality control checks'
	echo 'Designed for running on a slurm scheduling node'
	echo
	echo 'Required parameters:'
	echo '-i   input directory containing fastq.gz files'
	echo '-o   output directory to put the generated files'
	echo '-t   TSO barcode file (in fasta format) containing index names and sequences'
	echo
	echo 'Syntax: sbatch smart3seq_wrapper.sh -i [input_dir] -o [output_dir] -t [barcode_file]'
	echo
}

# get the input options from the command line arguments
while getopts ":i:o:t:h" flag; do
	case "${flag}" in
		i) INDIR=${OPTARG} ;;
		o) OUTDIR=${OPTARG} ;;
		t) TSO_FILE=${OPTARG} ;;
		h) Help; exit ;;
		:) echo "$0: Must supply an argument to -${OPTARG}"; exit 1;;
		?) echo "Invalid option: -${OPTARG}"; Help; exit 1 ;;
	esac
done
shift "$(($OPTIND -1))"

# convert the obtained paths to absolute paths
INDIR=$(realpath "${INDIR}")
OUTDIR=$(realpath "${OUTDIR}")
TSO_FILE=$(realpath "${TSO_FILE}")

# hardcode location of genome directory (generated previously using STAR genomeGenerate)
GENOMEDIR='/n/data1/hms/microbiology/jost/lab/genomes/human/gencode/star_index/'
# hardcode location of downloaded GTF reference file
GTFFILE='/n/data1/hms/microbiology/jost/lab/genomes/human/gencode/gencode.v38.annotation.gtf'
# hardcode location of downloaded BED reference file
BEDFILE='/n/data1/hms/microbiology/jost/lab/nolan/test/gencode.v38.annotation2.bed'

# find the source of the other scripts we need to run
# they should be placed in the same directory as this script
SCRIPT_PATH=$(dirname $(realpath $(scontrol show job "${SLURM_JOBID}" | awk -F= '/Command=/{print $2}'| awk '{print $1}')))


#################################################
# MAIN PROGRAM
#################################################

echo "Submitting smart3seq analysis job at $(date)"
echo
echo -e "Analyzing fastq.gz files in: \t ${INDIR}"
echo -e "Processed files will be placed in: \t ${OUTDIR}"
echo -e "Demultiplexing using barcodes in: \t ${TSO_FILE}"
echo -e "Aligning using genome directory in: \t ${GENOMEDIR}"
echo -e "Using GTF file in: \t ${GTFFILE}"
echo -e "Using BED file in: \t ${BEDFILE}"
echo -e "Using scripts located in: \t ${SCRIPT_PATH}"
echo

# make a directory to store the slurm logs
mkdir -p "${OUTDIR}/slurmlogs"

# count the number of input fastq files
NUM_FILES=$(find "${INDIR}" -maxdepth 1 -type f -size +0 -name '*.fastq.gz' -print0 | grep -zc .)
echo "${NUM_FILES} fastq.gz files found"

# submit the processing script and capture the jobid
echo "Submitting processing script as an array of ${NUM_FILES} jobs..."
echo 'Command submitted:'

# start xtrace mode to echo the next command before executing
set -x
PROCESSING_JOBID=$(sbatch --parsable -p short -c 10 -t 0-06:00 --mem=50G --mail-type=FAIL,ARRAY_TASKS \
	--array=1-"${NUM_FILES}" -o "${OUTDIR}/slurmlogs/smart3seq_processing_%A_%a.out" \
	"${SCRIPT_PATH}/smart3seq_processing.sh" \
	-i "${INDIR}" -o "${OUTDIR}" -t "${TSO_FILE}" -s "${GENOMEDIR}" -g "${GTFFILE}" -b "${BEDFILE}")
set +x
echo "Submitted batch job ${PROCESSING_JOBID}"
echo

# submit the summary script to run after the processing script finishes
echo "Submitting summary script to begin after completion of jobid: ${PROCESSING_JOBID}"
echo 'Command submitted:'

# start xtrace mode to echo the next command before executing
set -x
sbatch -p short -c 4 -t 0-01:00 --mem=10G --mail-type=FAIL,END \
	--dependency=afterok:"${PROCESSING_JOBID}" -o "${OUTDIR}/slurmlogs/smart3seq_summary_parallel_%j.out" \
	"${SCRIPT_PATH}/smart3seq_summary.sh" \
	-o "${OUTDIR}" -g "${GTFFILE}"
set +x
echo

echo 'smart3seq job submission complete'

##### END OF SCRIPT #####