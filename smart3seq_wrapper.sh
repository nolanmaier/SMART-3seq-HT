#! /bin/bash

#SBATCH -c 1 							# request 1 cores for job
#SBATCH -t 0-00:05						# job will run for max 5 minutes
#SBATCH -p short						# use the short partition
#SBATCH --mem=1G						# request 100 Mib memory for job
#SBATCH -o smart3seq_wrapper_%j.out		# name the slurm log like this
#SBATCH --mail-type=FAIL				# send an email on job failure

#################################################
# SMART3SEQ WRAPPER SCRIPT
# this script captures input from the user and submits additional scripts to process smart3seq data
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
	echo 'Processes smart3seq fastq files'
	echo 'Performs demultiplexing, trimming, alignment, UMI-deduplication, and some quality control checks'
	echo 'Designed for running on a slurm scheduling node'
	echo
	echo 'Required parameters:'
	echo '-i   input directory containing fastq.gz files'
	echo '-o   output directory to put the generated files'
	echo '-t   sample sheet in .csv format containing TSO barcodes'
	echo
	echo 'Syntax: sbatch smart3seq_wrapper.sh -i [input_dir] -o [output_dir] -t [barcode_file]'
	echo
}

# get the input options from the command line arguments
while getopts ":i:o:t:h" flag; do
	case "${flag}" in
		i) indir=${OPTARG} ;;
		o) outdir=${OPTARG} ;;
		t) tso_samplesheet=${OPTARG} ;;
		h) Help; exit ;;
		:) echo "$0: Must supply an argument to -${OPTARG}"; exit 1;;
		?) echo "Invalid option: -${OPTARG}"; Help; exit 1 ;;
	esac
done
shift "$(($OPTIND -1))"

# "SECONDS" is a special bash variable that tracks the time
SECONDS=0

# convert the obtained paths to absolute paths
indir=$(realpath "${indir}")
outdir=$(realpath "${outdir}")
tso_samplesheet=$(realpath "${tso_samplesheet}")

# find the path to the output TSO barcodes fasta files
tso_barcodes="${outdir}/TSO_barcodes"

# hardcode location of genome directory (generated previously using STAR genomeGenerate)
genomdir='/n/data1/hms/microbiology/jost/lab/genomes/human/gencode/star_index/'
# hardcode location of downloaded GTF reference file
gtffile='/n/data1/hms/microbiology/jost/lab/genomes/human/gencode/gencode.v38.annotation.gtf'
# hardcode location of downloaded BED reference file
bedfile='/n/data1/hms/microbiology/jost/lab/nolan/test/gencode.v38.annotation2.bed'

# find the source of the other scripts we need to run
# they should be placed in the same directory as this script
script_path=$(dirname $(realpath $(scontrol show job "${SLURM_JOBID}" | awk -F= '/Command=/{print $2}'| awk '{print $1}')))

# load the biogrids module
module load biogrids/latest
# explicitly set the versions of all tools used
# these modules and exported variables will be inherited by the downstream scripts
export CUTADAPT_X=4.4
export FASTQC_X=0.11.9
export STAR_X=2.7.9a
export SAMTOOLS_X=1.17
export UMITOOLS_X=1.1.2
export RSEQC_X=4.0.0
export QUALIMAP_X=2.2.1
export SUBREAD_X=2.0.3
export MULTIQC_X=1.14
echo

#################################################
# MAIN PROGRAM
#################################################

echo "Submitting smart3seq analysis job at $(date)"
echo
echo -e "Analyzing fastq.gz files in: \t ${indir}"
echo -e "Processed files will be placed in: \t ${outdir}"
echo -e "Demultiplexing using barcodes in: \t ${tso_samplesheet}"
echo -e "Aligning using genome directory in: \t ${genomdir}"
echo -e "Using GTF file in: \t ${gtffile}"
echo -e "Using BED file in: \t ${bedfile}"
echo -e "Using scripts located in: \t ${script_path}"
echo

# use the python script to convert the .csv sample sheet into .fasta files of TSO barcodes
# this should use the biogrids install of python loaded above, but any standard python3 install will work
echo 'generating fasta files for TSO barcodes'
echo -e "using python installed at: \t $(which python3)"
echo "$(python3 --version)"
python3 "${script_path}/barcode2fasta.py" "${tso_samplesheet}" "${tso_barcodes}"
echo 'fasta file generation completed'
echo


# make a directory to store the slurm logs
mkdir -p "${outdir}/slurmlogs"

# count the number of input fastq files
num_files=$(find "${indir}" -maxdepth 1 -type f -size +0 -name '*dT*.fastq.gz' -print0 | grep -zc .)
echo "${num_files} fastq.gz files found"

# submit the processing script and capture the jobid
echo "Submitting processing script as an array of ${num_files} jobs..."
echo 'Command submitted:'

# start xtrace mode to echo the next command before executing
set -x
processing_jobid=$(sbatch --parsable -p short -c 10 -t 0-03:00 --mem=50G --mail-type=FAIL,ARRAY_TASKS \
	--array=1-"${num_files}" -o "${outdir}/slurmlogs/smart3seq_processing_%A_%a.out" \
	"${script_path}/smart3seq_processing.sh" \
	-i "${indir}" -o "${outdir}" -t "${tso_barcodes}" -s "${genomdir}" -g "${gtffile}" -b "${bedfile}")
set +x
echo "Submitted batch job ${processing_jobid}"
echo

# submit the summary script to run after the processing script finishes
echo "Submitting summary script to begin after completion of jobid: ${processing_jobid}"
echo 'Command submitted:'

# start xtrace mode to echo the next command before executing
set -x
sbatch -p short -c 4 -t 0-01:00 --mem=40G --mail-type=FAIL,END \
	--dependency=afterok:"${processing_jobid}" -o "${outdir}/slurmlogs/smart3seq_summary_%j.out" \
	"${script_path}/smart3seq_summary.sh" \
	-o "${outdir}" -g "${gtffile}"
set +x
echo

# report the time
echo "smart3seq job submission completed in ${SECONDS} seconds"

##### END OF SCRIPT #####