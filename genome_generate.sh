#! /bin/bash

#SBATCH -c 10 							# cores requested
#SBATCH -t 0-04:00						# max time requested
#SBATCH -p short						# slurm partition used
#SBATCH --mem=40G						# memory requested
#SBATCH -o genome_generate_%j.out		# name the slurm log like this
#SBATCH --mail-type=FAIL				# send an email on job failure

#################################################
# SMART3SEQ GENOME GENERATION SCRIPT
# this script makes genome annotation files for SMART3-SEQ-HT
#################################################


#################################################
# SETUP
#################################################

# tell bash to exit immediately if something in the code fails
###set -e
###set -o pipefail

# define a function to display help information for the script
Help()
{
	echo 'Generates genome annotation files for SMART-3SEQ-HT'
	echo 'Will produce a GTF file, BED file, and STAR index'
	echo
	echo 'Required parameters:'
	echo '-f   fasta genome file (may be gzip compressed)'
	echo '-g   GFF3 annotation file (may be gzip compressed)'
	echo '-o   output directory'
	echo
	echo 'Syntax: sbatch genome_generate.sh -i [input_dir] -g [gff_file]'
	echo
}

# get the input options from the command line arguments
while getopts ":f:g:o:h" flag; do
	case "${flag}" in
		f) fasta_file=${OPTARG} ;;
		g) gff_file=${OPTARG} ;;
		o) out_dir=${OPTARG} ;;
		h) Help; exit ;;
		:) echo "$0: Must supply an argument to -${OPTARG}"; exit 1 ;;
		?) echo "Invalid option: -${OPTARG}"; Help; exit 1 ;;
	esac
done
shift "$(($OPTIND -1))"

# "SECONDS" is a special bash variable that tracks the time
SECONDS=0

# get the number of cores
cores_avail="${SLURM_CPUS_PER_TASK}"

# convert the obtained paths to absolute paths
fasta_file=$(realpath "${fasta_file}")
gff_file=$(realpath "${gff_file}")
out_dir=$(realpath "${out_dir}")

# source the biogrids environment
# and explicitly set the versions of all tools used
source /programs/biogrids.shrc
export AGAT_X=1.2.0
export STAR_X=2.7.9a
echo

#################################################
# MAIN PROGRAM
#################################################
echo "Starting genome generation at $(date)"
echo
echo -e "Using genome FASTA file: \t ${fasta_file}"
echo -e "Using GFF file: \t ${gff_file}"
echo

# get the base file name
out_prefix="$(basename ${gff_file} .gz)"
out_prefix="${out_prefix%.gff3}"
mkdir -p "${out_dir}"

# make a .gtf annotatgtfion file from the input .gff file
gtf_file="${out_dir}/${out_prefix}.gtf"
echo -e "Making GTF file: \t ${gtf_file}"
agat_convert_sp_gff2gtf.pl --gff "${gff_file}" --out "${gtf_file}"
echo

# make a .bed annotation file from the input .gff file
bed_file="${out_dir}/${out_prefix}.bed"
echo -e "Making BED file: \t ${bed_file}"
agat_convert_sp_gff2bed.pl --gff "${gff_file}" --out "${bed_file}"
echo

# make a STAR index from the chromosome fasta files and the .gff file
genome_dir="${out_dir}/${out_prefix}_STAR"
echo -e "Making STAR index at: \t ${genome_dir}"
STAR --runThreadN "${cores_avail}" --runMode genomeGenerate \
--genomeDir "${genome_dir}" --genomeFastaFiles <(zcat -f "${fasta_file}") \
--sjdbGTFfile <(zcat -f "${gff_file}") -sjdbGTFtagExonParentTranscript Parent
echo

# report the time
# report the total time
script_duration=$SECONDS
echo "Script completed at $(date)"
echo "Total script duration was: $(printf '%02dh:%02dm:%02ds\n' $(($script_duration/3600)) $(($script_duration%3600/60)) $(($script_duration%60)))"

##### END OF SCRIPT #####