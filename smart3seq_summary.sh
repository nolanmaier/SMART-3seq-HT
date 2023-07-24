#! /bin/bash

#################################################
# SMART3SEQ SUMMARY SCRIPT
# this script runs summary tools after all the individual fastq files have been processed
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
  echo 'Syntax: smart3seq_summary.sh -o [output directory] -g [gtf reference file]'
}

# get the input options from the command line arguments
while getopts ":o:g:h" flag; do
  case "${flag}" in
    o) outdir=${OPTARG} ;;
    g) gtffile=${OPTARG} ;;
    h) Help; exit ;;
    :) echo "$0: Must supply an argument to -${OPTARG}"; exit 1;;
    \?) echo "Invalid option: -${OPTARG}"; Help; exit 1 ;;
  esac
done
shift "$(($OPTIND -1))"


# "SECONDS" is a special bash variable that tracks the time
SECONDS=0

# print out the command submitted to initialize this script
echo "smart3seq summary script started at $(date)}"
echo 'Command Submitted:'
echo $(scontrol show job "${SLURM_JOBID}" | awk -F= '/Command=/{print $2;exit;}')
echo

cores_avail="${SLURM_CPUS_PER_TASK}"

# print versions of tools used for documentation
echo 'Programs used:'
featureCounts -v
multiqc --version
echo


#################################################
# FEATURECOUNTS to summarize reads by gene
#################################################
section_start="${SECONDS}"

# make output folder
mkdir -p "${outdir}/featurecounts"

echo 'Summarizing genecounts from:'
echo "${outdir}/deduplicated_bam"
echo

# run featurecounts to summarize exon counts for each gene
# note that "exon" encompasses 5'UTR, CDS, and 3' UTR
# read must match the correct strand of the genefin
# will use all deduplicated bam files (discard those with no reads and those without indices)
featureCounts -T "${cores_avail}" -t 'exon' -g 'gene_id' -s 1 -a "${gtffile}" \
  -o "${outdir}/featurecounts/genecounts.tsv" \
  $(find "${outdir}/deduplicated_bam" -maxdepth 1 -type f -size +0 -name "*.dedup.bam" ! -name "*TSONaN*" ! -name "*Undetermined*")
echo "featureCounts finished in "$((${SECONDS}-${section_start}))" seconds"
echo

#################################################
# MULTIQC to generate aggregate report on pipeline
#################################################
section_start="${SECONDS}"

# make output folder
mkdir -p "${outdir}/multiqc"

# make custom multiqc config file
custom_config="${outdir}/multiqc/multiqc_config.yaml"
echo "Writing custom MultiQC config file to: ${custom_config}"
cat >"${custom_config}" <<'EOF'
fn_ignore_files:
  - "LoadGenome*"
  - "RemoveGenome*"
  - "*TSONaN*"
  - "*Undetermined*"
sample_names_ignore:
  - "*TSONaN*"
  - "*Undetermined*"
run_modules:
  - cutadapt
  - fastqc
  - star
  - umitools
  - rseqc
  - qualimap
  - featureCounts
module_order:
  - cutadapt:
      name: "Cutadapt (demultiplexing)"
      anchor: "cutadapt_demultiplexing"
      info: "This section of the report shows the demultiplexing results from cutadapt"
      target: ""
      path_filters:
        - "*.cutadapt_demulti.log"
  - cutadapt:
      name: "Cutadapt (trimming)"
      anchor: "cutadapt_trimming"
      info: "This section of the report shows the trimming results from cutadapt"
      target: ""
      path_filters:
        - "*.cutadapt.log"
  - fastqc
  - star
  - umitools
  - rseqc
  - qualimap
  - featureCounts
sp:
  cutadapt:
    fn: '*.cutadapt*.log'
  umitools:
    fn: '*.umitools.log'
  rseqc/bam_stat:
    fn: '*.bamStat.txt'
    max_filesize: 500000
  rseqc/gene_body_coverage:
    skip: true
  rseqc/inner_distance:
    skip: true
  rseqc/junction_annotation:
    skip: true
  rseqc/junction_saturation:
    skip: true
  rseqc/read_gc:
    skip: true
  rseqc/read_distribution:
    fn: '*.readDistribution.txt'
    max_filesize: 500000
  rseqc/read_duplication_pos:
    skip: true
  rseqc/infer_experiment:
    fn: '*.inferExperiment.txt'
    max_filesize: 500000
  rseqc/tin:
    skip: true
remove_sections:
  - cutadapt_demultiplexing_cutadapt_trimmed_sequences
table_columns_visible:
  Cutadapt (demultiplexing):
    percent_trimmed: False
  FastQC:
    avg_sequence_length: False
    percent_duplicates: False
    median_sequence_length: True
  QualiMap:
    reads_aligned: False
table_columns_placement:
  Cutadapt (trimming):
    percent_trimmed: 200
  FastQC:
    total_sequences: 300
    median_sequence_length: 310
    percent_gc: 320
  STAR:
    uniquely_mapped: 400
    uniquely_mapped_percent: 410
  UMI-tools:
    output_reads: 500
    percent_passing_dedup: 510
  featureCounts:
    Assigned: 600
    percent_assigned: 610
  QualiMap:
    5_3_bias: 700
custom_plot_config:
  cutadapt_filtered_reads_plot:
    cpswitch_c_active: False
EOF
echo

# run multiqc on output directory, overwrite previous reports
multiqc "${outdir}" -o "${outdir}/multiqc" -c "${custom_config}" -f -z
echo
echo "MultiQC finished in "$((${SECONDS}-${section_start}))" seconds"
echo


# report the time
script_duration=$SECONDS
echo "Script completed at $(date)"
echo "Total script duration was: $(printf '%02dh:%02dm:%02ds\n' $(($script_duration/3600)) $(($script_duration%3600/60)) $(($script_duration%60)))"

##### END OF SCRIPT #####