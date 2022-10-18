#! /bin/bash

#################################################
# SMART3SEQ SUMMARY SCRIPT
# this script runs summary tools after all the individual fastq files have been processed
#################################################


#################################################
# SETUP
#################################################

# define a function to display help information for the script
Help()
{
  echo 'Syntax: smart3seq_summary.sh -o [output directory] -g [gtf reference file]'
}

# get the input options from the command line arguments
while getopts ":o:g:h" flag; do
  case "${flag}" in
    o) OUTDIR=${OPTARG} ;;
    g) GTFFILE=${OPTARG} ;;
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

# load biogrids module to access analysis software
# tools used: featurecounts, multiQC
module load biogrids/latest
echo

# print versions of tools used for documentation
echo 'Programs used:'
featureCounts -v
multiqc --version
echo


#################################################
# FEATURECOUNTS to summarize reads by gene
#################################################

# make output folder
mkdir -p "${OUTDIR}/featurecounts"

echo 'Summarizing genecounts from:'
echo "${OUTDIR}/deduplicated_bam"
echo

# run featurecounts to summarize exon counts for each gene
# note that "exon" encompasses 5'UTR, CDS, and 3' UTR
# read must match the correct strand of the genefin
# will use all deduplicated bam files (discard those with no reads and those without indices)
featureCounts -T ${SLURM_CPUS_PER_TASK} -t exon -g gene_id -s 1 -a "${GTFFILE}" \
  -o "${OUTDIR}/featurecounts/genecounts.txt" \
  $(find "${OUTDIR}/deduplicated_bam" -maxdepth 1 -type f -size +0 -name "*.dedup.bam" ! -name "*TSONaN*" ! -name "*Undetermined*")
echo


#################################################
# MULTIQC to generate aggregate report on pipeline
#################################################

# make output folder
mkdir -p "${OUTDIR}/multiqc"

# make custom multiqc config file
CUSTOM_CONFIG="${OUTDIR}/multiqc/multiqc_config.yaml"
cat >"${CUSTOM_CONFIG}" <<'EOF'
fn_ignore_files:
  - "LoadGenome*"
  - "RemoveGenome*"
  - "*TSONaN*"
  - "*Undetermined*"
sample_names_ignore:
  - "*TSONaN*"
  - "*Undetermined*"
module_order:
  - cutadapt:
      name: "Cutadapt (demultiplexing)"
      anchor: "cutadapt_demultiplexing"
      info: "This section of the report shows the demultiplexing results from cutadapt"
      target: ""
      path_filters:
        - "*.cutadapt_demulti.log"
  - fastqc:
      name: "FastQC (raw)"
      anchor: "fastqc_raw"
      info: "This section of the report shows FastQC results before trimming."
      target: ""
      path_filters:
        - "*.raw_fastqc.zip"
  - cutadapt:
      name: "Cutadapt (trimming)"
      anchor: "cutadapt_trimming"
      info: "This section of the report shows the trimming results from cutadapt"
      target: ""
      path_filters:
        - "*.cutadapt.log"
  - fastqc:
      name: "FastQC (trimmed)"
      anchor: "fastqc_trimmed"
      info: "This section of the report shows FastQC results after trimming."
      target: ""
      path_filters:
        - "*.trimmed_fastqc.zip"
  - star
  - samtools:
      name: "Samtools (aligned)"
      anchor: "samtools_aligned"
      info: "This section of the report shows samtools statistics before deduplication"
      target: ""
      path_filters:
        - "*.stats.aligned.out"
  - samtools:
      name: "Samtools (dedup)"
      anchor: "samtools_dedup"
      info: "This section of the report shows samtools statistics after deduplication"
      target: ""
      path_filters:
        - "*.stats.dedup.out"
  - rseqc
  - qualimap
  - featureCounts
remove_sections:
  - cutadapt_demultiplexing_cutadapt_trimmed_sequences
table_columns_visible:
  Cutadapt (demultiplexing):
    percent_trimmed: False
  Samtools (aligned):
    error_rate: False
    reads_mapped_percent: False
    raw_total_sequences: False
    non-primary_alignments: False
    reads_mapped: False
  Samtools (dedup):
    error_rate: False
    reads_mapped_percent: False
    raw_total_sequences: False
    non-primary_alignments: False
    reads_mapped: False
  FastQC (raw):
    percent_duplicates: False
    percent_gc: False
  FastQC (trimmed):
    avg_sequence_length: True
    percent_duplicates: False
  QualiMap:
    reads_aligned: False
table_columns_placement:
  FastQC (raw):
    total_sequences: 100
  Cutadapt (trimming):
    percent_trimmed: 200
  FastQC (trimmed):
    total_sequences: 300
    avg_sequence_length: 310
    percent_gc: 320
  STAR:
    uniquely_mapped: 400
    uniquely_mapped_percent: 410
  featureCounts:
    Assigned: 500
    percent_assigned: 510
  QualiMap:
    5_3_bias: 600
table_columns_name:
  FastQC (raw):
    total_sequences: "Raw Reads"
  FastQC (trimmed):
    total_sequences: "Trimmed Reads"
EOF

# run multiqc on output directory, overwrite previous reports
multiqc "${OUTDIR}" -o "${OUTDIR}/multiqc" -c "${CUSTOM_CONFIG}" -f -z --interactive


# report the time
echo
SCRIPT_DURATION=$SECONDS
echo "Script completed at $(date)"
echo "Total script duration was: $(printf '%02dh:%02dm:%02ds\n' $(($SCRIPT_DURATION/3600)) $(($SCRIPT_DURATION%3600/60)) $(($SCRIPT_DURATION%60)))"

##### END OF SCRIPT #####