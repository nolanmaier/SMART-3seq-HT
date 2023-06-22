# SMART-3seq-HT
This repository contains scripts to demultiplex and analyze SMART-3seq-HT sequencing reads.

### Quick-Start Guide (TL:DR version)
1. clone this repo to your folder: `git clone https://github.com/nolanmaier/SMART-3seq-HT`
2. make a sample demultiplexing *.csv* file and upload it to your folder
3. run the wrapper script using slurm: `sbatch SMART-3seq-HT/smart3seq_wrapper.sh -i path/to/input_fastq_dir -o path/to/output_dir -t sample_demultiplexing.csv`

### Requirements
This repository is designed and tested for use on the Harvard Medical School [O2 computing cluster](https://harvardmed.atlassian.net/wiki/spaces/O2/overview). 

As written, it minimally requires a [SLURM](https://slurm.schedmd.com/documentation.html) job scheduler and a [BioGrids](https://biogrids.org/) software stack available as a module.

Additionally, the 

### Inputs
The pipeline takes as input three 

### Sample Demultiplexing CSV
The pipeline requires a comma separated file containing all samples and their associated oligo-dT and TSO barcodes. Reference the example in the repository. I have not yet found a good way to generate this *.csv* file programatically. So I just do it by hand. 😓

The minimum columns required are: **dT_Index**, **TSO_Index**, **TSO_Offset**, **TSO_Sequence** (the column names must match exactly).

I usually include extra columns (that are not parsed by the pipeline) describing the samples (i.e treatment names, treatment concentrations, etc) so that I can use this same *.csv* file during data analysis to decode the samples into treatments and controls. It is important to **NOT** use commas within fields in your *.csv* or the *.csv* will not be read correctly. If you must use commas (such as for sample names) enclose the entire field in double quotes e.g.:	✅ "3,4-hexanediol"	❌ 3,4-hexanediol.

NM011_sampledecode.csv
Login to the HMS O2 HPC cluster if have not already (this is different from the transfer server we used above)
I like to use the OpenOnDemand portal: https://o2portal.rc.hms.harvard.edu/
From there you can click Clusters > O2 Cluster Terminal and input your O2 password
Request an interactive session and wait a moment until your session starts
The Research Computing Group does not want you running analyses on the login nodes
srun --pty -p interactive -t 0-04:00 /bin/bash
Make a directory on the O2 server to store the scripts associated with this processing run using the mkdir command (or you can use the OpenOnDemand web interface).
mkdir -pv /n/data1/hms/microbiology/jost/lab/nolan/rna_sequencing_data/NM011/20230519_processing/scripts
Navigate to the scripts directory for this experiment that you created above
cd /n/data1/hms/microbiology/jost/lab/nolan/rna_sequencing_data/NM011/20230519_processing/scripts
Upload the sample decode .csv file from above and the four scripts below that constitute the pipeline
I tend to use the Upload button on the OpenOnDemand web interface: https://o2portal.rc.hms.harvard.edu/
smart3seq_wrapper.sh
barcode2fasta.py
smart3seq_processing.sh
smart3seq_summary.sh
The pipeline requires genome/transcriptome annotation files for your organism of interest. The locations of these files are hardcoded in the smart3seq_wrapper.sh script. If you are running this script without access to the Jost lab O2 folder, or are analyzing non-human data, you will need to generate new annotation files and change the relevant file paths in this script.
Submit the smart3seq_wrapper.sh script as a slurm sbatch job pointing it to the correct subfolders
The script takes three arguments: 
-i an input directory containing the raw fastq.gz files transferred from Basespace
-o an output directory to store the (many) output files generated
-t the sample decoding .csv you made above containing the TSO barcodes
Here I have used absolute paths to the correct folders, storing the majority of the path in a variable
experiment_dir="/n/data1/hms/microbiology/jost/lab/nolan/rna_sequencing_data/NM011"
sbatch ${experiment_dir}/20230519_processing/scripts/smart3seq_wrapper.sh -i ${experiment_dir}/raw_fastq -o ${experiment_dir}/20230519_processing -t ${experiment_dir}/20230519_processing/scripts/NM011_sampledecode.csv
You can monitor the progress of the pipeline using the O2squeue command
You should originally see a short 5 min job queued. This is the smart3seq_wrapper.sh script which submits the other scripts. If you need to update the memory requirements or requested cpu time of the pipeline, you can edit the parameters in this script using the sbatch options.
Once the wrapper completes you will see one job in the queue for each fastq.gz file in your input directory, These jobs are all running the smart3seq_processing.sh script and will run in parallel to each other.
There will be one additional job in the queue marked as STATE: DEPENDENCY. This is the smart3seq_summary.sh script. This job will not run until all of the smart3seq_processing.sh jobs finish. Once the smart3seq_summary.sh script finishes, the pipeline is complete. You should receive an email notification upon successful pipeline completion or if any of the jobs fail.

Pipeline Output
The pipeline will create a bunch of sub-directories within the output directory you indicated above. These sub-directories correspond to the various steps in the pipeline. The pipeline uses analysis software assembled and maintained by BioGrids and pre-installed on the O2 cluster. If you are not running this on the O2 cluster, these programs will need to be installed and the scripts updated to reflect the new install locations.
TSO_barcodes
Demultiplexing of the oligo-dT barcodes was done on the instrument or in Basespace above yielding one fastq.gz file per oligo-dT. The read processing pipeline will demultiplex each of these files further according to the TSO barcodes yielding one fastq.gz file per oligo-dT/TSO pair. This step requires one .fasta file per oligo-dT, containing the sequences of the TSO barcodes that correspond to that oligo-dT. The barcode2fasta.py script uses the sample decoding .csv to generate these .fasta files. The .fasta files also encode information about the TSO offset which is used in the pipeline for read trimming.
slurmlogs
This directory contains the log files from the bash scripts above. If the pipeline fails, this is the first place to look for explanations. Even if the pipeline completes, it is useful to scan through at least one "processing" log and one "summary" log to make sure the pipeline ran as expected.
demultiplexed_fastq
Each input fastq.gz file is first validated by fqtools to make sure it is formatted correctly before beginning. 
The pipeline uses Cutadapt to demultiplex each input fastq.gz file into multiple demultiplexed fastq.gz files stored in this directory. Cutadapt pulls the TSO barcode from the 5' end of each read and matches it against the FASTA files that we generated above to demultiplex the reads. The TSO number and the TSO offset (0 or 3 bp) is stored in the filename of the output fastq.gz.
trimmed_fastq
Another round of Cutadapt is performed on the demultiplexed fastq.gz. This round does multiple things and outputs trimmed fastq.gz files into this sub-directory:
The offset is removed from the 5' end
The UMI is removed and appended to the read name
The constant TSO adapter is removed from the 5' end
Low quality bases, Illumina adapters, and poly-A sequences are removed from the 3' end
Any reads now shorter than a minimum length (default=20bp) are removed
The trimmed fastq.gz files are validated by fqtools to make sure no formatting corruption occured.
fastqc
FastQC is run on the trimmed fastq.gz files. This outputs a bunch of summary data about the reads in each file and is useful for QC.
aligned_bam
The reads are now aligned to the human genome using STAR. The output .bam files are stored in this sub-directory. STAR has many options to tune the alignment. I have not yet done an exhaustive test of these parameters. 
These aligned .bam files are then validated and indexed by samtools and a .bai index file is generated.
deduplicated_bam
The aligned reads are UMI deduplicated using UMI-tools to remove PCR duplicated reads. The deduplicated .bam files are stored here 
These deduplicated .bam files are then validated and indexed by samtools and a .bai index file is generated.
rseqc
RSEQC is a set of python scripts that generate nice summary data from the deduplicated .bam files and is useful for QC.
qualimap
Qualimap is another summary tool run on the deduplicated .bam files and is useful for QC.
featurecounts
This program is called in the smart3seq_summary.sh script.
featureCounts (part of Subread) is used to count the reads for each gene in each sample. The output .tsv file is the count table that we will use as the basis for our further analyses.
multiqc
This program is called in the smart3seq_summary.sh script.
MultiQC scans through all of log files generated by the tools above and summarizes the output into a useful .html report. We will check the report below. The .zip file contains the data for all of the graphs in the report if you need to investigate something deeper or want to make your own graphs.

