# SMART-3seq-HT
This repository contains scripts to demultiplex and analyze SMART-3seq-HT sequencing reads. If you are looking for a wet-bench protocol for generating SMART-3seq-HT libraries, consult this [Benchling notebook](https://benchling.com/s/prt-nFxom6I089vjB7FXlLiM?m=slm-qGI9cuwzvx6r4ctuojSr)

## Quick-Start Guide (TL;DR version)
1. clone this repo to your folder: `git clone https://github.com/nolanmaier/SMART-3seq-HT`
2. make a sample demultiplexing *.csv* file and upload it to your folder
3. run the wrapper script using slurm: `sbatch SMART-3seq-HT/smart3seq_wrapper.sh -i path/to/input_fastq_dir -o path/to/output_dir -t sample_demultiplexing.csv`

## Repository Contents
- `samplesheet_helper.py` python3 script that compiles multiple csv files into a samplesheet csv
- `smart3seq_wrapper.sh` bash script that submits the rest of pipeline. Most parameters that an end-user will need to tweak (such as sbatch variables and genome locations) are in this script.
- `smart3seq_processing.sh` bash script that demultiplexes, trims, aligns, deduplicates, and performs quality control on a input *fastq.gz* file using [Cutadapt](https://cutadapt.readthedocs.io/en/stable/index.html), [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [STAR](https://github.com/alexdobin/STAR), [SAMtools](http://www.htslib.org/doc/samtools.html), [UMI-tools](https://umi-tools.readthedocs.io/en/latest/), [RSEQC](https://rseqc.sourceforge.net/), and [Qualimap](http://qualimap.conesalab.org/).
- `smart3seq_summary.sh` bash script that summarizes deduplicated *.bam* files into a counts table using [featureCounts](https://subread.sourceforge.net/featureCounts.html) and generates a quality control report using [MultiQC](https://multiqc.info/docs/). 
- `barcode2fasta.py` python3 script that generates *.fasta* files containing TSO barcodes for each oligo-dT
- `sample_demultiplexing.csv` a minimal example sample demultiplexing *.csv* to be used as a reference when generating experiment specific sample demultiplexing *.csv* files.
- `SMART-3seq-HT_barcodes.xlsx` an excel file containing the sequences and barcodese of all of the oligo-dT and TSO. Used as a reference when generating experiment specific sample demultiplexing *.csv* files.

## Requirements
This repository is designed and tested for use on the Harvard Medical School [O2 computing cluster](https://harvardmed.atlassian.net/wiki/spaces/O2/overview). 

As written, it minimally requires a [SLURM](https://slurm.schedmd.com/documentation.html) job scheduler and a [BioGrids](https://biogrids.org/) software stack available as a module. You can check if these are available by running the following commands: `sbatch -V` and `module load biogrids/latest`

## Inputs
The wrapper script takes three named inputs:

`-i` the directory containing the input *.fastq.gz* files. The filenames must match the pattern: `*dT[0-9][0-9]*.fastq.gz` where the two-digit number following *dT* should match the **dT_Index** in the sample demultiplexing *.csv* (see below). The following are valid example filenames: *dT03.fastq.gz*, 	*LIB056_S1_dT12.fastq.gz*, *dT20TRA00243971.fastq.gz*

`-o` the directory where the output files should be generated

`-t` the sample demultiplexing *.csv* file (see below)

## Sample Demultiplexing CSV
The pipeline requires a comma separated file containing all samples and their associated oligo-dT and TSO barcodes. Reference the example in the repository. The `samplesheet_helper.py` script can be used to generate this csv file. Benchling plate layouts can be exported to csv files. These csv files can be used as input for the helper script to generate the samplesheet.

The minimum columns required are below (the column names must match exactly): 
- **dT_Index**
- **TSO_Index**
- **TSO_Offset**
- **TSO_Sequence**

I usually include optional extra columns (that are not parsed by the pipeline) describing the samples (i.e treatment names, treatment concentrations, etc) so that I can use this same *.csv* file during data analysis to decode the samples into treatments and controls. It is important to **not use commas** within fields in your *.csv* or the *.csv* will not be read correctly. If you must use commas (such as for sample names) enclose the entire field in double quotes e.g.:	✅ "3,4-hexanediol"	❌ 3,4-hexanediol.

## Genome/Transcriptome annotations
The pipeline also requires three different genome/transcriptome annotation files. The location of these files are hardcoded in the `smart3seq_wrapper.sh` script (lines 60-67). If you are not using the human genome or do not have access permissions for the Jost lab directories on O2, you will need to get your own versions of these files and change the appropriate lines in the wrapper script.

- `genomedir` the directory containing the STAR genome index for your organism of interest. This index should be generated beforehand using `STAR --runMode genomeGenerate --genomeDir /path/to/genomeDir --genomeFastaFiles /path/to/genome/fasta --sjdbGTFfile /path/to/annotations.gtf`. The input *.fasta* files should be **genomic** sequences (one chromosome per file) and the *.gtf* file should be the corresponding **transcript** annotations. [GENCODE](https://www.gencodegenes.org/) has a simple download page for human and mouse. See the [STAR](https://github.com/alexdobin/STAR) documentation for more information.
- `gtffile` the transcript annotation file for your organism of interest in *.gtf* format. The *.gtf* file should be the same one used above to generate the splice-junction database for the STAR genome generation. [GENCODE](https://www.gencodegenes.org/) has a simple download page for human and mouse.
- `bedfile` the transcript annotation file for your organims of interest in *.bed* format. This should be the same genome assembly and annotation as the *.gtf* file and can be a bit tricky to find. For human GENCODE annotations the [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1650209460_0G3ATy4TF3Sa3DTp0Cts3VMkDxe3&clade=mammal&org=Human&db=hg38&hgta_group=genes&hgta_track=fixSeqLiftOverPsl&hgta_table=0&hgta_regionType=genome&position=chr2%3A25%2C160%2C915-25%2C168%2C903&hgta_outputType=bed&hgta_outFileName=) provides a download portal. In theory, *.gtf* files can be [directly converted](https://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/gtf2bed.html) to *.bed* format but I have not explored this deeply.

## Detailed Instructions
1. Make a sample demultiplexing *.csv* as outlined above
2. Verify that the correct genome and transcriptome annotations are available and hardcoded correctly in the `smart3seq_wrapper.sh` script
3. Login to the HMS O2 HPC cluster. I like to use the [OpenOnDemand portal](https://o2portal.rc.hms.harvard.edu/). From there you can access a shell by clicking `Clusters` > `O2 Cluster Terminal`
4. Running analyses on the login nodes is frowned upon, so request an interactive session and wait a moment until your session starts
`srun --pty -p interactive -t 0-01:00 /bin/bash`
5. Copy the scripts and sample demultiplexing *.csv* to your directory on the O2 server. This can be done from your local machine using the `Upload` button on the OpenOnDemand web interface. Alternatively, use git to clone this repository from GitHub `git clone https://github.com/nolanmaier/SMART-3seq-HT`
6. Submit the `smart3seq_wrapper.sh` script as a SLURM `sbatch` job with the appropriate inputs. For example: `sbatch SMART-3seq-HT/smart3seq_wrapper.sh -i path/to/input_fastq_dir -o path/to/output_dir -t sample_demultiplexing.csv`

## Pipeline Monitoring
The HMS Research Computing Group has great [guides](https://harvardmed.atlassian.net/wiki/spaces/O2/pages/1586793632/Using+Slurm+Basic) for running and monitoring jobs on O2 using SLURM.You can [monitor](https://harvardmed.atlassian.net/wiki/spaces/O2/pages/1586793632/Using+Slurm+Basic#Monitoring-Jobs) the progress of the pipeline using the `O2squeue` command.

You should originally see a short 5 minute job queued. This is the `smart3seq_wrapper.sh` script which submits the other scripts. If you need to update the memory requirements or requested cpu time of the pipeline, you can edit the parameters in this script using the [sbatch options](https://harvardmed.atlassian.net/wiki/spaces/O2/pages/1586793632/Using+Slurm+Basic#sbatch-options-quick-reference).

Once the wrapper completes you will see one job in the queue for each *.fastq.gz* file in your input directory, These jobs are all running the `smart3seq_processing.sh` script and will run in parallel to each other as a [SLURM array](https://harvardmed.atlassian.net/wiki/spaces/O2/pages/1586793632/Using+Slurm+Basic#Job-Arrays). There will be one additional job in the queue marked as `STATE: DEPENDENCY`. This is the `smart3seq_summary.sh` script which relies on the the previous jobs using [SLURM dependencies](https://harvardmed.atlassian.net/wiki/spaces/O2/pages/1586793632/Using+Slurm+Basic#Job-Dependencies). This job will not run until all of the `smart3seq_processing.sh` jobs finish. Once the `smart3seq_summary.sh` script finishes, the pipeline is complete. You should receive an email notification upon successful pipeline completion or if any of the jobs fail.

## Output Files
The pipeline will create a bunch of sub-directories within the output directory you indicated above. These sub-directories correspond to the various steps in the pipeline.

### slurmlogs
This directory contains the log files from all of the SLURM jobs. **If the pipeline fails, this is the first place to look for explanations.** Even if the pipeline completes, it is useful to scan through at least one "processing" log and one "summary" log to make sure the pipeline ran as expected and there are no obvious errors.

### TSO_barcodes
Demultiplexing of the oligo-dT barcodes should have previously been done (usually by the sequencer itself) yielding one *.fastq.gz* file per oligo-dT. The read processing pipeline will demultiplex each of these files further according to the TSO barcodes yielding one *.fastq.gz* file per oligo-dT/TSO pair. This step requires one *.fasta* file per oligo-dT, containing the sequences of the TSO barcodes that correspond to that oligo-dT. The `barcode2fasta.py` script is a simple python script that uses the sample demultiplexing *.csv* to generate these *.fasta* files and places them in the TSO_barcodes directory. The *.fasta* files also encode information about the TSO offset which is used in the pipeline for read trimming.

### demultiplexed_fastq
The pipeline uses [Cutadapt](https://cutadapt.readthedocs.io/en/stable/index.html) to validate then demultiplex each input *.fastq.gz* file into multiple demultiplexed *.fastq.gz* files stored in this directory. Cutadapt pulls the TSO barcode from the 5' end of each read and matches it against the *.fasta* files in the TSO_barcodes directory. The TSO number and the TSO offset (0 or 3 bp) is stored in the filename of the output *.fastq.gz*

### trimmed_fastq
[Cutadapt](https://cutadapt.readthedocs.io/en/stable/index.html) is used again on the demultiplexed *.fastq.gz* and performs many operations in order:
1. The 3bp offset is removed from the 5' end of the read (if offset=3) or from the 3' end of the read (if offset=0; this maintains an equal length for all reads).
2. The 5bp UMI is removed from the 5' end and appended to the read name
3. The constant TSO adapter sequence is removed from the 5' end
4. Low quality bases are removed from the 3' end
5. Illumina adapter sequences are removed from the 3' end
6. poly-A sequences are removed from the 3' end
7. Trailing N bases are removed from the 5' and 3' ends
8. Any reads now shorter than a 10bp are removed entirely
9. The remaining reads are written to *.fastq.gz* files in the trimmed_fastq directory
10. The trimmed *.fastq.gz* files are validated to make sure everything is formatted correctly

### fastqc
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is run on the trimmed *.fastq.gz* files and a bunch of summary data is output into this directory. This data is useful for quality control.

### aligned_bam
The reads are now aligned to the human genome using [STAR](https://github.com/alexdobin/STAR). The output *.bam* files are stored in this directory. These aligned *.bam* files are then validated and indexed by [SAMtools](http://www.htslib.org/doc/samtools.html) and a *.bai* index file is generated.

### deduplicated_bam
The aligned reads are UMI deduplicated using [UMI-tools](https://umi-tools.readthedocs.io/en/latest/) to remove PCR duplicated reads. The deduplicated *.bam* files are stored in this directory. These deduplicated *.bam* files are then validated and indexed by [SAMtools](http://www.htslib.org/doc/samtools.html) and a *.bai* index file is generated.

### rseqc
[RSEQC](https://rseqc.sourceforge.net/) is a set of python scripts that analyzes the deduplicated *.bam* files. The summary data stored in this directory is useful for quality control.

### qualimap
[Qualimap](http://qualimap.conesalab.org/) is a Java tool that analyzes the deduplicated *.bam* files. The summary data stored in this directory is useful for quality control.

### featurecounts
[featureCounts](https://subread.sourceforge.net/featureCounts.html) (part of [Subread](https://subread.sourceforge.net/)) is used to count the reads for each gene in each sample. The output *.tsv* file in this directory is **the count table that you will use as the basis for further analyses** (using for example DESeq2 or edgeR). This program is called in the `smart3seq_summary.sh` script (once per pipeline).

### multiqc
[MultiQC](https://multiqc.info/docs/) scans through all of log files generated by the tools above and summarizes the output into a useful *.html* report stored in this directoy. **You should definitely open this report and scan through it to make sure that the pipeline ran as expected and identify any outlier samples.** The *.zip* file in this directory contains the data for all of the graphs in the report if you need to investigate something deeper or want to make your own graphs. This program is called in the `smart3seq_summary.sh` script (once per pipeline).

