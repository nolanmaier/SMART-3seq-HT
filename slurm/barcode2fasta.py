#!/usr/bin/env python3

''' 
To demultiplex SMART-3Seqv2 data, we need a fasta file containing relevant TSO barcodes for each dT barcode
This python script reads a csv file containing sample information and outputs the fasta files
The csv file must have the following columns (case-sensitive):
	dT_Index
	TSO_Index
	TSO_Offset
	TSO_Sequence


Usage:
python3 barcode2fasta.py sample_info.csv output_directory
'''

# import modules
import os, platform, sys, csv, shutil, glob
from pathlib import Path

print(f"running {__file__}")
print(f"using python version: {platform.python_version()}")

# read sample list from path given at command line
sample_csv = Path.cwd() / sys.argv[1]
print(f"reading sample information from {sample_csv}")

# make a folder to store the output fasta files
outdir = Path.cwd() / sys.argv[2]
# remove the directory and its contents if it already exists to avoid appending to an existing file
try:
    shutil.rmtree(outdir)
except FileNotFoundError:
    pass
# make the empty output directory
os.makedirs(outdir)
print(f"exporting .fasta files to {outdir}")

# initialize empty lists to hold the data
dT_Index = []
TSO_Index = []
TSO_Offset = []
TSO_Sequence = []
# read the csv file row by row
with open(sample_csv, mode='r', encoding='utf-8-sig') as csvfile:
    reader = csv.DictReader(csvfile)
    sample_count = 0
    # extract the relevant column information to the appropriate lists
    for row in reader:
        dT_Index.append(int(row["dT_Index"]))
        TSO_Index.append(int(row["TSO_Index"]))
        TSO_Offset.append(int(row["TSO_Offset"]))
        TSO_Sequence.append(str(row["TSO_Sequence"]))
        sample_count += 1
# zip the lists together so that each tuple is one entry
data = list(zip(dT_Index, TSO_Index, TSO_Offset, TSO_Sequence))
print(f"{sample_count} samples found")

# append each entry to the appropriate fasta file
for entry in data:
    outname = outdir / f"dT{format(entry[0], '02')}_barcodes.fasta"
    with open(outname, mode='a') as outfile:
        outfile.write(f">TSO{format(entry[1], '02')}_OFF{entry[2]}\n{entry[3]}\n")

# count the barcodes in each fasta file
output_file_list = sorted(glob.glob(str(outdir / 'dT*_barcodes.fasta')))
print(f"{len(output_file_list)} fasta files written")
for f in output_file_list:
    f_path = Path(f)
    with open(f_path, mode='r') as fasta_file:
        barcode_count, remainder = divmod(len(fasta_file.readlines()), 2)
        if remainder != 0:
            raise AssertionError(f"wrong number of lines detected in {f_path.name}")
        print(f"{barcode_count} barcodes in {f_path.name}")

print("DONE")