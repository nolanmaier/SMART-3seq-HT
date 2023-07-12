#!/usr/bin/env python3

# import libraries
import argparse, csv
from pathlib import Path

# write a detailed program description
description_message = '''
This script will generate a samplesheet csv file for the SMART-3Seq-HT processing pipeline

All arguments should be filenames or paths to csv files containing a 96- or 384-well plate layout
Such csv files can be downloaded from Benchling with a single click
The dT-Index csv and TSO-Index csv are required
Any additional arguments are optional columns that will be added to the samplesheet
'''

epilog_message='''
optional positional arguments:
    [column_name]=[filename]    [column name] is the requested name of the column in the final sample sheet 
                                [csv filename] is the filename or path to the csv file containing that information
                                e.g. "trx_name=treatments.csv"
                                multi-line entries will be split into multiple columns
'''
script_name = "samplesheet_helper.py"

# assert the correct number and format of command line arguments
parser = argparse.ArgumentParser(prog=script_name, 
                                 formatter_class=argparse.RawDescriptionHelpFormatter, 
                                 description=description_message,
                                 epilog=epilog_message)
parser.add_argument("dT_Index", type=Path, help="path to a csv file containing oligo-dT indices in a plate layout")
parser.add_argument("TSO_Index", type=Path, help="path to a csv file containing TSO indices in a plate layout")
parser.add_argument("-o", "--outfile", type=Path, default="samplesheet.csv", help="(optional) name of the output csv file (default: 'samplesheet.csv')")
parser.add_argument("-t", "--tso_sequences", default=None, help="(optional) name of the csv file containing custom TSO sequences and offsets (default: standard SMART-3Seq-HT sequences)")

# parse the command line arguments
args, unknownargs = parser.parse_known_args()
print(f"starting {script_name}")

# hardcode the default TSO barcodes and offset into a dictionary
default_TSO_seq = {
    1: ["GGTAGCTA", 0],
    2: ["GTATCAGC",	3],
    3: ["CCATTAAT", 0],
    4: ["TGGGACAT",	3],
    5: ["CGTTCGCA",	0],
    6: ["TCACCTTG",	3],
    7: ["ATTGTCTA",	0],
    8: ["CGACTCAG",	3],
    9: ["AATGGTCG",	0],
    10: ["CTGTTCGT", 3],
    11: ["CTACTATA", 0],
    12: ["CGTTTAGT", 3],
    13: ["TCCGATGT", 0],
    14: ["GATTCCAC", 3],
    15: ["TGGCGACC", 0],
    16: ["CTTTTGGC", 3],
    17: ["ACGACAAT", 0],
    18: ["GAAATTAC", 3],
    19: ["CGACTACC", 0],
    20: ["GTTCAGCT", 3],
    21: ["TTATGGTT", 0],
    22: ["GCAAGTCG", 3],
    23: ["CAGCCTGC", 0],
    24: ["TGAATCGT", 3],
    25: ["CTGTCATA", 0],
    26: ["GAATAGGC", 3],
    27: ["CGATGGTA", 0],
    28: ["CATCCCTC", 3],
    29: ["GCTAACCA", 0],
    30: ["TTTAAGTG", 3],
    31: ["ATCGGAAT", 0],
    32: ["TAACAGTC", 3]
    }

def extract_custom_TSO_sequences(csv_in):
    ''' function to extract custom TSO sequences and offsets
        accepts a Path to a csv file with three columns
        header must be: ["TSO_Index", "TSO_Sequence", "TSO_Offset"]
        values should be: [int, str, int]
        generates a dictionary with TSO_Index as keys
        and TSO_Sequence, TSO_Offset in a list as the values'''
    TSO_dict = {}
    assert csv_in.is_file(), f"arguments must be paths to valid files, '{csv_in}' does not exist"
    with open(csv_in, "r") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            TSO_dict[int(row["TSO_Index"])] = [str(row["TSO_Sequence"]), int(row["TSO_Offset"])]
    return TSO_dict

def csv_to_wells(csv_in):
    ''' function to read a plate layout csv
        generates a dictionary with well labels as keys (e.g. 'A1') 
        and cell entries as values '''
    well_count = 0
    well_dict = {}
    # open and read the csv
    assert csv_in.is_file(), f"arguments must be paths to valid files, '{csv_in}' does not exist"
    with open(csv_in, "r") as csvfile:
        reader = csv.DictReader(csvfile)
        # iterate through each row of the csv
        for row in reader:
            # extract the row label (e.g. 'A')
            prefix = row['']
            # iterate over each well in the row
            for column in row:
                # skip any column that does not have a column label (the first column of the csv)
                if column:
                    well_count += 1
                    # skip any well that is empty
                    if row[column]:
                        well_dict[prefix+column] = row[column]
    return well_count, well_dict

# convert the command line arguments into file paths
csv_paths = {}
csv_paths["dT_Index"] = args.dT_Index
csv_paths["TSO_Index"] = args.TSO_Index

# check that the extra arguments don't use a restricted label
for argument in unknownargs:
    assert argument.partition("=")[0] not in args, \
        f"received invalid optional argument '{argument}', optional arguments cannot have a restricted name in {list(vars(args).keys())}'"
    # split the extra arguments into a label and a path
    csv_paths[argument.partition("=")[0]] = Path(argument.partition("=")[2])

# initialize a dict to store the number of wells in each csv
counts = {}
# initialize a dict to store the well data as a dictionary of dictionaries
data = {}
# iterate through each of the csv files
for key, filepath in csv_paths.items():
    # read each csv and append the appropriate data to the storage dict
    print(f"reading '{key}' from '{filepath}'")
    counts[key], data[key] = csv_to_wells(filepath)

# check the plate format of each csv matches an expected layout and each other
for key, value in counts.items():
    assert value in [96, 384], f"96- or 384-well plate format expected, read {value} wells in '{key}'"
    assert value == counts['dT_Index'], \
        f"all files should contain the same plate format, 'dT_Index' has {counts['dT_Index']}-well format but '{key}' has {value}-well format"
print(f"{counts['dT_Index']}-well plate format detected")

# check that each csv has the same number of samples and the same wells populated
for key, value in data.items():
    assert len(data[key]) == len(data['dT_Index']), \
        f"all files should contain the same number of samples, 'dT_Index' has {len(data['dT_Index'])} samples but '{key}' has {len(data[key])} samples"
    for well in value.keys():
        assert well in data['dT_Index'].keys(), \
            f"all files should contain sample information for the same wells, well {well} is in '{key}' but not in 'dT_Index"
print(f"{len(data['dT_Index'])} samples detected")

# if the user want custom TSO sequences, extract the sequences from the csv
if args.tso_sequences:
    print(f"using custom TSO sequences and offsets provided by the user in {args.tso_sequences}")
    TSO_decoder = extract_custom_TSO_seq(Path(args.tso_sequences))
# if the user want to use the default TSO sequences, use the hardcoded dictionary
else:
    print("using default TSO sequences and offsets")
    TSO_decoder = default_TSO_seq

# initialize a dictionary to store the cleaned, concatenated data in a format for export
outdict = {}

for well, value in data['dT_Index'].items():
    # add the well number to the output dictionary
    outdict[well] = {}
    outdict[well]['well'] = well
    # add the dT_Index, making sure to clean up the entry to only the numbers
    clean_index = int(''.join(c for c in value if c.isdigit()))
    outdict[well]['dT_Index'] = clean_index

for well, value in data['TSO_Index'].items():
    # add the TSO_Index, making sure to clean up the entry to only the numbers and that the index is in the barcode list
    clean_index = int(''.join(c for c in value if c.isdigit()))
    assert clean_index in TSO_decoder.keys(), \
        f"TSO_Index must be in hardcoded barcode list, could not find {value} or {clean_index} from well {well}"
    outdict[well]['TSO_Index'] = clean_index
    # also add the appropriate sequence and offset from the decoder dictionary
    outdict[well]['TSO_Sequence'] = TSO_decoder[clean_index][0]
    outdict[well]['TSO_Offset'] = TSO_decoder[clean_index][1]

# add any additional information from the other csv files
for key in set(data.keys()) - {'dT_Index', 'TSO_Index'}:
    # csv entries may be multi-line, check what the maximum number of lines is
    max_lines = 0
    for value in data[key].values():
        max_lines = max(max_lines, len(value.splitlines()))
    if max_lines > 1:
        print(f"detected {max_lines} lines per entry in '{key}', splitting into multiple columns: {[key+'_'+str(i+1) for i in range(max_lines)]}")
        # initialize an empty list of the max length for each well
        for well, value in data[key].items():
            split_data = [""] * max_lines
            # split the lines into the list
            for i, x in enumerate(value.splitlines()):
                split_data[i] = x
            # write the lines to the output dictionary, appending an integer to the name if needed
            for i in range(max_lines):
                outdict[well][key+'_'+str(i+1)] = split_data[i]
    else:
        for well, value in data[key].items():
            outdict[well][key] = value

# write the output csv
print(f"writing samples to '{args.outfile}'")
with open(args.outfile, 'w') as csvfile:
    fieldnames = list(outdict[next(iter(outdict))].keys())
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for entry in outdict.values():
        writer.writerow(entry)
print("DONE")