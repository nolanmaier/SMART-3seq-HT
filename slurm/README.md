README.md

# login to slurm

# initialize interactive job
srun --pty -p interactive -t 0-04:00 /bin/bash

# move to the correct directory
cd /path/to/directory

# copy scripts folder and data folder to directory

# run script
sbatch ./slurm/smart3seq_wrapper.sh -i ./test_data -o ./out -t ./test_data/TSO_barcodes.fasta