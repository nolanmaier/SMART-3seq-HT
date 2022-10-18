README.md

# login to slurm

# initialize interactive job
srun --pty -p interactive -t 0-04:00 /bin/bash

# move to the correct directory
test_dir='/n/data1/hms/microbiology/jost/lab/nolan/test/smartslurm_test/SMART-3seq-HT'
cd ${test_dir}

# copy scripts folder and data folder to directory

# run script
sbatch ${test_dir}/slurm/smart3seq_wrapper.sh -i ${test_dir}/test_data -o ${test_dir}/out -t ${test_dir}/test_data/TSO_barcodes.fasta 2>&1 | tee output`date +\%Y_\%m_\%d_\%H:\%M`.log