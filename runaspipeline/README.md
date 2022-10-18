README.md

# login to slurm

# initialize interactive job
srun --pty -p interactive -t 0-04:00 /bin/bash

# move to the correct directory
test_dir='/n/data1/hms/microbiology/jost/lab/nolan/test/smartslurm_test/SMART-3seq-HT'
cd ${test_dir}

# copy scripts folder and data folder to directory

# load the rcbio module
module load rcbio/1.3.3

# test pipeline builder
runAsPipeline "${test_dir}/runaspipeline/smart3seq_wrapper.sh -i ${test_dir}/test_data -o ${test_dir}/out2 -t${test_dir}/test_data/TSO_barcodes.fasta" "sbatch -p short -t 0-00:10 -c 1" noTmp

# run pipeline builder
runAsPipeline "${test_dir}/runaspipeline/smart3seq_wrapper.sh -i ${test_dir}/test_data -o ${test_dir}/out2 -t ${test_dir}/test_data/TSO_barcodes.fasta" "sbatch -p short -t 0-00:10 -c 1" noTmp run 2>&1 | tee output`date +\%Y_\%m_\%d_\%H:\%M`.log