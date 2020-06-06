#!/bin/bash
#SBATCH -N 8
#SBATCH --tasks-per-node=32
#SBATCH -C haswell
#SBATCH -q debug
#SBATCH -J caldisks
#SBATCH -o final/caldisks.%j.stdout
#SBATCH -e final/caldisks.%j.error
#SBATCH -t 7:00
#DW jobdw capacity=10GB access_mode=striped type=scratch
#DW stage_in source=/global/cscratch1/sd/erk24/data/real-data/celegans40x_allfastq.fastq destination=$DW_JOB_STRIPED/celegans40x_allfastq.fastq type=file
#DW stage_in source=/global/cscratch1/sd/erk24/data/real-data/ecfull100x.fastq destination=$DW_JOB_STRIPED/ecfull100x.fastq type=file

# NOTE: the files must be located in the scratch buffer to work with Cori's burst buffer. This enables parallel I/O
# More info about the Burst buffer and scratch buffer here: https://www.nersc.gov/assets/Uploads/Burst-Buffer-tutorial.pdf
# The capacity requested should coordinate with the size of the files

echo "nodes =" $SLURM_JOB_NUM_NODES
echo "tasks/node =" $SLURM_NTASKS_PER_NODE
echo "Using k = 13 and sketch size = 150"
echo

TASKS=$(($SLURM_JOB_NUM_NODES*$SLURM_NTASKS_PER_NODE))

FILE_A=${DW_JOB_STRIPED}celegans40x_allfastq.fastq
FILE_B=${DW_JOB_STRIPED}ecfull100x.fastq

module swap PrgEnv-intel PrgEnv-gnu
module load openmpi
# rm caldiskstest # caldiskstest does not rebuild properly if not deleted. This should be uncommented during development
cd include/sketch/BigInt/
make
cd ../../../
./run_tests.sh 

mpirun -np $TASKS caldiskstest 13 150 $FILE_A $FILE_B
mpirun -np $TASKS caldiskstest 13 150 $FILE_A $FILE_B
mpirun -np $TASKS caldiskstest 13 150 $FILE_A $FILE_B
