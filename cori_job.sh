#!/bin/bash
#BATCH -N 4
#SBATCH -C haswell
#SBATCH -q debug
#SBATCH -J caldisks
#SBATCH -o caldisks.%j.stdout
#SBATCH -e caldisks.%j.error
#SBATCH --mail-user=erk24@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH -t 30:00

module swap PrgEnv-intel PrgEnv-gnu
module load openmpi
rm caldiskstest
cd include/sketch/BigInt/
make
cd ../../../
./run_tests.sh 
mpirun -np 64 caldiskstest celegans500.fastq

