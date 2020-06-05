#!/bin/bash
#BATCH -N 1
#SBATCH -C haswell
#SBATCH -q debug
#SBATCH -J caldisks
#SBATCH -o final/caldisks.%j.stdout
#SBATCH -e final/caldisks.%j.error
#SBATCH -t 30:00

module swap PrgEnv-intel PrgEnv-gnu
module load openmpi
rm caldiskstest
cd include/sketch/BigInt/
make
cd ../../../
./run_tests.sh 

# mpirun -np  1 caldiskstest 13 150 /project/projectdirs/m1982/bella-data/real-data/celegans40x_allfastq.fastq /project/projectdirs/m1982/bella-data/real-data/ecfull100x.fastq
# mpirun -np  2 caldiskstest 13 150 /project/projectdirs/m1982/bella-data/real-data/celegans40x_allfastq.fastq /project/projectdirs/m1982/bella-data/real-data/ecfull100x.fastq
# mpirun -np  4 caldiskstest 13 150 /project/projectdirs/m1982/bella-data/real-data/celegans40x_allfastq.fastq /project/projectdirs/m1982/bella-data/real-data/ecfull100x.fastq
mpirun -np  8 caldiskstest 13 150 /project/projectdirs/m1982/bella-data/real-data/celegans40x_allfastq.fastq /project/projectdirs/m1982/bella-data/real-data/ecfull100x.fastq
mpirun -np 16 caldiskstest 13 150 /project/projectdirs/m1982/bella-data/real-data/celegans40x_allfastq.fastq /project/projectdirs/m1982/bella-data/real-data/ecfull100x.fastq
mpirun -np 32 caldiskstest 13 150 /project/projectdirs/m1982/bella-data/real-data/celegans40x_allfastq.fastq /project/projectdirs/m1982/bella-data/real-data/ecfull100x.fastq

# mpirun -np 1 caldiskstest 7 50 test_a.txt test_b.txt
# mpirun -np 2 caldiskstest 7 50 test_a.txt test_b.txt
# mpirun -np 4 caldiskstest 7 50 test_a.txt test_b.txt
# mpirun -np 8 caldiskstest 7 50 test_a.txt test_b.txt
# mpirun -np 16 caldiskstest 7 50 test_a.txt test_b.txt


# mpirun -np 1 caldiskstest 21 150 celegans_e.fastq ecsample4000.fastq
# mpirun -np 2 caldiskstest 21 150 celegans_e.fastq ecsample4000.fastq
# mpirun -np 4 caldiskstest 21 150 celegans_e.fastq ecsample4000.fastq


#mpirun -np 1 caldiskstest 21 150 /project/projectdirs/m1982/bella-data/real-data/celegans40x_allfastq.fastq /project/projectdirs/m1982/bella-data/real-data/ecfull100x.fastq
# mpirun -np 4 caldiskstest 21 150 /project/projectdirs/m1982/bella-data/real-data/celegans40x_allfastq.fastq /project/projectdirs/m1982/bella-data/real-data/ecfull100x.fastq
#mpirun -np 8 caldiskstest 21 150 /project/projectdirs/m1982/bella-data/real-data/celegans40x_allfastq.fastq /project/projectdirs/m1982/bella-data/real-data/ecfull100x.fastq
#mpirun -np 16 caldiskstest 21 150 /project/projectdirs/m1982/bella-data/real-data/celegans40x_allfastq.fastq /project/projectdirs/m1982/bella-data/real-data/ecfull100x.fastq
#mpirun -np 4 caldiskstest 21 150 /project/projectdirs/m1982/bella-data/real-data/ecfull100x.fastq /project/projectdirs/m1982/bella-data/real-data/celegans40x_allfastq.fastq
#mpirun -np 5 caldiskstest 21 150 /project/projectdirs/m1982/bella-data/real-data/ecfull100x.fastq /project/projectdirs/m1982/bella-data/real-data/celegans40x_allfastq.fastq
#mpirun -np 1 caldiskstest 21 150 /project/projectdirs/m1982/bella-data/real-data/celegans40x_allfastq.fastq /project/projectdirs/m1982/bella-data/real-data/ecfull100x.fastq
