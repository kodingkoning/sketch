#!/bin/bash
module swap PrgEnv-intel PrgEnv-gnu
module load openmpi
rm caldiskstest
cd include/sketch/BigInt/
make
cd ../../../
./run_tests.sh && mpirun -np 2 caldiskstest ecsample1.fastq

