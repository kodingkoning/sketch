#!/bin/bash
rm caldiskstest
cd include/sketch/BigInt/
make
cd ../../../
./run_tests.sh && mpirun -np 2 caldiskstest ecsample1.fastq

