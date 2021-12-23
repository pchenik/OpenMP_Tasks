#!/bin/sh
#module add mpi/openmpi-local
module add hpc-pract-2021 
#module add mhpc
hostname
export OMP_NUM_THREADS=16
./main

