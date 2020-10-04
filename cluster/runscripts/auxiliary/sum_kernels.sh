#!/bin/bash

#SBATCH --job-name=sum_kernels
#SBATCH --nodes=4
#SBATCH --ntasks=144
#SBATCH --cpus-per-task=1
#SBATCH --account=nesi00263
#SBATCH --partition=nesi_research
#SBATCH --hint=nomultithread
#SBATCH --time 0:05:00

currentdir=`pwd`

echo
echo "  ./bin/xsum_kernels"
echo

srun -n 144 ./bin/xsum_kernels

echo
echo "done"


