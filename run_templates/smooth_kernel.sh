#!/bin/bash

#SBATCH --job-name=smooth_kernel
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --account=nesi00263
#SBATCH --partition=NeSI
#SBATCH --hint=nomultithread
#SBATCH --time 0:15:00

currentdir=`pwd`

echo
echo "  ./bin/xsmooth_sem 5000. 1000. beta_kernel OUTPUT_FILES/DATABASES_MPI/ OUTPUT_FILES/DATABASES_MPI/ .false"
echo

srun -n 4 ./bin/xsmooth_sem 5000. 1000. beta_kernel OUTPUT_FILES/DATABASES_MPI/ OUTPUT_FILES/DATABASES_MPI/ .false

echo
echo "done"


