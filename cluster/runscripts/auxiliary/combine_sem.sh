#!/bin/bash

#SBATCH --job-name=combine_sem
#SBATCH --nodes=4
#SBATCH --ntasks=144
#SBATCH --cpus-per-task=1
#SBATCH --account=nesi00263
#SBATCH --partition=nesi_research
#SBATCH --hint=nomultithread
#SBATCH --time 0:15:00

currentdir=`pwd`

echo
echo "srun -n 144 ./bin/xcombine_sem alpha_kernel,beta_kernel,rho_kernel kernels_list.txt OUTPUT_SUM/"
echo

srun -n 144 ./bin/xcombine_sem alpha_kernel,beta_kernel,rho_kernel kernels_list.txt OUTPUT_SUM_TEST/ 

echo
echo "done"


