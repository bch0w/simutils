#!/bin/bash

#SBATCH --job-name=combine_vol_data_vtk
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=nesi00263
#SBATCH --partition=nesi_research
#SBATCH --hint=nomultithread
#SBATCH --time 0:01:30

currentdir=`pwd`

srun -n 1 ./bin/xcombine_vol_data_vtk 0 143 beta_kernel ./OUTPUT_SUM/ ./OUTPUT_SUM/ 0

echo
echo "done"


