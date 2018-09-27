#!/bin/bash

#SBATCH --job-name=combine_vol
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=nesi00263
#SBATCH --partition=NeSI
#SBATCH --hint=nomultithread
#SBATCH --time 0:15:00

currentdir=`pwd`

echo
echo "  ./bin/xcombine_vol_data_vtk 0 3 beta_kernel_smooth OUTPUT_FILES/DATABASES_MPI/ OUTPUT_FILES/ 0"
echo

##srun -n 1 ./bin/xcombine_vol_data_vtk 0 3 alpha_kernel_smooth OUTPUT_FILES/DATABASES_MPI/ OUTPUT_FILES/ 0
srun -n 1 ./bin/xcombine_vol_data_vtk 0 3 beta_kernel_smooth OUTPUT_FILES/DATABASES_MPI/ OUTPUT_FILES/ 0

echo
echo "done"


