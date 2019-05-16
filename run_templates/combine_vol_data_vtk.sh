#!/bin/bash

#SBATCH --job-name=combine_vol_data_vtk
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=nesi00263
#SBATCH --clusters=maui
#SBATCH --partition=nesi_research
#SBATCH --hint=nomultithread
#SBATCH --time 0:00:30

# EXAMPLE CALL sbatch combine_vol_data_vtk 2018p130600 beta_kernel_smooth
# TO DO add NPROC figure outer here

KERNEL=$1
if [ -z "$1" ]
then
	echo "KERNEL NEEDS TO BE SPECIFIED e.g. hess_kernel, beta_kernel_smooth"
	exit
fi

echo "`date`"
currentdir=`pwd`
DIR_IN="./COMBINE_MODEL/"
# DIR_IN="./OUTPUT_SUM/"
DIR_OUT=${DIR_IN}

#srun -n nproc ./bin/xcombine_vol_data_vtk proc_start proc_end kernel dir_in dir_out gpu_accel
srun -n 1 ./bin/xcombine_vol_data_vtk 0 143 ${KERNEL} ${DIR_IN} ${DIR_OUT} 0

echo
echo "done"


