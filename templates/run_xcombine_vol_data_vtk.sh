#!/bin/bash

#SBATCH --job-name=combine_vol_data_vtk
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=nesi00263
#SBATCH --clusters=maui
#SBATCH --partition=nesi_research
#SBATCH --hint=nomultithread
#SBATCH --time=0:00:30
#SBATCH --output=sum.log

# EXAMPLE CALL sbatch combine_vol_data_vtk 2018p130600 beta_kernel_smooth
# TO DO add NPROC figure outer here

echo "`date`"

QUANTITY=$1
if [ -z "$1" ]
then
	echo "QUANTITY NEEDS TO BE SPECIFIED e.g. vs, hess_kernel, beta_kernel_smooth"
	exit
fi

NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
NPROC_END=`expr $NPROC - 1`

currentdir=`pwd`
DIR_IN="./SUM/"
DIR_OUT=${DIR_IN}

#srun -n nproc ./bin/xcombine_vol_data_vtk proc_start proc_end kernel dir_in dir_out hi_res
srun -n 1 ./bin/xcombine_vol_data_vtk 0 $NPROC_END ${QUANTITY} ${DIR_IN}/ ${DIR_OUT}/ 0

echo
echo "done"


