#!/bin/bash

#SBATCH --job-name=combine_vol_data_vtk
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=nesi00263
#SBATCH --clusters=maui
#SBATCH --partition=nesi_research
#SBATCH --time=0:00:30
#SBATCH --output=combine_vol_data_vtk_%j.out

# Set the compiler option
# COMPILER=SPECFEM3D/20190730-CrayCCE-19.04
# COMPILER=SPECFEM3D/20190730-CrayIntel-19.04
module load gcc/8.3.0
COMPILER=SPECFEM3D/20190730-CrayGNU-19.04

module load ${COMPILER}

# Quantity needs to be specified by the user
QUANTITY=$1
if [ -z "$1" ]
then
	echo "QUANTITY REQUIRED (e.g. vs, hess_kernel, beta_kernel_smooth)"
	exit
fi

# Dynamically get the number of processors from the Par_file
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
NPROC_START=0
NPROC_END=`expr $NPROC - 1`

# Set the paths for Specfem to search
DIR_IN="./SUM/"
DIR_OUT=${DIR_IN}

# Example Call
# srun -n nproc xcombine_vol_data_vtk proc_start proc_end kernel dir_in dir_out hi_res

# Run the Exectuable
echo ${COMPILER}
echo "xcombine_vol_data_vtk ${NPROC_START} ${NPROC_END} for ${QUANTITY}"
echo
echo "`date`"
time srun -n 1 xcombine_vol_data_vtk ${NPROC_START} ${NPROC_END} ${QUANTITY} ${DIR_IN}/ ${DIR_OUT}/ 0
echo
echo "finished at: `date`"


