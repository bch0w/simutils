#!/bin/bash -e

#SBATCH --job-name=xsmooth_sem
#SBATCH --nodes=4
#SBATCH --ntasks=144
#SBATCH --cpus-per-task=1
#SBATCH --account=nesi00263
#SBATCH --clusters=maui
#SBATCH --partition=nesi_research
#SBATCH --time=00:30:00
#SBATCH --output=smooth_sem_%j.out

# Set the compiler option
COMPILER=SPECFEM3D/20190730-CrayCCE-19.04
COMPILER=SPECFEM3D/20190730-CrayGNU-19.04
COMPILER=SPECFEM3D/20190730-CrayIntel-19.04

module load ${COMPILER}

# Kernel to smooth and smoothing parameters must be specified by user
KERNEL="vs"
SGMAH=40000.
SGMAV=1000.
DIR_IN="SMOOTH/"
DIR_OUT=${DIR_IN}
USE_GPU=".false"

# Get the number of processors from Par_file, ignore comments
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

echo ${COMPILER}
echo "xsmooth_sem ${KERNEL} w/ sigma_h=${SGMAH}, sigma_v=${SGMAV}"
echo "${NPROC} processors, GPU option: ${USE_GPU}"
echo
echo "`date`"
# EXAMPLE CALL:
# srun -n NPROC xmooth_sem SIGMA_H SIGMA_V KERNEL_NAME INPUT_DIR OUTPUT_DIR USE_GPU
time srun -n ${NPROC} xsmooth_sem ${SGMAH} ${SGMAV} ${KERNEL} ${DIR_IN} ${DIR_OUT} ${USE_GPU}

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "finished at: `date`"
