#!/bin/bash -e

#SBATCH --job-name=xsmooth_sem
#SBATCH --nodes=2
#SBATCH --ntasks=80
#SBATCH --cpus-per-task=1
#SBATCH --account=nesi00263
#SBATCH --clusters=maui
#SBATCH --partition=nesi_research
#SBATCH --hint=nomultithread
#SBATCH --time=01:30:00
#SBATCH --output=smooth_sem_%j.out

# Kernel to smooth and smoothing parameters must be specified by user
KERNEL="beta_kernel"
SGMAH=2000.
SGMAV=1000.
DIR_IN="SMOOTH/"
DIR_OUT=${DIR_IN}
USE_GPU=".false"

# Get the number of processors from Par_file, ignore comments
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
NPROC_START=0
NPROC_END=`expr $NPROC - 1`

echo "xsmooth_sem ${KERNEL} w/ sigma_h=${SGMAH}, sigma_v=${SGMAV}"
echo "${NPROC} processors, GPU option: ${USE_GPU}"
echo
echo "`date`"
# EXAMPLE CALL:
# srun -n NPROC xmooth_sem SIGMA_H SIGMA_V KERNEL_NAME INPUT_DIR OUTPUT_DIR USE_GPU
time srun -n ${NPROC} ./bin/xsmooth_sem ${SGMAH} ${SGMAV} ${KERNEL} ${DIR_IN} ${DIR_OUT} ${USE_GPU}

# combine into vtk file 
time srun -n 1 ./bin/xcombine_vol_data_vtk \
        ${NPROC_START} ${NPROC_END} ${KERNEL}_smooth ${DIR_IN}/ ${DIR_OUT}/ 0


# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "finished at: `date`"

