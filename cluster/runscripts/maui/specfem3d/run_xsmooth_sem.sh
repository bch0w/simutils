#!/bin/bash -e

#SBATCH --job-name=xsmooth_sem
#SBATCH --nodes=2
#SBATCH --ntasks=80
#SBATCH --cpus-per-task=1
#SBATCH --account=nesi00263
#SBATCH --clusters=maui
#SBATCH --partition=nesi_research
#SBATCH --hint=nomultithread
#SBATCH --time=02:30:00
#SBATCH --output=smooth_sem_%j.out

# note: runtime is dependent on the size of the smoothing Gaussian, i.e.
#   larger radii = longer runtime. from experience, if sigma < 10km, set 
#   time ~= 30min for a start. add 30 min for every additional 5km to sigma.
#   trial and error required
# 
# ==============================================================================
# INPUT PARAMETERS
# ==============================================================================
# KERNEL: the suffix of the .bin files, e.g. proc000000_vs.bin would be 'vs'
# SGMAH: horizontal standard deviation of the Gaussian in meters
# SGMAV: vertical standard deviation of the Gaussian in meters
# > note: Gaussian full width = sqrt(8) * sigma
# DIR_IN: directory to look for input .bin files
# DIR_OUT: director to output the smoothed .bin files
# USE_GPU: use GPU acceleration to speed up calculations. not available on Maui

KERNEL="vs_kernel"
SGMAH=20000.
SGMAV=10000.
DIR_IN="SMOOTH/"
DIR_OUT=${DIR_IN}
USE_GPU=".false"

# ==============================================================================
# Example usage:
#
# $ cd path/to/specfem/workdir
# $ mkdir SMOOTH
# $ cd SMOOTH
# $ ln -s ../OUTPUT_FILES/DATABASES_MPI/*vs.bin .
# $ cd ..
# $ sbatch run_xsmooth_sem.sh
#
# ==============================================================================

# Get the number of processors from Par_file, ignore comments
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

echo "xsmooth_sem ${KERNEL} w/ sigma_h=${SGMAH}, sigma_v=${SGMAV}"
echo "${NPROC} processors, GPU option: ${USE_GPU}"
echo
echo "`date`"

# EXAMPLE CALL:
# srun -n NPROC xmooth_sem SIGMA_H SIGMA_V KERNEL_NAME INPUT_DIR OUTPUT_DIR USE_GPU
time srun -n ${NPROC} ./bin/xsmooth_sem ${SGMAH} ${SGMAV} ${KERNEL} ${DIR_IN} ${DIR_OUT} ${USE_GPU}

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "finished at: `date`"

