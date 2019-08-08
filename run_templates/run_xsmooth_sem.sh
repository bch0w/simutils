#!/bin/bash -e

#SBATCH --job-name=xsmooth_sem
#SBATCH --nodes=4
#SBATCH --ntasks=144
#SBATCH --cpus-per-task=1
#SBATCH --account=nesi00263
#SBATCH --clusters=maui
#SBATCH --partition=nesi_research
#SBATCH --hint=nomultithread
#SBATCH --time=0:30:00
#SBATCH --output=smooth.log

KERNEL="beta_kernel"

echo "`date`"

SGMAH=20000.
SGMAV=3500.
DIR_IN="SMOOTH/"

echo "smoothing ${KERNEL} w/ sigma_h=${SGMAH}, sigma_v=${SGMAV}"
echo
currentdir=`pwd`

# EXAMPLE CALL:
# srun -n NPROC ./bin/xmooth_sem SIGMA_H SIGMA_V KERNEL_NAME INPUT_DIR OUTPUT_DIR USE_GPU
srun -n 144 ./bin/xsmooth_sem ${SGMAH} ${SGMAV} ${KERNEL} ${DIR_IN} ${DIR_IN} .false

echo
echo "done `date`"
