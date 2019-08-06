#!/bin/bash -e

#SBATCH --job-name=smooth_kernel
#SBATCH --nodes=4
#SBATCH --ntasks=144
#SBATCH --cpus-per-task=1
#SBATCH --account=nesi00263
#SBATCH --clusters=maui
#SBATCH --partition=nesi_research
#SBATCH --hint=nomultithread
#SBATCH --time=0:55:00
#SBATCH --output=smooth.log

# EXAMPLE CALL:
# sbatch simutils/run_templates/smooth_kernel.sh 2018p130600 hess_kernel
# default time was 55 minutes but this timed out for a few kernel smooths so
# time was bumped up to 1:30 for 5000x1000 smoothing

KERNEL=$1
if [ -z "$1" ]
then
	echo "KERNEL NEEDS TO BE *_kernel"
	exit
fi

echo "`date`"

SGMAH=40000.
SGMAV=8000.
DIR_IN="INPUT_SMOOTH/"

echo "smoothing ${KERNEL} w/ sigma_h=${SGMAH}, sigma_v=${SGMAV}"
echo
currentdir=`pwd`

# NOTE smoothing time for alpha_kernel with CPU, 144 tasks, 4 nodes
#	   was ~3000 sec or 50 min

# EXAMPLE CALL:
# srun -n NPROC ./bin/xmooth_sem SIGMA_H SIGMA_V KERNEL_NAME INPUT_DIR OUTPUT_DIR USE_GPU
srun -n 144 ./bin/xsmooth_sem ${SGMAH} ${SGMAV} ${KERNEL} ${DIR_IN} ${DIR_IN} .false

echo
echo "done `date`"
