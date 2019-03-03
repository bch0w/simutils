#!/bin/bash -e

#SBATCH --job-name=smooth_alpha
#SBATCH --nodes=4
#SBATCH --ntasks=144
#SBATCH --cpus-per-task=1
#SBATCH --account=nesi00263
#SBATCH --partition=nesi_research
#SBATCH --hint=nomultithread
#SBATCH --time 0:55:00

# EXAMPLE CALL:
# sbatch simutils/run_templates/smooth_kernel.sh 2018p130600 hess_kernel
# default time was 55 minutes but this timed out for a few kernel smooths so
# time was bumped up to 1:30


echo "`date`"

SGMAH=5000.
SGMAV=1000.
DIR_IN="OUTPUT_SUM"
KERNEL="alpha_kernel"


echo "smoothing ${KERNEL} for w/ sigma_h=${SGMAH}, sigma_v=${SGMAV}"
echo
currentdir=`pwd`


# EXAMPLE CALL:
# srun -n NPROC ./bin/xmooth_sem SIGMA_H SIGMA_V KERNEL_NAME INPUT_DIR OUTPUT_DIR USE_GPU
srun -n 144 ./bin/xsmooth_sem ${SGMAH} ${SGMAV} ${KERNEL} ${DIR_IN} ${DIR_IN} .false

echo
echo "done `date`"
