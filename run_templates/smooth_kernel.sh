#!/bin/bash -e

#SBATCH --job-name=smooth_kernel
#SBATCH --nodes=4
#SBATCH --ntasks=144
#SBATCH --cpus-per-task=1
#SBATCH --account=nesi00263
#SBATCH --partition=nesi_research
#SBATCH --hint=nomultithread
#SBATCH --time 1:15:00

# EXAMPLE CALL:
# sbatch simutils/run_templates/smooth_kernel.sh 2018p130600 hess_kernel
# default time was 55 minutes but this timed out for a few kernel smooths so
# time was bumped up to 1:30

EVENT_ID=$1
if [ -z "$1" ]
then
	echo "EVENT ID REQUIRED"
	exit
fi
KERNEL=$2
if [ -z "$2" ]
then
	echo "KERNEL NEEDS TO BE alpha_kernel, beta_kernel, rho_kernel or hess_kernel"
	exit
fi

echo "`date`"


SGMAH=5000.
SGMAV=1000.
DIR_IN="INPUT_KERNELS/${EVENT_ID}"

echo "smoothing ${KERNEL} for ${EVENT_ID} w/ sigma_h=${SGMAH}, sigma_v=${SGMAV}"
echo
currentdir=`pwd`

# NOTE smoothing time for alpha_kernel with CPU, 144 tasks, 4 nodes
#	   was ~3000 sec or 50 min

# EXAMPLE CALL:
# srun -n NPROC ./bin/xmooth_sem SIGMA_H SIGMA_V KERNEL_NAME INPUT_DIR OUTPUT_DIR USE_GPU
srun -n 144 ./bin/xsmooth_sem ${SGMAH} ${SGMAV} ${KERNEL} ${DIR_IN} ${DIR_IN} .false

echo
echo "done `date`"
