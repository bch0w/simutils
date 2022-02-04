#!/bin/bash -e

#SBATCH --job-name=setup_model_cartesian
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1
#SBATCH --clusters=maui
#SBATCH --account=gns03247
#SBATCH --partition=nesi_research
#SBATCH --hint=nomultithread
#SBATCH --time=00:10:00
#SBATCH --output=setup_model_cartesian_%j.out


# Get the number or processors and Database directory form the Par_file
# ignore comments in the line
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

# Make the Database directory
mkdir -p ${BASEMPIDIR}

# Decomposes mesh using files contained in ./MESH
echo "setup_model_cartesian on ${NPROC} processors"
echo "changing the velocity model on GLL points"
echo
echo "`date`"
time srun -n ${NPROC} ./bin/setup_model_cartesian
echo
echo "finished at: `date`"

