#!/bin/bash -e

#SBATCH --job-name=node_stretch
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1
#SBATCH --clusters=maui
#SBATCH --account=gns03247
#SBATCH --partition=nesi_research
#SBATCH --hint=nomultithread
#SBATCH --time=00:02:00
#SBATCH --output=node_stretch_%j.out


# Get the number or processors and Database directory form the Par_file
# ignore comments in the line
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

# Make the Database directory
mkdir -p ${BASEMPIDIR}

# Decomposes mesh using files contained in ./MESH
echo "node_stretching_parallel on ${NPROC} processors"
echo "adding surface and interior topography"
echo
echo "`date`"
time srun -n ${NPROC} ./bin/node_stretching_parallel
echo
echo "finished at: `date`"

