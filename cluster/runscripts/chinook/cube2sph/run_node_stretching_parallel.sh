#!/bin/sh

#SBATCH --job-name=node_stretch
#SBATCH --ntasks=56
#SBATCH --tasks-per-node=28
#SBATCH --partition=t2small
#SBATCH --time=00:30:00
#SBATCH --output=node_stretch_%j.out


ulimit -s unlimited
ulimit -l unlimited


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

