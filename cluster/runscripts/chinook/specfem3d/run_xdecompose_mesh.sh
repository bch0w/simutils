#!/bin/sh 
 
#SBATCH --job-name=xdecompose_mesh
#SBATCH --ntasks=1 
#SBATCH --tasks-per-node=24
#SBATCH --partition=t2small
#SBATCH --time=00:00:45 
#SBATCH --output=decompose_mesh_%j.out

ulimit -s unlimited
ulimit -l unlimited

# Set the directory to search for external mesh files
MESH="./MESH"

# Get the number or processors and Database directory form the Par_file
# ignore comments in the line
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

# Make the Database directory
mkdir -p ${BASEMPIDIR}

# Decomposes mesh using files contained in ./MESH
echo "xdecompose_mesh"
echo
echo "`date`"
time ./bin/xdecompose_mesh ${NPROC} ${MESH} ${BASEMPIDIR}
echo
echo "finished at: `date`"

