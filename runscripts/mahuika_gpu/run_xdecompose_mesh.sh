#!/bin/bash -e

#SBATCH --job-name=xdecompose_mesh
#SBATCH --clusters=mahuika
#SBATCH --account=nesi00263
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=100G
#SBATCH --time=00:02:00
#SBATCH --output=decompose_mesh_%j.out

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

