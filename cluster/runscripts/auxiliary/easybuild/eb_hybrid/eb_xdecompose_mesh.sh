#!/bin/bash -e

#SBATCH --job-name=xdecompose_mesh
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --clusters=maui
#SBATCH --account=nesi00263
#SBATCH --partition=nesi_research
#SBATCH --time=00:01:00
#SBATCH --output=decompose_mesh_%j.out

# Set the compiler option
# COMPILER=SPECFEM3D/20190730-CrayCCE-19.04
# COMPILER=SPECFEM3D/20190730-CrayIntel-19.04
module load gcc/8.3.0
COMPILER=SPECFEM3D/20190730-CrayGNU-19.04

module load ${COMPILER}

# Set the directory to search for external mesh files
MESH="./MESH"

# Get the number or processors and Database directory form the Par_file
# ignore comments in the line
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

# Make the Database directory 
mkdir -p ${BASEMPIDIR}

# Decomposes mesh using files contained in ./MESH
echo ${COMPILER}
echo "xdecompose_mesh"
echo
echo "`date`"
xdecompose_mesh ${NPROC} ${MESH} ${BASEMPIDIR}
echo
echo "finished at: `date`"
