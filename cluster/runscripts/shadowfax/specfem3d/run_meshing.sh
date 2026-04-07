#!/bin/bash

#SBATCH --job-name=xmeshfem3D
#SBATCH --ntasks=32
#SBATCH --partition=defq
#SBATCH --time=00:20:00
#SBATCH --output=mesher_%j.out

# Get the number of processors from the Par_file, ignore comments
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

# Remake the output directories assuming we are starting fresh
rm -rf OUTPUT_FILES
mkdir -p OUTPUT_FILES

rm -rf ${BASEMPIDIR}
mkdir -p ${BASEMPIDIR}

echo "xmeshfem3D on ${NPROC} processors"
echo
echo "`date`"
time mpiexec -n ${NPROC} ./bin/xmeshfem3D
echo
echo "finished at: `date`"

# This is a MPI simulation
echo "`date`"
echo "xgenerate_databases ${NPROC} processors"
echo
time mpiexec -n ${NPROC} ./bin/xgenerate_databases

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "finished at: `date`"

