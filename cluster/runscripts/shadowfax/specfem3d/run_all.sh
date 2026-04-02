#!/bin/bash

#SBATCH --job-name=xmeshfem3D
#SBATCH --ntasks=8
#SBATCH --partition=defq
#SBATCH --time=01:00:00
#SBATCH --output=meshfem3D_%j.out

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
time mpiexec -n ${NPROC} ./bin/xgenerate_databases
time mpiexec -n ${NPROC} ./bin/xspecfem3D
