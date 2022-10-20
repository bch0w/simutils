#!/bin/bash

#SBATCH --job-name=xmeshfem3D
#SBATCH --ntasks=4
#SBATCH --partition=debug
#SBATCH --time=00:00:45
#SBATCH --output=meshfem3D_%j.out

ulimit -s unlimited
ulimit -l unlimited
umask 022

# Get the number of processors from the Par_file, ignore comments
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

# Make the Database directory
mkdir -p ${BASEMPIDIR}
mkdir -p MESH

echo "xmeshfem3D on ${NPROC} processors"
echo
echo "`date`"
time mpiexec -n ${NPROC} ./bin/xmeshfem3D
echo
echo "finished at: `date`"

