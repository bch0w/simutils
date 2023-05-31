#!/bin/bash -e

#SBATCH --job-name=xspecfem3D
#SBATCH --ntasks=80
#SBATCH --partition=t1small
#SBATCH --time 00:30:00
#SBATCH --output=specfem3D_%j.out

ulimit -s unlimited
ulimit -l unlimited
umask 022

# Get the number of processors from Par_file, ignore comments
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

# Make the Database directory
mkdir -p $BASEMPIDIR

# This is a MPI simulation
echo "xspecfem3d ${NPROC} processors"
echo
time mpiexec -n ${NPROC} ./bin/xspecfem3D

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "finished at: `date`"

