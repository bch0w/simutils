#!/bin/sh

#SBATCH --job-name=xgenerate_databases
#SBATCH --ntasks=40
#SBATCH --partition=t2small
#SBATCH --time=00:02:00
#SBATCH --output=generate_databases_%j.out

ulimit -s unlimited
ulimit -l unlimited
umask 022

# get the number of processors, ignoring comments in the Par_file
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

# Make the Database directory
mkdir -p ${BASEMPIDIR}

# This is a MPI simulation
echo "`date`"
echo "xmeshfem3D ${NPROC} processors"
echo
time mpiexec -n ${NPROC} ./bin/xmeshfem3D
echo "xgenerate_databases ${NPROC} processors"
echo
time mpiexec -n ${NPROC} ./bin/xgenerate_databases

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "finished at: `date`"

