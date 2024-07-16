#!/bin/bash

#SBATCH --job-name=xmeshfem3D
#SBATCH --ntasks=1
#SBATCH --partition=debug
#SBATCH --time=00:15:00
#SBATCH --output=meshfem3D_%j.out

ulimit -s unlimited
ulimit -l unlimited
umask 022

module load intel

# Get the number of processors from the Par_file, ignore comments
BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `
NPROC_XI=`grep ^NPROC_XI DATA/Par_file | cut -d = -f 2 `
NPROC_ETA=`grep ^NPROC_ETA DATA/Par_file | cut -d = -f 2` 
NCHUNKS=`grep ^NCHUNKS DATA/Par_file | cut -d = -f 2 ` 
NPROC=$(( $NCHUNKS * $NPROC_XI * $NPROC_ETA )) 

# Make the Database directory
mkdir -p OUTPUT_FILES
mkdir -p ${BASEMPIDIR}

echo "xmeshfem3D on ${NPROC} processors"
echo
echo "`date`"
time mpiexec -n ${NPROC} ./bin/xmeshfem3D
echo
echo "finished at: `date`"

