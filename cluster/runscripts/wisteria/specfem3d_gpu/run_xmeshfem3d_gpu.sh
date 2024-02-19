#!/bin/bash -e

#PJM -N test_meshfem
#PJM -L rscgrp=share-debug
#PJM -L gpu=2
#PJM --mpi proc=2
#PJM -L elapse=00:01:00
#PJM -g gc62

module load cuda/12.2
module load gcc
module load ompi-cuda

echo "running example: `date`"
currentdir=`pwd`

# sets up directory structure in current example directory
echo
echo "   setting up example..."
echo

# cleans output files
mkdir -p OUTPUT_FILES
rm -rf OUTPUT_FILES/*

# stores setup
cp DATA/Par_file OUTPUT_FILES/
cp DATA/CMTSOLUTION OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/

BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `
mkdir -p $BASEMPIDIR


# get the number of processors, ignoring comments in the Par_file
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2 | cut -d \# -f 1`
echo "The simulation will run on NPROC = " $NPROC " MPI tasks"

#
echo
echo " xmeshfem3D ..."
echo

mpiexec -machinefile $PJM_O_NODEINF -n \
	$PJM_MPI_PROC -npernode $NPROC ./bin/xmeshfem3D 

echo "finished at `date`"

