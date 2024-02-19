#!/bin/bash -e

#PJM -N test_combine_vtk
#PJM -L rscgrp=share-debug
#PJM -L gpu=4
#PJM --mpi proc=4
#PJM -L elapse=0:02:00
#PJM -g gc62

module load cuda
module load gcc
module load ompi-cuda

# Example Call
# srun -n nproc xcombine... proc_start proc_end kernel dir_in dir_out hi_res

# Quality can either be low (0) or high (1)
QUALITY=0

# Dynamically get the number of processors from the Par_file
#NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
NPROC=2
NPROC_START=0
NPROC_END=1

# Set the paths for Specfem to search
DIR_IN="OUTPUT_FILES/DATABASES_MPI"
DIR_OUT="OUTPUT_FILES/DATABASES_MPI"

echo
QUANTITY="vs"
mpiexec -machinefile $PJM_O_NODEINF -n $PJM_MPI_PROC -npernode $NPROC ./bin/xcombine_vol_data_vtk $NPROC_START $NPROC_END $QUANTITY $DIR_IN/ $DIR_OUT/ $QUALITY
echo

echo
QUANTITY="vp"
mpiexec -machinefile $PJM_O_NODEINF -n $PJM_MPI_PROC -npernode $NPROC ./bin/xcombine_vol_data_vtk $NPROC_START $NPROC_END $QUANTITY $DIR_IN/ $DIR_OUT/ $QUALITY
echo

echo "finished at: `date`"
