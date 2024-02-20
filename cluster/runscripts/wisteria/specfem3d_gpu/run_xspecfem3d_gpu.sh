#!/bin/bash -e

#PJM -N test_specfem
#PJM -L rscgrp=share-debug
#PJM -L gpu=2
#PJM --mpi proc=2
#PJM -L elapse=00:10:00
#PJM -g gc62

module load cuda/12.2
module load gcc
module load ompi-cuda

# get the number of processors, ignoring comments in the Par_file
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2 | cut -d \# -f 1`
echo "The simulation will run on NPROC = " $NPROC " MPI tasks"

# runs simulation
echo
echo "  running solver on $NPROC processors..."
echo "start from `date`"
echo

mpiexec -machinefile $PJM_O_NODEINF -n \
        $PJM_MPI_PROC -npernode $NPROC ./bin/xspecfem3D

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "finished at `date`"

