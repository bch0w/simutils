#!/bin/bash -e

#PJM -N test_xsmooth
#PJM -L rscgrp=share-debug
#PJM -L gpu=4
#PJM --mpi proc=4
#PJM -L elapse=00:10:00
#PJM -g gc62

module load cuda
module load gcc
module load ompi-cuda

echo "running xsmooth_sem : `date`"
echo

# Kernel to smooth and smoothing parameters must be specified by user
SGMAH=20000. #10000.
SGMAV=10000. #7000.
DIR_IN="OUTPUT_FILES/DATABASES_MPI/"
DIR_OUT="OUTPUT_FILES/SMOOTH/"
NPROC=1
USE_GPU=".true."

mkdir -p ./OUTPUT_FILES/SMOOTH
rm -rf ./OUTPUT_FILES/SMOOTH/*

# EXAMPLE
# srun -n NPROC xmooth_sem SIGMA_H SIGMA_V KERNEL_NAME INPUT_DIR OUTPUT_DIR USE_GPU

KERNEL="rho"
mpiexec -machinefile $PJM_O_NODEINF -n \
        $PJM_MPI_PROC -npernode $NPROC ./bin/xsmooth_sem \
        $SGMAH $SGMAV $KERNEL $DIR_IN $DIR_OUT $USE_GPU 

echo
echo "checkpoint: `date`"
echo

KERNEL="vs"
mpiexec -machinefile $PJM_O_NODEINF -n \
        $PJM_MPI_PROC -npernode $NPROC ./bin/xsmooth_sem \
        $SGMAH $SGMAV $KERNEL $DIR_IN $DIR_OUT $USE_GPU

echo
echo "checkpoint: `date`"
echo

KERNEL="vp"
mpiexec -machinefile $PJM_O_NODEINF -n \
        $PJM_MPI_PROC -npernode $NPROC ./bin/xsmooth_sem \
        $SGMAH $SGMAV $KERNEL $DIR_IN $DIR_OUT $USE_GPU


echo
echo "finished at: `date`"
echo
