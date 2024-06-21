#!/bin/bash

#SBATCH --job-name=xinterpolate_model
#SBATCH --ntasks=40
#SBATCH --partition=t1small
#SBATCH --time=00:05:00
#SBATCH --output=interpolate_%j.out

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


# See src/tomography/postprocess_sensitivity_kernels/interpolate_model.F90
# for model paramter input descriptions
MODEL_PARAMETER=${1}
# Mesh that you want to interpolate
OLD_MODEL_DIR="MODEL_OLD/"  
OLD_TOPO_DIR="MODEL_OLD/"  # location of: proc*_reg?_solver_data.bin
# Mesh that you will interpolate onto
NEW_MODEL_DIR="MODEL_NEW/"
NEW_MODEL_TOPO="MODEL_NEW/"

# MIDPOINT: 1: yes, good for NEX_old=NEX_new; 0: no, good for NEX_old!=NEX_new
MIDPOINT_SEARCH=1 

echo "xinterpolate_model on ${NPROC} processors"
echo
echo "`date`"
time mpiexec -n ${NPROC} ./bin/xinterpolate_model ${MODEL_PARAMTER} ${OLD_MODEL_DIR} ${OLD_TOPO_DIR} ${NEW_MODEL_DIR} ${NEW_MODEL_TOPO} ${MIDPOINT_SEARCH}
echo
echo "finished at: `date`"

