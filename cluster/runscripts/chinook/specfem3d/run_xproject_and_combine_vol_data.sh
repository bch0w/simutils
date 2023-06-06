#!/bin/bash

#SBATCH --job-name=xproject
#SBATCH --ntasks=80
#SBATCH --partition=t1small
#SBATCH --time=00:10:00
#SBATCH --output=project_and_combine_%j.out

ulimit -s unlimited
ulimit -l unlimited
umask 022

INPUT_DIR="OUTPUT_FILES/SMOOTHED_DATABASES_MPI/"
OUTPUT_DIR="OUTPUT_FILES/"
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

for DATA_FILENAME in vp vs vp rho qmu qkappa
do
    echo "xproject_and_combine_vol_data_on_regular_grid ${DATA_FILENAME} on ${NPROC} processors"
    time mpiexec -n ${NPROC} ./bin/xproject_and_combine_vol_data_on_regular_grid ${DATA_FILENAME} ${INPUT_DIR}/ ${OUTPUT_DIR}/
done


