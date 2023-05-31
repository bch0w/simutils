#!/bin/bash

#SBATCH --job-name=xproject
#SBATCH --ntasks=80
#SBATCH --partition=t1small
#SBATCH --time=00:02:00
#SBATCH --output=project_and_combine_%j.out

ulimit -s unlimited
ulimit -l unlimited
umask 022

DATA_FILENAME="vs"
INPUT_DIR="OUTPUT_FILES/SMOOTHED_DATABASES_MPI/"
OUTPUT_DIR="OUTPUT_FILES/"
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

echo "writing proj grid"
python write_proj_grid_fd.py

echo "xproject_and_combine_vol_data_on_regular_grid on ${NPROC} processors"
time mpiexec -n ${NPROC} ./bin/xproject_and_combine_vol_data_on_regular_grid ${DATA_FILENAME} ${INPUT_DIR}/ ${OUTPUT_DIR}/

echo "running convert and plot"
python convert_and_plot_projection.py

