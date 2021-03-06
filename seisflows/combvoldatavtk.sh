#!/bin/bash

#SBATCH --job-name=combine_vol_data_vtk
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=nesi00263
#SBATCH --clusters=maui
#SBATCH --partition=nesi_research
#SBATCH --time=0:15:00
#SBATCH --hint=nomultithread
#SBATCH --output=comb.log

# Example Call
# srun -n nproc xcombine... proc_start proc_end kernel dir_in dir_out hi_res

# Quantity needs to be specified by the user, otherwise all things will be run`
QUANTITY=$1

# Dynamically get the number of processors from the Par_file
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
NPROC_START=0
NPROC_END=`expr $NPROC - 1`

# Set the paths for Specfem to search
DIR_IN="./SUM"
DIR_OUT="./SUM"
mkdir -p ${DIR_OUT}

# If no quantity is specified, grab everything in the sum folder
if [ -z "$1" ]
then
    for f in ${DIR_IN}/proc000000_*.bin;
    do
        # Cut out the kernel quantity name from the filename by removing 
        # the proc000000_ portion of the 
        QUANTITY=`echo ${f:17} | cut -d'.' -f 1`

        # Run the binary
        echo "xcombine_vol_data_vtk ${NPROC_START} ${NPROC_END} for ${QUANTITY}"
        echo
        echo "`date`"
        time srun -n 1 ./bin/xcombine_vol_data_vtk \
                ${NPROC_START} ${NPROC_END} ${QUANTITY} ${DIR_IN}/ ${DIR_OUT}/ 0
        echo
        echo "finished at: `date`"
    done
# Else just run the quantity that was specified
else
    echo "xcombine_vol_data_vtk ${NPROC_START} ${NPROC_END} for ${QUANTITY}"
    echo
    echo "`date`"
    time srun -n 1 ./bin/xcombine_vol_data_vtk \
            ${NPROC_START} ${NPROC_END} ${QUANTITY} ${DIR_IN}/ ${DIR_OUT}/ 0
    echo
    echo "finished at: `date`"
fi
