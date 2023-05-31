#!/bin/sh                                                                        
                                                                                 
#SBATCH --job-name=combine_vol_data_vtk
#SBATCH --ntasks=1
#SBATCH --partition=debug                                                        
#SBATCH --time=00:05:00                                                          
#SBATCH --output=combine_vol_data_vtk_%j.out
                                                                                 
                                                                                 
ulimit -s unlimited                                                              
ulimit -l unlimited    

# Rubric:
# * This exectuable needs to be run inside a SPECFEM3D working directory
# * You will need to symlink (preferred) or copy the DATABASE binary files
#   that you want to visualize into a the $DIR_IN directory, which by default
#   is named SUM. 
# * Once you run this batch script, a .vtk file will be generated inside 
#   $DIR_OUT, by default the same SUM directory. 
# * You can symlink multiple parameters, e.g. Vp, Vs, Rho, and .vtk files
#   will be generated for each of these parameters
#
# An example of this workflow is as follows:
#
# $ cd /path/to/specfem_work_dir
# $ mkdir SUM
# $ cd SUM
# $ ln -s ../OUTPUT_FILES/DATABASES_MPI/proc*_vs.bin .
# $ ln -s ../OUTPUT_FILES/DATABASES_MPI/proc*_vp.bin .
# $ cd ..
# $ sbatch run_xcombine_vol_data_vtk.sh

# Example call for the SPECFEM binary
# srun -n nproc xcombine... proc_start proc_end kernel dir_in dir_out hi_res

# Quantity needs to be specified by the user, otherwise all things will be run`
QUANTITY=$1

# Resolution of the resultant VTK file. 
# 0 for low-res, outputting points at the element corners
# 1 for hi-res, outputting points for each GLL point, takes longer and 
#   much larger file sizes
RES=$2

# Dynamically get the number of processors from the Par_file
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
NPROC_START=0
NPROC_END=`expr $NPROC - 1`

# Set the paths for Specfem to search
DIR_IN="./SUM"
DIR_OUT=${DIR_IN}
mkdir -p ${DIR_OUT}

# If no quantity is specified, grab everything in the sum folder
if [ -z "$1" ]
then
    for f in ${DIR_IN}/proc000000_*.bin;
    do
        # Cut out the kernel quantity name from the filename
        QUANTITY=`echo ${f:17} | cut -d'.' -f 1`

        # Run the binary
        echo "xcombine_vol_data_vtk ${NPROC_START} ${NPROC_END} for ${QUANTITY}"
        echo
        echo "`date`"
        time srun -n 1 ./bin/xcombine_vol_data_vtk \
                ${NPROC_START} ${NPROC_END} ${QUANTITY} ${DIR_IN}/ ${DIR_OUT}/ ${RES}
        echo
        echo "finished at: `date`"
    done
# Else just run the quantity that was specified
else
    echo "xcombine_vol_data_vtk ${NPROC_START} ${NPROC_END} for ${QUANTITY}"
    echo
    echo "`date`"
    time srun -n 1 ./bin/xcombine_vol_data_vtk \
            ${NPROC_START} ${NPROC_END} ${QUANTITY} ${DIR_IN}/ ${DIR_OUT}/ ${RES}
    echo
    echo "finished at: `date`"
fi
