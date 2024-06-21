#!/bin/bash -e

#SBATCH --job-name=xsmooth_laplacian_sem
#SBATCH --ntasks=40
#SBATCH --partition=debug
#SBATCH --time=00:05:00
#SBATCH --output=smooth_lap_sem_%j.out

# ==============================================================================
# INPUT PARAMETERS
# ==============================================================================
# KERNEL: the suffix of the .bin files, e.g. proc000000_vs.bin would be 'vs'
# SGMAH: horizontal standard deviation of the Gaussian in km
# SGMAV: vertical standard deviation of the Gaussian in km
# > note: Gaussian full width = sqrt(8) * sigma
# DIR_IN: directory to look for input .bin files
# DIR_OUT: director to output the smoothed .bin files

KERNEL=$1
SGMAH=10.
SGMAV=5.
DIR_IN="DATABASES_MPI/"
DIR_OUT=${DIR_IN}

# ==============================================================================
# Example usage:
#
# $ cd path/to/specfem/workdir
# $ mkdir SMOOTH
# $ cd SMOOTH
# $ ln -s ../OUTPUT_FILES/DATABASES_MPI/*vs.bin .
# $ cd ..
# $ sbatch run_xsmooth_sem.sh
#
# ==============================================================================

# Get the number of processors from Par_file, ignore comments
ulimit -s unlimited                                                              
ulimit -l unlimited                                                              
umask 022                                                                        
                                                                                 
                                                                                 
BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `                     
                                                                                 
# script to run the mesher and the solver                                        
# read DATA/Par_file to get information about the run                            
# compute total number of nodes needed                                           
NPROC_XI=`grep ^NPROC_XI DATA/Par_file | cut -d = -f 2 `                         
NPROC_ETA=`grep ^NPROC_ETA DATA/Par_file | cut -d = -f 2`                        
NCHUNKS=`grep ^NCHUNKS DATA/Par_file | cut -d = -f 2 `                           
                                                                                 
# total number of nodes is the product of the values read                        
NPROC=$(( $NCHUNKS * $NPROC_XI * $NPROC_ETA ))  

echo "xsmooth_laplacian_sem ${KERNEL} w/ sigma_h=${SGMAH}, sigma_v=${SGMAV}"
echo
echo "`date`"

# EXAMPLE CALL:
# srun -n NPROC xmooth_laplacian_sem SIGMA_H SIGMA_V KERNEL_NAME INPUT_DIR OUTPUT_DIR 
time mpirun -n ${NPROC} ./bin/xsmooth_laplacian_sem ${SGMAH} ${SGMAV} ${KERNEL} ${DIR_IN} ${DIR_OUT} 

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "finished at: `date`"

