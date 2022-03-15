#!/bin/sh                                                                        
                                                                                 
#SBATCH --job-name=write_cmt_solution_file
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=24
#SBATCH --partition=debug
#SBATCH --time=00:05:00
#SBATCH --output=write_cmt_solution_file_%j.out
                                                                                 
ulimit -s unlimited                                                              
ulimit -l unlimited 

# User-defined parameters
OLD_FILE='./DATA/CMTSOLUTION'
NEW_FILE='./DATA/CMTSOLUTION_C2S'
TOPO='.true.'
ELLP='.true.'


# Get the number or processors and Database directory form the Par_file
# ignore comments in the line
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

# Make the Database directory
mkdir -p ${BASEMPIDIR}

# Decomposes mesh using files contained in ./MESH
echo "write_cmt_solution_file"
echo
echo "`date`"
# write_cmt_solution_file OLD_FILE NEW_FILE TOPOGRAPHY (bool) ELLIPTICITY (bool)
time ./bin/write_cmt_solution_file ${OLD_FILE} ${NEW_FILE} ${TOPO} ${ELLP}
echo
echo "finished at: `date`"

