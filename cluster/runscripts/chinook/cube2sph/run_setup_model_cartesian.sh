#!/bin/sh                                                                        
                                                                                 
#SBATCH --job-name=setup_model_cartesian
#SBATCH --ntasks=48                                                              
#SBATCH --tasks-per-node=24                                                      
#SBATCH --partition=debug                                                        
#SBATCH --time=00:010:00                                                          
#SBATCH --output=setup_model_cartesian_%j.out
                                                                                 
ulimit -s unlimited                                                              
ulimit -l unlimited 


# Get the number or processors and Database directory form the Par_file
# ignore comments in the line
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

# Make the Database directory
mkdir -p ${BASEMPIDIR}

# Decomposes mesh using files contained in ./MESH
echo "setup_model_cartesian on ${NPROC} processors"
echo "changing the velocity model on GLL points"
echo
echo "`date`"
time srun -n ${NPROC} ./bin/setup_model_cartesian
echo
echo "finished at: `date`"

