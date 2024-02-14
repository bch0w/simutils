#!/bin/bash -e

#SBATCH --job-name=xspecfem3D
#SBATCH --ntasks=30
#SBATCH --partition=t1small
#SBATCH --time 00:15:00
#SBATCH --output=specfem3D_%j.out

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
numnodes=$(( $NCHUNKS * $NPROC_XI * $NPROC_ETA ))                                
                                                                                 

# Make the Database directory
mkdir -p $BASEMPIDIR

# This is a MPI simulation
echo "xspecfem3d ${numnodes} processors"
echo
time mpirun -np ${numnodes} ./bin/xspecfem3D

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "finished at: `date`"

