#!/bin/sh                                                                        
                                                                                 
#SBATCH --job-name=write_stations_file
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=24
#SBATCH --partition=debug
#SBATCH --time=00:05:00
#SBATCH --output=write_stations_file_%j.out
                                                                                 
ulimit -s unlimited                                                              
ulimit -l unlimited 

# User-defined parameters
OLD_FILE='./DATA/STATIONS'
NEW_FILE='./DATA/STATIONS_C2S'
ROT_FILE='./DATA/rotations_nu'
TOPO='.true.'
ELLP='.true.'


# Get the number or processors and Database directory form the Par_file
# ignore comments in the line
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

# Make the Database directory
mkdir -p ${BASEMPIDIR}

# Decomposes mesh using files contained in ./MESH
echo "write_stations_file"
echo
echo "`date`"
# write_cmt_solution_file OLD_FILE NEW_FILE ROT_FILE TOPOGRAPHY (bool) ELLIPTICITY (bool)
time ./bin/write_stations_file ${OLD_FILE} ${NEW_FILE} ${TOPO} ${ELLP}
echo
echo "finished at: `date`"

