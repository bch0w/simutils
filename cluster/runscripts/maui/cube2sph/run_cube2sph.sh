#!/bin/bash -e

#SBATCH --job-name=cube2sph
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1
#SBATCH --clusters=maui
#SBATCH --account=gns03247
#SBATCH --partition=nesi_research
#SBATCH --hint=nomultithread
#SBATCH --time=00:02:00
#SBATCH --output=cube2sph_%j.out


# USER SET PARAMETERS, TRIAL AND ERROR ON ROTATE
# NZ North Extended trial
# CENTER_LAT=-37.95
# CENTER_LON=175.5
# ROTATE_ANGLE=10.0
# Alaska TLiu example
CENTER_LAT=62.5
CENTER_LON=-151.0
ROTATE_ANGLE=20.0

# Get the number or processors and Database directory form the Par_file
# ignore comments in the line
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

# Make the Database directory
mkdir -p ${BASEMPIDIR}

# Decomposes mesh using files contained in ./MESH
echo "cube2sph on ${NPROC} processors"
echo
echo "`date`"
time srun -n ${NPROC} ./bin/cube2sph ${CENTER_LAT} ${CENTER_LON} ${ROTATE_ANGLE}
echo
echo "finished at: `date`"

