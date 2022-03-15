#!/bin/bash -e

#SBATCH --job-name=wf_cube2sph
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1
#SBATCH --clusters=maui
#SBATCH --account=gns03247
#SBATCH --partition=nesi_research
#SBATCH --hint=nomultithread
#SBATCH --time=00:02:00
#SBATCH --output=wf_cube2sph_%j.out


# USER SET PARAMETERS, TRIAL AND ERROR ON ROTATE
MESH="./MESH"
CENTER_LAT=-37.95
CENTER_LON=175.5
ROTATE_ANGLE=60.0

# Get the number or processors and Database directory form the Par_file
# ignore comments in the line
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

# Make the Database directory
rm -rf ${BASEMPIDIR}
mkdir -p ${BASEMPIDIR}

echo "xdecompose_mesh"
echo
time ./bin/xdecompose_mesh ${NPROC} ${MESH} ${BASEMPIDIR}     
echo
echo "cube2sph on ${NPROC} processors"
echo
time srun -n ${NPROC} ./bin/cube2sph ${CENTER_LAT} ${CENTER_LON} ${ROTATE_ANGLE}
echo
echo "xgenerate_databases ${NPROC} processors"      
echo
time srun -n ${NPROC} ./bin/xgenerate_databases   
echo
cd 
echo "finished at: `date`"

