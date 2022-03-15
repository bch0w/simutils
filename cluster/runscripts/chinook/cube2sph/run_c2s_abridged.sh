#!/bin/sh

#SBATCH --job-name=cube2sph
#SBATCH --ntasks=56
#SBATCH --tasks-per-node=28
#SBATCH --partition=t2small
#SBATCH --time=00:30:00
#SBATCH --output=c2_%j.out

module purge
module load slurm
module load data/HDF5/1.10.6-pic-intel-2019b   
module load data/netCDF/4.7.4-pic-intel-2019b
module load data/netCDF-Fortran/4.5.2-pic-intel-2019b 

ulimit -s unlimited
ulimit -l unlimited

# NOTE: You need to run this AFTER xdecompose_mesh has been run serially


# USER SET PARAMETERS, TRIAL AND ERROR ON ROTATE
# NZ North Extended trial
# CENTER_LAT=-37.95
# CENTER_LON=175.5
# ROTATE_ANGLE=10.0
# Alaska TLiu example
CENTER_LAT=67.5
CENTER_LON=-152.5
ROTATE_ANGLE=0.0

# Get the number or processors and Database directory form the Par_file
# ignore comments in the line
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `
REFDIR="./DATABASES_MPI_REF"

# Make the Database directory
mkdir -p ${BASEMPIDIR}

# Run Cube2Sph to bend the cartesian mesh
echo
echo "cube2sph on ${NPROC} processors"
time srun -n ${NPROC} ./bin/cube2sph ${CENTER_LAT} ${CENTER_LON} ${ROTATE_ANGLE}

# Run generate databases to create database files from the mesh
echo
echo "xgenerate_databases ${NPROC} processors"  
time srun -n ${NPROC} ./bin/xgenerate_databases 
echo

# New DATA directory for applying velocity models to new mesh
rm ./DATA
ln -s DATA_utils DATA

# REF directory is required by C2S programs
mkdir -p ${REFDIR}
mv ${BASEMPIDIR}/* ${REFDIR}




echo
echo "node_stretching_parallel on ${NPROC} processors" 
time srun -n ${NPROC} ./bin/node_stretching_parallel  
echo
echo "setup_model_cartesian on ${NPROC} processors" 
time srun -n ${NPROC} ./bin/setup_model_cartesian   
echo
rm ./DATA
ln -s DATA_c2d DATA
echo
echo "xgenerate_databases ${NPROC} processors"  
time srun -n ${NPROC} ./bin/xgenerate_databases 
echo
echo "finished at: `date`"
