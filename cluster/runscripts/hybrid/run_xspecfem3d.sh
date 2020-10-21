#!/bin/bash -e

#SBATCH --job-name=xspecfem3D
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=4
#SBATCH --hint=nomultithread
#SBATCH --clusters=maui
#SBATCH --account=nesi00263
#SBATCH --partition=nesi_research
#SBATCH --time 00:05:00
#SBATCH --output=specfem3D_%j.out

# Set options to enable OpenMP/MPI Hybryd
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=true
export OMP_PLACES=cores

# Get the number of processors from Par_file, ignore comments
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

# Make the Database directory
mkdir -p $BASEMPIDIR

# This is a MPI simulation
echo "xspecfem3d ${NPROC} processors"
echo
time srun -n ${NPROC} ./bin/xspecfem3D

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "finished at: `date`"

