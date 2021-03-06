#!/bin/bash -e

#SBATCH --job-name=xmeshfem3D
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=4
#SBATCH --hint=nomultithread
#SBATCH --clusters=maui
#SBATCH --account=nesi00263
#SBATCH --partition=nesi_research
#SBATCH --time=00:00:45
#SBATCH --output=meshfem3D_%j.out

# Set options to enable OpenMP/MPI Hybryd
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=true
export OMP_PLACES=cores

# Get the number of processors from the Par_file, ignore comments
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

# Make the Database directory
mkdir -p ${BASEMPIDIR}

echo "xmeshfem3D on ${NPROC} processors"
echo
echo "`date`"
time srun -n ${NPROC} ./bin/xmeshfem3D
echo
echo "finished at: `date`"

