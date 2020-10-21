#!/bin/bash -e

#SBATCH --job-name=xgenerate_databases
#SBATCH --nodes=10
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=10
#SBATCH --hint=nomultithread
#SBATCH --clusters=maui
#SBATCH --account=nesi00263
#SBATCH --partition=nesi_research
#SBATCH --time=00:15:00
#SBATCH --output=generate_databases_%j.out

# Set options to enable OpenMP/MPI Hybryd
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=true
export OMP_PLACES=cores

# Set the compiler option
# COMPILER=SPECFEM3D/20190730-CrayCCE-19.04
# COMPILER=SPECFEM3D/20190730-CrayIntel-19.04
module load gcc/8.3.0
COMPILER=SPECFEM3D/20190730-CrayGNU-19.04

module load ${COMPILER}

# get the number of processors, ignoring comments in the Par_file
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

# Make the Database directory 
mkdir -p ${BASEMPIDIR}

# This is a MPI simulation
echo ${COMPILER}
echo "xgenerate_databases ${NPROC} processors"
echo
echo "`date`"
time srun -n ${NPROC} xgenerate_databases

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "finished at: `date`"
