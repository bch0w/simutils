#!/bin/bash -e

#SBATCH --job-name=xspecfem3D
#SBATCH --nodes=10
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=10
#SBATCH --clusters=maui
#SBATCH --account=nesi00263
#SBATCH --partition=nesi_research
#SBATCH --time 00:02:00
#SBATCH --output=specfem3D_%j.out

# Set the compiler option
# COMPILER=SPECFEM3D/20190730-CrayCCE-19.04
# COMPILER=SPECFEM3D/20190730-CrayIntel-19.04
module load gcc/8.3.0
COMPILER=SPECFEM3D/20190730-CrayGNU-19.04

module load ${COMPILER}

# Get the number of processors from Par_file, ignore comments
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

# Make the Database directory
mkdir -p $BASEMPIDIR

# This is a MPI simulation
echo ${COMPILER}
echo "xspecfem3d ${NPROC} processors"
echo
time srun -n ${NPROC} xspecfem3D

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "finished at: `date`"

