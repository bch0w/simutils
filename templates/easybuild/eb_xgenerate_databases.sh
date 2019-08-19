#!/bin/bash -e

#SBATCH --job-name=xgenerate_databases
#SBATCH --nodes=4
#SBATCH --ntasks=144
#SBATCH --cpus-per-task=1
#SBATCH --clusters=maui
#SBATCH --account=nesi00263
#SBATCH --partition=nesi_research
#SBATCH --time=00:15:00
#SBATCH --output=generate_databases_%j.out

# Set the compiler option
COMPILER=SPECFEM3D/20190730-CrayCCE-19.04
# COMPILER=SPECFEM3D/20190730-CrayGNU-19.04
# COMPILER=SPECFEM3D/20190730-CrayIntel-19.04

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
