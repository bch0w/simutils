#!/bin/bash

#SBATCH --job-name=???
#SBATCH --nodes=x
#SBATCH --ntasks=y
#SBATCH --cpus-per-task=z
#SBATCH --hint=nomultithread
#SBATCH --account=???
#SBATCH --clusters=???
#SBATCH --partition=???
#SBATCH --time=0:01:00
#SBATCH --output=???_%j.out

# From email correspondence with NeSI engineer Alex Pletzer:
#
# To set the above values of x, y, and z, ensure that:
# y * z == c * x
# where c == 40, and c is the the number of cores per node defined by 
# the computing architechture of the current NeSI systems
#
# An example configuration is: x=2, y=10, z=8
# If we set NPROC=40 in Par_file: one example configuration is x=10, y=40, z=10

# Set options to enable OpenMP/MPI Hybryd
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=true
export OMP_PLACES=cores

# Choose the compiler option, delete or comment out the others
COMPILER=SPECFEM3D/20190730-CrayCCE-19.04
COMPILER=SPECFEM3D/20190730-CrayGNU-19.04
COMPILER=SPECFEM3D/20190730-CrayIntel-19.04

# Load the compiler
module load ${COMPILER}

# Run your exectuable
echo ${COMPILER}
echo
echo "`date`"
srun <executable> [options]
echo
echo "finished at: `date`"


