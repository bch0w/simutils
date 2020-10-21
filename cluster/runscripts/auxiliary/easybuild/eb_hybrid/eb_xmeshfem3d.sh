#!/bin/bash -e

#SBATCH --job-name=xmeshfem3D
#SBATCH --nodes=2
#SBATCH --ntasks=88
#SBATCH --cpus-per-task=1
#SBATCH --clusters=maui
#SBATCH --account=nesi00263
#SBATCH --partition=nesi_research
#SBATCH --time=00:05:00
#SBATCH --output=meshfem3D_%j.out

# Set the compiler option
# COMPILER=SPECFEM3D/20190730-CrayCCE-19.04
# COMPILER=SPECFEM3D/20190730-CrayIntel-19.04
module load gcc/8.3.0
COMPILER=SPECFEM3D/20190730-CrayGNU-19.04

module load ${COMPILER}

# Get the number of processors from the Par_file, ignore comments
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

# Make the Database directory
mkdir -p ${BASEMPIDIR}

echo ${COMPILER}
echo "xmeshfem3D on ${NPROC} processors"
echo
echo "`date`"
time srun -n ${NPROC} xmeshfem3D
echo
echo "finished at: `date`"
