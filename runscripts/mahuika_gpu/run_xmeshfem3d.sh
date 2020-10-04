#!/bin/bash -e

#SBATCH --job-name=xmeshfem3D
#SBATCH --clusters=mahuika
#SBATCH --account=nesi00263
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=100G
#SBATCH --time=00:00:45
#SBATCH --output=meshfem3D_%j.out

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

