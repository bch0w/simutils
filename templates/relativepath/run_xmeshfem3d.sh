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


BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `
mkdir -p $BASEMPIDIR

NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

# decomposes mesh using the pre-saved mesh files in MESH-default
echo
echo "  running xmeshfem3D..."
echo
srun -n $NPROC ./bin/xmeshfem3D

echo "finished at: `date`"
