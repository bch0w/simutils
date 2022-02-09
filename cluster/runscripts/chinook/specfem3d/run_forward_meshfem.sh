#!/bin/bash -e

#SBATCH --job-name=forward
#SBATCH --nodes=1
#SBATCH --ntasks=80
#SBATCH --cpus-per-task=1
#SBATCH --clusters=maui
#SBATCH --account=nesi00263
#SBATCH --partition=nesi_research
#SBATCH --time=03:00:00
#SBATCH --output=meshfem3D_%j.out

# Get the number of processors from the Par_file, ignore comments
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

# Make the Database directory
mkdir -p ${BASEMPIDIR}

echo "`date`"
echo "xmeshfem3D on ${NPROC} processors"
echo
time srun -n ${NPROC} ./bin/xmeshfem3D
echo
echo "xgenerate_databases ${NPROC} processors"
echo
time srun -n ${NPROC} ./bin/xgenerate_databases_meshfem
echo
echo "xspecfem3D ${NPROC} processors"
echo
time srun -n ${NPROC} ./bin/xspecfem3D
echo
echo "finished at: `date`"

