#!/bin/bash -e

#SBATCH --job-name=xspecfem3d
#SBATCH --nodes=4
#SBATCH --ntasks=144
#SBATCH --cpus-per-task=1
#SBATCH --clusters=maui
#SBATCH --account=nesi00263
#SBATCH --partition=nesi_research
#SBATCH --time 00:10:00
#SBATCH --output=simulation.log

echo "running simulation: `date`"
currentdir=`pwd`

NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `
mkdir -p $BASEMPIDIR

# runs simulation
if [ "$NPROC" -eq 1 ]; then
	# This is a serial simulation
	echo
	echo "  running solver..."
	echo
	./bin/xspecfem3D
else
	# This is a MPI simulation
	echo
	echo "  running solver on $NPROC processors..."
	echo
	time srun -n $NPROC ./bin/xspecfem3D
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "done"
echo `date`

