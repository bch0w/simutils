#!/bin/bash -e

#SBATCH --job-name=${EVENT_ID}_fwd
#SBATCH --nodes=4
#SBATCH --ntasks=144
#SBATCH --cpus-per-task=1
#SBATCH --account=nesi00263
#SBATCH --partition=nesi_research
#SBATCH --hint=nomultithread
#SBATCH --time 01:15:00

echo "running simulation: `date`"
currentdir=`pwd`

# stores setup
cp DATA/Par_file OUTPUT_FILES/
cp DATA/CMTSOLUTION OUTPUT_FILES/CMTSOLUTION
cp RUNFORWARD_*.sh OUTPUT_FILES/

# get the number of processors, ignoring comments in the Par_file
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

