#!/bin/bash -e

#SBATCH --job-name=mesh_and_db
#SBATCH --nodes=4
#SBATCH --ntasks=144
#SBATCH --cpus-per-task=1
#SBATCH --account=nesi00263
#SBATCH --partition=nesi_research
#SBATCH --hint=nomultithread
#SBATCH --time 00:08:00

echo "decomposing mesh and generating databases: `date`"
currentdir=`pwd`

# get the number of processors, ignoring comments in the Par_file
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `
mkdir -p $BASEMPIDIR

# decomposes mesh using the pre-saved mesh files in MESH-default
echo
echo "  decomposing mesh..."
echo
./bin/xdecompose_mesh $NPROC ./MESH $BASEMPIDIR
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# runs database generation
if [ "$NPROC" -eq 1 ]; then
	# This is a serial simulation
	echo
	echo "  running database generation..."
	echo
	./bin/xgenerate_databases
else
	# This is a MPI simulation
	echo
	echo "  running database generation on $NPROC processors..."
	echo
	srun -n $NPROC ./bin/xgenerate_databases
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo "finished at: `date`"
