#!/bin/bash -e

#SBATCH --job-name=xgenerate_databases
#SBATCH --nodes=4
#SBATCH --ntasks=144
#SBATCH --cpus-per-task=1
#SBATCH --clusters=maui
#SBATCH --account=nesi00263
#SBATCH --partition=nesi_research
#SBATCH --hint=nomultithread
#SBATCH --time=00:15:00
#SBATCH --output=generate_databases.log


# for nz_north tomo files, 8 minute gen db
# for nz_x1200_y600, 12 minute gen db
echo "generating databases: `date`"
currentdir=`pwd`

# get the number of processors, ignoring comments in the Par_file
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `
mkdir -p $BASEMPIDIR

# This is a MPI simulation
echo
echo "  running database generation on $NPROC processors..."
echo
srun -n $NPROC ./bin/xgenerate_databases

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo "finished at: `date`"
