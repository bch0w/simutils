#!/bin/bash -e

#SBATCH --job-name=generate_dbs
#SBATCH --nodes=4
#SBATCH --ntasks=144
#SBATCH --cpus-per-task=1
#SBATCH --account=nesi00263
#SBATCH --partition=nesi_research
#SBATCH --hint=nomultithread
#SBATCH --time 00:07:15

echo "generating databases: `date`"
currentdir=`pwd`

# get the number of processors, ignoring comments in the Par_file
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

echo $NPROC

mkdir -p $BASEMPIDIR

srun -n $NPROC ./bin/xgenerate_databases

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo "finished at: `date`"
