#!/bin/bash

#SBATCH --job-name=test_homog_halfs
#SBATCH --nodes=4
#SBATCH --ntasks=144
#SBATCH --cpus-per-task=1
#SBATCH --account=nesi00263
#SBATCH --partition=NeSI
#SBATCH --hint=nomultithread
#SBATCH --time 02:15:00

echo "running kernel simulation: `date`"
currentdir=`pwd`

# get the number of processors, ignoring comments in the Par_file
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

# make sure Par_file set to backward
./utils/change_simulation_type.pl -b

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
sleep 2


echo "finished successfully"





