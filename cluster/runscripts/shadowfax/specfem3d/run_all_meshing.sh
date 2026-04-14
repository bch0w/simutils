#!/bin/bash -e

#SBATCH --job-name=simblastmesh
#SBATCH --ntasks=16
#SBATCH --partition=defq
#SBATCH --time=02:30:00
#SBATCH --output=LOGS/simblastmesh_%j.out

# Get the number of processors from the Par_file, ignore comments
NPROC=${SLURM_NTASKS}
DATABASES=OUTPUT_FILES/DATABASES_MPI
mkdir -p ${DATABASES}

echo "${NPROC} TASKS"
echo "starting at: `date`"

# Run xmeshfem3D if it hasn't been run yet
if [ ! -f "${DATABASES}/proc000000_Database" ]; then
	echo "======== Running xmeshfem3D ========"
	time mpiexec -n ${NPROC} ./bin/xmeshfem3D
fi

# Run xgenenerate_databases if it hasn't been run yet
if [ ! -f "${DATABASES}/proc000000_vs.bin" ]; then
	echo "======== Running xgenerate_databases ========"	
	time mpiexec -n ${NPROC} ./bin/xgenerate_databases
fi

echo "finished at: `date`"
