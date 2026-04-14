#!/bin/bash -e

#SBATCH --job-name=simblast
#SBATCH --ntasks=16
#SBATCH --partition=defq
#SBATCH --time=02:30:00
#SBATCH --output=LOGS/simblast_%j.out

# Establish some run parameters
# Expand CMTSOLUTION symlink which will be tagged with the event name
TAG=readlink -f DATA/CMTSOLUTION | grep -oP '[^/]*$' | grep -oP '(?<=_)[^_]*$'
NSTEP=`grep ^NSTEP DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
NPROC=${SLURM_NTASKS}
OUTPUT_FILES="OUTPUT_FILES"
DATABASES="${OUTPUT_FILES}/DATABASES_MPI"
RESULTS="RESULTS/${TAG}/"

mkdir -p ${DATABASES}
mkdir -p ${RESULTS}

# Copy over some data
cp DATA/CMTSOLUTION ${RESULTS}
cp DATA/Par_file ${RESULTS}

echo "RUNNING ${TAG} WITH ${NPROC} TASKS"
echo "starting at: `date`"

# Run xmeshfem3D if it hasn't been run yet
if [ ! -f "${DATABASES}/proc000000_Database" ]; then
	echo "======== Running xmeshfem3D ========"
	time mpiexec -n ${NPROC} ./bin/xmeshfem3D
	cp ${OUTPUT_FILES}/output_meshfem3D.txt ${RESULTS}
fi

# Run xgenenerate_databases if it hasn't been run yet
if [ ! -f "${DATABASES}/proc000000_vs.bin" ]; then
	echo "======== Running xgenerate_databases ========"	
	time mpiexec -n ${NPROC} ./bin/xgenerate_databases
	cp ${OUTPUT_FILES}/output_generate_databases.txt ${RESULTS}
fi

# Run the solver if no waveform files exist
echo "======== Running xspecfem3D ========"
if ! compgen -G "${OUTPUT_FILES}/*.semv" > /dev/null; then
	time mpiexec -n ${NPROC} ./bin/xspecfem3D
	cp ${OUTPUT_FILES}/output_solver.txt ${RESULTS}
	mv ${OUTPUT_FILES}/*semv ${RESULTS}
	echo "finished at: `date`"
fi

# Make ParaView movie files if needed
echo "======== Running xcreate_move_shakemap_AVS_DX_GMT ========"
if ! compgen -G "${OUTPUT_FILES}/moviedata*" > /dev/null; then
    ./bin/xcreate_movie_shakemap_AVS_DX_GMT <<- EOF
	2
	1
	${NSTEP}
	2
	4
	EOF
	mv ${OUTPUT_FILES}/moviedata* ${RESULTS}
	mv ${OUTPUT_FILES}/*inp ${RESULTS}
fi


