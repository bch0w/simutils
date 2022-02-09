#!/bin/bash -e

#SBATCH --job-name=xsem_model_slice
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1
#SBATCH --clusters=maui
#SBATCH --account=nesi00263
#SBATCH --partition=nesi_research
#SBATCH --time=01:00:00
#SBATCH --output=sem_model_slice_%j.out

XYZ_FILE="shallow.xyz"
MODEL_DIR="model_0017/"
DATA_NAME="vs"
OUTFILE="model_0017_vs.out"

echo "xsem_model_slice"
echo
echo "`date`"
time srun -n 40 xsem_model_slice ${XYZ_FILE} ${MODEL_DIR} ${DATA_NAME} ${OUTFILE}
echo
echo "finished at: `date`"

