#!/bin/bash -e
#SBATCH --job-name=python_job
#SBATCH --time=02:00:00
#SBATCH --account=nesi00263
#SBATCH --clusters=maui_ancil
#SBATCH --partition=nesi_prepost
#SBATCH --cpus-per-task=1
#ABATCH --ntasks=1

module load Anaconda3/5.2.0-GCC-7.1.0
source activate tomo

python set_values_above_topo.py

echo "finished"


