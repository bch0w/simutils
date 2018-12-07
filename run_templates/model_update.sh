#!/bin/bash

#SBATCH --job-name=model_update
#SBATCH --nodes=4
#SBATCH --ntasks=144
#SBATCH --cpus-per-task=1
#SBATCH --account=nesi00263
#SBATCH --partition=nesi_research
#SBATCH --hint=nomultithread
#SBATCH --time 0:01:00

srun -n 144 ./bin/xmodel_update 0.03

echo
echo "done"


