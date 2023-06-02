#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=debug
#SBATCH --time=00:10:00
#SBATCH --job-name="python"
#SBATCH --output="%j_python.log"

module purge
module load slurm

ulimit -l unlimited

python $1
