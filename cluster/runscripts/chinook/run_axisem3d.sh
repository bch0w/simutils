#!/bin/sh

#SBATCH --job-name=axisem3d
#SBATCH --ntasks=160
#SBATCH --nodes=4
#SBATCH --partition=t2standard
#SBATCH --time=10:00:00
#SBATCH --output=log_axisem3d_%j.out

ulimit -s unlimited
ulimit -l unlimited
umask 022

conda activate axisem3d
time mpiexec -n ${SLURM_NTASKS} ./axisem3d

echo "finished at `date`"
