#!/bin/sh

#SBATCH --job-name=axisem3d
#SBATCH --ntasks=64
#SBATCH --mem=50G
#SBATCH --time=01:00:00
#SBATCH --output=log_%A.out

mpiexec -n 64 ./axisem3d
