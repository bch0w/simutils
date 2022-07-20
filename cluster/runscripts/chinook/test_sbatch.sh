#!/bin/sh

#SBATCH --job-name=test
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=24
#SBATCH --partition=debug
#SBATCH --time=00:00:05
#SBATCH --output=%j.out
#SBATCH --parsable

echo "hello"
