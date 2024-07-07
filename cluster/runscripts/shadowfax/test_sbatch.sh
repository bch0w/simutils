#!/bin/sh

#SBATCH --job-name=test
#SBATCH --ntasks=1
#SBATCH --array=0-2
#SBATCH --partition=defq
#SBATCH --time=00:00:05
#SBATCH --output=%A_%a.out
#SBATCH --parsable

echo "hello"

lscpu
