#!/bin/bash

#SBATCH --job-name=python
#SBATCH --ntasks=1
#SBATCH --partition=debug
#SBATCH --time=00:10:00
#SBATCH --output=python_%j.out

ulimit -s unlimited
ulimit -l unlimited
umask 022

echo
echo "`date`"
time python convert_and_plot_projection.py
echo
echo "finished at: `date`"

