#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=24
#SBATCH --partition=debug
#SBATCH --time=00:10:00
#SBATCH --job-name="%j_pyatoa_windowing.log"

module purge
module load slurm

ulimit -l unlimited

eval "$(conda shell.bash hook)"

conda activate /import/c1/ERTHQUAK/bhchow/REPOS/miniconda3/envs/adjtomo
python ytian_pysep_synthetic_windowing.py
