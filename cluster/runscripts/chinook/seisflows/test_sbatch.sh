#!/bin/sh

#SBATCH --job-name=submit_seisflows
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --partition=debug
#SBATCH --time=00:01:00
#SBATCH --output=seisflows_%j.out

WORK="/import/c1/ERTHQUAK/bhchow/scratch"
IMAGE="/import/c1/ERTHQUAK/bhchow/REPOSITORIES/containers/pyatoa_centos7.sif"

module load singularity
singularity exec -c --bind ${WORK},/home:/home1 ${IMAGE} seisflows -w /import/c1/ERTHQUAK/bhchow/scratch submit
