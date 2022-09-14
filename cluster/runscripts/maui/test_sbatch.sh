#!/bin/bash -e

#SBATCH --job-name=test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=%j.out
#SBATCH --clusters=maui
#SBATCH --array=0-1
#SBATCH --account=gns03247
#SBATCH --partition=nesi_research
#SBATCH --parsable

echo "hello"

