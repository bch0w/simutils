#!/bin/bash
#SBATCH --job-name=test
#SBATCH --time=00:01:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=nesi00263                  ### change your projectID here
#SBATCH --clusters=maui_ancil                 ### therewith you can also submit from other systems
#SBATCH --partition=nesi_prepost
#SBATCH --export=None

module load Anaconda3/5.2.0-GCC-7.1.0
module list

sacct -nL -o jobid,state -j 193254
ssh w-mauivlab01.maui.nesi.org.nz
sacct -nL -o jobid,state -j 193254

python -c "import numpy"
echo "success"
