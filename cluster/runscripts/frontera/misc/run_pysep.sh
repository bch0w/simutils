#!/bin/sh
#SBATCH -J example           # Job name
#SBATCH -o example_%j.out       # Name of stdout output file
#SBATCH -p development          # Queue (partition) name
#SBATCH -N 1                # Total # of nodes (must be 1 for serial)
#SBATCH -n 1                # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 00:01:00         # Run time (hh:mm:ss)


cd /home1/08540/bchow/scratch/

echo "LOADING SINGULARITY"
# module load tacc-singularity  

# The above doesnt work, throws the following error
# /opt/apps/tacc-singularity/3.7.2/etc/bash_completion.d/singularity: line 128: syntax error near unexpected token `<'
# /opt/apps/tacc-singularity/3.7.2/etc/bash_completion.d/singularity: line 128: `        done < <(compgen -W "${out[*]}" -- "$cur")'

echo "RUNNING SINGULARITY"
singularity run ../container-base_nightly.sif pysep -p mtuq_workshop_2022 -e 2020-04-04T015318_SOUTHERN_CALIFORNIA.yaml 
~
