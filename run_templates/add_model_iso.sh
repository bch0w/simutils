#!/bin/bash

#SBATCH --job-name=add_model_iso
#SBATCH --nodes=4
#SBATCH --ntasks=144
#SBATCH --cpus-per-task=1
#SBATCH --account=nesi00263
#SBATCH --partition=nesi_research
#SBATCH --hint=nomultithread
#SBATCH --time 0:00:15

# xadd_model_iso can be used to update isotropic model files with smoothed and 
# summed kernels. the algorithm uses a steepest descent method with step length
# determined by the given maximum update percentage

# command line call:
# run -np 4 bin/xadd_model_iso step_factor [INPUT_KERNELS_DIR/] [OUTPUT-MODEL-DIR/]
# 	step_factor: step length to scale the gradient, e.g. 0.03 for a +-3% update
#	INPUT-KERNELS-DIR/: optional directory that holds summed kernels (e.g. proc****alpha_kernel.bin), defaults to INPUT_GRADIENT/
#	OUTPUT-MODEL-DIR/: optional directory that will hold new model files, defaults to OUTPUT_MODEL/

srun -n 144 ./bin/xadd_model_iso 0.02

echo
echo "done"


