#!/bin/bash

#SBATCH --job-name=smooth_kernel
#SBATCH --nodes=4
#SBATCH --ntasks=144
#SBATCH --cpus-per-task=1
#SBATCH --account=nesi00263
#SBATCH --partition=nesi_research
#SBATCH --hint=nomultithread
#SBATCH --time 1:00:00

currentdir=`pwd`

# NOTE smoothing time for alpha_kernel with CPU, 144 tasks, 4 nodes
#	   was ~3000 sec or 50 min

# srun -n NPROC ./bin/xmooth_sem SIGMA_H SIGMA_V KERNEL_NAME INPUT_DIR OUTPUT_DIR
#srun -n 144 ./bin/xsmooth_sem 5000. 1000. alpha_kernel OUTPUT_SUM/ OUTPUT_SUM/ .false
#srun -n 144 ./bin/xsmooth_sem 5000. 1000. beta_kernel OUTPUT_SUM/ OUTPUT_SUM/ .false
srun -n 144 ./bin/xsmooth_sem 5000. 1000. rho_kernel OUTPUT_SUM/ OUTPUT_SUM/ .false

echo
echo "done"


