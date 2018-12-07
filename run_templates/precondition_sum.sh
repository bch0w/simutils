#!/bin/bash

#SBATCH --job-name=precond_sum
#SBATCH --nodes=4
#SBATCH --ntasks=144
#SBATCH --cpus-per-task=1
#SBATCH --account=nesi00263
#SBATCH --partition=nesi_research
#SBATCH --hint=nomultithread
#SBATCH --time 00:01:00

currentdir=`pwd`
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

# xsum_preconditioned_kernels: sums transverse isotropic kernels while 
# preconditioning with an approximate Hessian outputted by the adjoint sim
# srun -n 144 ./bin/xsum_preconditioned_kernels
# DEFAULTS ARGS:
# input file: kernels_list.txt
#	lists all event kernel directories that should be summed together. does not
#	recognize the smooth appendix (i.e. will only look for alpha_kernel), so
#	all the kernels need to be renamed beforehand
# input directory: INPUT_KERNELS/
# output direcotry: OUTPUT_SUM/

echo "	summing kernels with preconditioning"
srun -n ${NPROC} ./bin/xsum_preconditioned_kernels

