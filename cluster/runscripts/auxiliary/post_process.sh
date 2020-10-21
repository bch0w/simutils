#!/bin/bash

#SBATCH --job-name=postprocess_20??p??????
#SBATCH --nodes=4
#SBATCH --ntasks=144
#SBATCH --cpus-per-task=1
#SBATCH --account=nesi00263
#SBATCH --partition=nesi_research
#SBATCH --hint=nomultithread
#SBATCH --time 0:15:00

# This script is responsible for all post-processing steps following an
# adjoint simulation. Steps will be explained along the way in comments.

currentdir=`pwd`
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
SIGMA_H=5000.
SIGMA_V=1000.
SUMDIR=OUTPUT_SUM/

echo "smoothing kernels"
# xsmooth_sem: Smooths kernels by convolution with a gaussian. 
# note: smoothing time for alpha_kernel with CPU, 144 PROCs, 4 nodes was 50 min 
# srun -np NPROC ./bin/xsmooth_sem SIGMA_H SIGMA_V KERNEL_NAME INPUT_DIR OUTPUT_DIR USE_GPU
# ARGUMENTS:
#	SIGMA_H		- horizontal smoothing radius
#	SIGMA_V		- vertical smoothing radius
#	KERNEL_NAME	- kernel to be smoothed, e.g. alpha_kernel
#	INPUT_DIR	- directory where kernels are located
#	OUTPUT_DIR	- directory where smoothed kernels are written
#	USE_GPU		- use GPUs for computation
echo "	smoothing alpha kernel"
srun -n 144 ./bin/xsmooth_sem ${SIGMA_H} ${SIGMA_V} alpha_kernel ${SUMDIR} ${SUMDIR} .false
echo "	smoothing beta kernel"
srun -n 144 ./bin/xsmooth_sem ${SIGMA_H} ${SIGMA_V} beta_kernel ${SUMDIR} ${SUMDIR} .false
echo "	smoothing rho kernel"
srun -n 144 ./bin/xsmooth_sem ${SIGMA_H} ${SIGMA_V} rho_kernel ${SUMDIR} ${SUMDIR} .false
echo "smoothing approximate Hessian"
srun -n 144 ./bin/xsmooth_sem ${SIGMA_H} ${SIGMA_V} hess_kernel ${SUMDIR} ${SUMDIR} .false

# ==============================================================================
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
echo "	renaming smoothed kernels for use in summation"
rename kernel.bin kernel_original.bin ${SUMDIR}/*kernel.bin
rename kernel_smooth.bin kernel.bin ${SUMDIR}/*kernel_smooth.bin

echo "	summing kernels with preconditioning"
srun -n 144 ./bin/xsum_preconditioned_kernels

# ==============================================================================
# xadd_model_iso can be used to update isotropic model files with smoothed and 
# summed kernels. the algorithm uses a steepest descent method with step length
# determined by the given maximum update percentage

# command line call:
# run -np 4 bin/xadd_model_iso step_factor [INPUT_KERNELS_DIR/] [OUTPUT-MODEL-DIR/]
# 	step_factor: step length to scale the gradient, e.g. 0.03 for a +-3% update
#	INPUT-KERNELS-DIR/: optional directory that holds summed kernels (e.g. proc****alpha_kernel.bin), defaults to INPUT_GRADIENT/
#	OUTPUT-MODEL-DIR/: optional directory that will hold new model files, defaults to OUTPUT_MODEL/

srun -n 144 ./bin/xadd_model_iso 0.03 INPUT_GRADIENT/ OUTPUT_MODEL/ 

echo
echo "done"

# Auxiliary function is summing without preconditioning is preferred. Usually no

# echo "summing sensitivity kernels"
# # srun -np NPROC bin/xcombine_sem KERNEL_NAMES INPUT_FILE OUTPUT_DIR
# # ARGUMENTS: 
# # 	KERNEL_NAMES - one or more kernel names separated by commas, 
# #				   e.g. alpha_kernel,beta_kernel,rho_kernel
# #	INPUT_FILE   - text file containing list of kernel directories
# #	OUTPUT_PATH  - directory where summed kernels is written
# #
# srun -n $NPROC ./bin/xcombine_sem alpha_kernel,beta_kernel_rho_kernel kernels_list.txt OUTPUT_SUM/
# 
