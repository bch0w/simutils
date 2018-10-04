#module list
module load slurm
module swap gcc/4.9.3 gcc/7.3.0 # gcc 7 required for skylake
module load craype-x86-skylake # skylake processor
echo modules loaded
#module list
