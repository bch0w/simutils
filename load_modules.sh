module list
module switch PrgEnv-cray PrgEnv-gnu
module switch gcc/7.3.0 gcc/7.1.0 # gcc 7.1.0 required for skylake
module load craype-x86-skylake # skylake processor
echo modules loaded
module list
