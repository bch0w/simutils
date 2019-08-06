# module switch PrgEnv-cray PrgEnv-gnu
# module switch gcc/8.3.0 gcc/7.1.0 # gcc 7.1.0 required for skylake
module load gcc/7.1.0
module load craype-x86-skylake # skylake processor
echo modules loaded
