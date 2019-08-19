# Before using Easybuild, manual selection of cluster modules was required
# This selection of compiler options was deemed the fastest by NeSI engineers 
# for the executable xspecfem3d, as of 2017, however further testing in 2019 
# further improved execution time, and the use of Easybuild recipes and tailored
# modules made this manual selection of modules outdated
module switch PrgEnv-cray PrgEnv-gnu
module switch gcc/8.3.0 gcc/7.1.0 # gcc 7.1.0 required for skylake
module load craype-x86-skylake # skylake processor
module list
echo modules loaded GNU, gcc-7.1.0, craype-x86-skylake
