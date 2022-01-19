# From GEOTOOLS/cluster_util/startup/README_specfem
# Should be run from inside a Specfem repository
module purge
module load slurm
module load compiler/ifort/2018.5.274-GCC-5.4.0-2.26  
module load openmpi/intel/2.1

./configure FC=ifort MPIFC=mpif90

