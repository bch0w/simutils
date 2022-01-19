# From GEOTOOLS/cluster_util/startup/README_specfem
# Should be run from inside a Specfem repository
module purge
module load slurm
module load compiler/GCC/5.4.0-2.26
module load mpi/OpenMPI/1.10.3-GCC-5.4.0-2.26

./configure FC=gfortran MPIFC=mpif90

