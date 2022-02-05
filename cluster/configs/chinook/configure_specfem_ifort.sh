# From GEOTOOLS/cluster_util/startup/README_specfem
# Should be run from inside a Specfem repository
module purge
module load slurm
# module load compiler/ifort/2018.5.274-GCC-5.4.0-2.26  
# module load openmpi/intel/2.1
module load data/HDF5/1.10.6-pic-intel-2019b
module load data/netCDF/4.7.4-pic-intel-2019b
module load data/netCDF-Fortran/4.5.2-pic-intel-2019b

./configure FC=ifort CC=icc MPIFC=mpif90 --with-mpi --enable-vectorization

