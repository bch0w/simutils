# For LuNGS machines, configure SPECFEM3D_Cartesian (probably also 2D/3D_GLOBE, 
# but untested) for CentOS 7 (last access 3/4/24). Thanks Aakash for testing.

module load mpi/openmpi-x86_64
./configure FC=gfortran CC=gcc MPIFC=mpif90 --with-mpi
