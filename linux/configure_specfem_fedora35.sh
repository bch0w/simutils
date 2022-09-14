# This was used for compiling SPECFEM3D on my Fedora35 Linux workstation
# with MPI enabled. I had to manually set where the MPI header and 
# Fortran compiler were as they weren't exposed on PATH
./configure MPI_INC=/usr/include/openmpi-x86_64 MPIFC=/usr/lib64/openmpi/bin/mpif90
