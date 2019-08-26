# Original configuration optimized for GNU-7.1.0
echo './configure FC=ftn CC=cc CXX=CC MPIFC=ftn MPICC=cc MPICXX=CC --with-mpi'

# MPI/OpenMP Hybrid w/ vectorization
echo './configure FC=ftn CC=cc CXX=CC MPIFC=ftn MPICC=cc MPICXX=CC --with-mpi --enable-openmp --enable-vectorization'
