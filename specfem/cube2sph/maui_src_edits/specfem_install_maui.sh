# !!! This was written with SPECFEM3D commit: b7ed7a33
# !!! Which is pretty outdated so the overwrites of some of these f90 scripts
# !!! may not be compatible with current (2022) versions. Tianshi said it is
# !!! OKAY to comment out all the replacements
#!/usr/bin/bash
module load gcc/8.3.0
module switch PrgEnv-cray/6.0.5 PrgEnv-gnu

# MODIFY module load command and specfem_dir for the path to specfem package
specfem_dir="/home/chowbr/gns03247/project/bchow/specfem/specfem3d_cube2sph"

#################  OKAY TO COMMENT BELOW  #######################
# cp specfem_patch/decompose_mesh.F90 specfem_patch/part_decompose_mesh.F90 ${specfem_dir}/src/decompose_mesh
cp specfem_patch/create_slice.f90 specfem_patch/create_slice_loc.f90 ${specfem_dir}/src/decompose_mesh
#cp Makefile_slice ${specfem_dir}
#################  OKAY TO COMMENT ABOVE  #######################

cd ${specfem_dir}
# ./configure FC=ifort CC=icc MPIFC=mpif90 --with-mpi -enable-vectorization
# ./configure FC=ftn CC=cc CXX=CC MPIFC=ftn MPICC=cc MPICXX=CC --with-mpi --enable-openmp --enable-vectorization

# MAUI CONFIGURATION
./configure FC=ftn CC=cc CXX=CC MPIFC=ftn MPICC=cc MPICXX=CC --with-mpi --enable-vectorization

make realclean
make all > compile_log
#make bin/xcreate_slice_loc -f Makefile_slice >> compile_log
#make bin/xcreate_slice -f Makefile_slice >> compile_log
