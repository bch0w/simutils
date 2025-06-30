# Removes unncessary files from a SPECFEM3D repository cloned from GitHub to
# save storage space. Particularly used for reducing space requirements of 
# a Docker Image.

# Tested w/ commit ee3b095 on 05/14/2024
# this MUST be run inside a SPECFEM3D repository
# brings down overall directory size from 1GB to <70Mb

# MAIN DIRECTORY
rm -rf obj
rm -rf doc
rm -rf external_libs
rm -rf .git
rm -rf CUBIT_GEOCUBIT

# Large dataset in src/
rm -rf src/inverse_problem_for_source

# Remove some larger utils
rm -rf utils/dynamic_rupture
rm -rf utils/Visualization
rm -rf utils/unused_routines

# Remove a few very large examples
rm -rf EXAMPLES/benchmarks
rm -rf EXAMPLES/applications
rm -rf EXAMPLES/reproducible_study/
rm -rf EXAMPLES/CPML_examples
rm -rf EXAMPLES/Gmsh_*
rm -rf EXAMPLES/decompose_*
rm -rf EXAMPLES/fault_examples
rm -rf EXAMPLES/inversion_examples
rm -rf EXAMPLES/noise_non_uniform
rm -rf EXAMPLES/noise_tomography
rm -rf EXAMPLES/oldstuff
rm -rf EXAMPLES/small_*
rm -rf EXAMPLES/waterlayered_*

# Remove large files from meshfem3D_examples
rm -rf EXAMPLES/meshfem3D_examples/cavity
rm -rf EXAMPLES/meshfem3D_examples/simple_model
rm -rf EXAMPLES/meshfem3D_examples/*/REF_SEIS*

