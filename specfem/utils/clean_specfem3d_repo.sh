# Removes unncessary files from a SPECFEM3D repository cloned from GitHub to
# save storage space. Particularly used for reducing space requirements of 
# a Docker Image.

# Tested w/ commit efd332e on 09/14/2022
# this MUST be run inside a SPECFEM3D repository
# brings down overall directory size to <250Mb

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
rm -rf EXAMPLES/oldstuff
rm -rf EXAMPLES/CPML_examples
rm -rf EXAMPLES/analytic_examples
rm -rf EXAMPLES/fault_examples
rm -rf EXAMPLES/meshfem3D_examples/sep_bathymetry
rm -rf EXAMPLES/homogeneous_poroelastic

# Remove large files from meshfem3D_examples
rm -rf EXAMPLES/meshfem3D_examples/cavity
rm -rf EXAMPLES/meshfem3D_examples/simple_model
rm -rf EXAMPLES/meshfem3D_examples/*/REF_SEIS*

