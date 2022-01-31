#!/bin/bash -e
#SBATCH --job-name=xcube2sph_p1
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --clusters=maui
#SBATCH --account=gns03247
#SBATCH --partition=nesi_research
#SBATCH --time=00:30:00
#SBATCH --output=log_cube2sph_%j.out

## MODIFY the module load command
module load gcc/8.3.0   
module switch PrgEnv-cray/6.0.5 PrgEnv-gnu   
module load cray-hdf5/1.10.2.0  
module load cray-netcdf/4.6.1.3

# This script has been split into two parts because Python can't be run from the
# main compute cluster and so has to be run on its own from an ancillary node
# and then the workflow continued
#################################
echo "running example: `date`"
current_dir=`pwd`

# sets up directory structure in current example directory
echo
echo "   setting up example..."
echo
## MODIFY the specfem_dir
specfem_dir="/scale_wlg_persistent/filesets/project/gns03247/bchow/specfem/cube2sph_specfem3d_b7ed7a33"
################################
cube2sph_dir="cube2sph_utils"
# cleans output files
rm -rf DATABASES_MPI*
rm -rf MESH*
rm -rf OUTPUT_FILES_*
mkdir -p OUTPUT_FILES
rm -rf OUTPUT_FILES/*
mkdir -p MESH

mkdir -p bin
cd bin/
rm -f *
cp $specfem_dir/bin/* ./
cd ../

# stores setup
cp DATA/Par_file_initmesh DATA/Par_file
cp DATA/meshfem3D_files/Mesh_Par_file_init DATA/meshfem3D_files/Mesh_Par_file
cp DATA/meshfem3D_files/Mesh_Par_file OUTPUT_FILES/
cp DATA/Par_file OUTPUT_FILES/

# get the number of processors, ignoring comments in the Par_file
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `
mkdir -p $BASEMPIDIR

# This is a MPI simulation
echo
echo "  running mesher on $NPROC processors..."
echo
srun -n $NPROC ./bin/xmeshfem3D
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi
#mv OUTPUT_FILES/ OUTPUT_FILES_initmesh
rm -rf OUTPUT_FILES

###########transform cube into cubed sphere, then decompose & generate database with reference model###
echo 
echo "  produce a HEX27 mesh then transform cube into cubed sphere  "
echo

cp -R MESH/ MESH-default
rm -rf DATABASES_MPI

mkdir -p OUTPUT_FILES
# stores setup
cp DATA/Par_file_ref DATA/Par_file
cp DATA/Par_file OUTPUT_FILES/

# prepare for convert HEX8 mesh to HEX27 mesh
cp MESH-default/* ${cube2sph_dir}/temp
