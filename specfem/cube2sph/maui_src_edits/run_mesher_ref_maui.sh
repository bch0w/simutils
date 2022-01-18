#!/bin/bash -e
#SBATCH --job-name=xcube2sph 
#SBATCH --nodes=1 
#SBATCH --ntasks=40 
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

#################################
echo "running example: `date`"
current_dir=`pwd`

# sets up directory structure in current example directory
echo
echo "   setting up example..."
echo
## MODIFY the specfem_dir
specfem_dir="/scale_wlg_nobackup/filesets/nobackup/gns03247/bchow/tomo/cube2sph/specfem3d_workdir"
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
cd ${cube2sph_dir}
python3 hex8tohex27.py
#bash change_name.sh
for file in temp/*_27; do
  name=`echo $file|awk -Ftemp/ '{print $2}' |awk -F_27 '{print $1}'` 
  cp temp/${name}_27 ${current_dir}/MESH-default/${name}
done
cd ${current_dir}

# get the number of processors, ignoring comments in the Par_file
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `
mkdir -p $BASEMPIDIR

# decomposes mesh using the pre-saved mesh files in MESH-default
echo
echo "  decomposing mesh..."
echo
./bin/xdecompose_mesh $NPROC ./MESH-default $BASEMPIDIR
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi
sbatch slurm_cube2sph_mesher.sh
echo "done"
echo `date`


