#!/bin/bash -e

#SBATCH --job-name=xcube2sph
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1
#SBATCH --clusters=maui
#SBATCH --account=gns03247
#SBATCH --partition=nesi_research
#SBATCH --time=00:30:00
#SBATCH --output=log_cube2sphp3_%j.out


# Machine specific module loadout
module load gcc/8.3.0  
module switch PrgEnv-cray/6.0.5 PrgEnv-gnu  
module load cray-hdf5/1.10.2.0 
module load cray-netcdf/4.6.1.3


# Set executable here
exc="srun -n" 

# Set directories here
current_dir=`pwd`
specfem_dir="/scale_wlg_nobackup/filesets/nobackup/gns03247/bchow/tomo/cube2sph/cube2sph_specfem3d_b7ed7a33"
# cube2sph_dir="cube2sph_utils"
cube2sph_dir="/scale_wlg_persistent/filesets/project/gns03247/bchow/repos/cube2sph_update/cube2sph_utils"

# Typical Specfem run script dir grabs
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

# ================================ RUN CUBE2SPH ================================
echo
echo "  cube2sph transform on $NPROC processors..."
echo
## MODIFY cube2sph transform parameters - mesh center and rotation
# cube2sph CENTER_LAT CENTER_LON ROTATE_ANGLE
# !!! ${exc} $NPROC ${cube2sph_dir}/bin/cube2sph 175.75 -57.75 20.0


# ================================ RUN XGENDBS =================================
echo
echo "  running database generation on $NPROC processors..."
echo
# !!! ${exc} $NPROC ./bin/xgenerate_databases


# ================================ RUN XCOMBVOL ================================
echo
echo "    generate vtk file for vs... "
echo
# !!! ./bin/xcombine_vol_data_vtk 0 $((NPROC-1)) vs DATABASES_MPI/ . 0
# mv vs.vtk vs_ref.vtk

mkdir -p ${cube2sph_dir}/DATABASES_MPI_REF
mkdir -p ${cube2sph_dir}/DATABASES_MPI

mv DATABASES_MPI/* ${cube2sph_dir}/DATABASES_MPI_REF
rm -rf OUTPUT_FILES
########change the model, then generate database again################

mkdir -p OUTPUT_FILES
# stores setup
cp DATA/Par_file_gll DATA/Par_file
cp DATA/Par_file ${cube2sph_dir}/DATA

# get the number of processors, ignoring comments in the Par_file
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

#BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `
#mkdir -p $BASEMPIDIR

cd ${cube2sph_dir}

echo
echo "    adding surface and interior topography..."
echo
#./bin/node_stretching
${exc} $NPROC ./bin/node_stretching_parallel
#cp nodes_coords_file_topo ../../MESH-default
cp DATABASES_MPI/*Database ${current_dir}/DATABASES_MPI

#mkdir -p DATABASES_MPI
echo
echo "    changing the velocity model on GLL points..."
echo 

${exc} $NPROC ./bin/setup_model_cartesian
cp DATABASES_MPI/*vp.bin DATABASES_MPI/*vs.bin DATABASES_MPI/*rho.bin ${current_dir}/DATABASES_MPI

#echo
#echo "    writing force and station files in cartesian style..."
#echo

#./bin/write_force_solution_file
#./bin/write_stations_file

#cp DATA/FORCESOLUTION_cartesian ../../DATA/FORCESOLUTION
#cp DATA/STATIONS_cartesian ../../DATA/STATIONS

cd ${current_dir}
# backup data files
cp DATA/meshfem3D_files/Mesh_Par_file OUTPUT_FILES/
cp DATA/Par_file OUTPUT_FILES/
#cp DATA/FORCESOLUTION OUTPUT_FILES/
#cp DATA/STATIONS OUTPUT_FILES/

#cp MESH-default/nodes_coords_file_topo MESH-default/nodes_coords_file
# decomposes mesh using the pre-saved mesh files in MESH-default
#echo
#echo "    decomposing mesh with topography..."
#echo
#./bin/xdecompose_mesh $NPROC ./MESH-default $BASEMPIDIR

echo
echo "  running database generation on $NPROC processors..."
echo
${exc} $NPROC ./bin/xgenerate_databases

echo
echo "    generate vtk file for vs... "
echo
./bin/xcombine_vol_data_vtk 0 $((NPROC-1)) vs DATABASES_MPI/ . 0

######## forward simulation for one source

#cp DATA/FORCESOLUTION_globe cube2sph_utils/stretch/DATA/FORCESOLUTION
#cp DATA/STATIONS_globe cube2sph_utils/stretch/DATA/STATIONS
#
#cd cube2sph_utils/stretch
#
#echo
#echo "    writing force and station files in cartesian style..."
#echo
#
#./bin/write_force_solution_file
#./bin/write_stations_file
#
#cp DATA/FORCESOLUTION_cartesian ../../DATA/FORCESOLUTION
#cp DATA/STATIONS_cartesian ../../DATA/STATIONS
#
#cd ../..
#cp DATA/FORCESOLUTION OUTPUT_FILES/
#cp DATA/STATIONS OUTPUT_FILES/
#
#echo
#echo "    start solver..."
#echo
#sbatch --wait go_solver.sh
#if [[ $? -ne 0 ]]; then exit 1; fi
#echo
#echo "see results in directory: OUTPUT_FILES/"
#echo
echo "done"
echo `date`


