#!/bin/bash -e
#SBATCH --job-name=xcube2sph 
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

cd ${cube2sph_dir}
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


