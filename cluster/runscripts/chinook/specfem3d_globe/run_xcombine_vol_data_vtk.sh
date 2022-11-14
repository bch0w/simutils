#!/bin/sh                                                                        
                                                                                 
#SBATCH --job-name=combine_vol_data_vtk
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --partition=debug                                                        
#SBATCH --time=00:02:00                                                          
#SBATCH --output=combine_vol_data_vtk_%j.log
                                                                                 
                                                                                 
ulimit -s unlimited                                                              
ulimit -l unlimited    

# USAGE
# xcombine_vol_data_vtk \
# slice_list filename input_topo_dir input_file_dir output_dir high/low-resolution [region]

# PARAMETERS
SLICE_LIST="all"
FILENAME=$1  # or filename, use $1 to use command line argumnent
INPUT_TOPO_DIR="DATABASES_MPI/"
INPUT_FILE_DIR="DATABASES_MPI/"
OUTPUT_DIR="OUTPUT_FILES/"
RES=0  # 0 for low-res, 1 for hi-res (hi-res is LARGE)

mkdir -p ${OUTPUT_DIR}

echo "xcombine_vol_data_vtk for ${FILENAME}"
echo
echo "`date`"
time srun -n 1 ./bin/xcombine_vol_data_vtk \
    ${SLICE_LIST} ${FILENAME} ${INPUT_TOPO_DIR}/ ${INPUT_FILE_DIR}/ ${OUTPUT_DIR}/ ${RES}
echo
echo "finished at: `date`"
