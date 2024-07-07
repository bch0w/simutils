#!/bin/sh                                                                        
                                                                                 
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
./bin/xcombine_vol_data_vtk \
    ${SLICE_LIST} ${FILENAME} ${INPUT_TOPO_DIR}/ ${INPUT_FILE_DIR}/ ${OUTPUT_DIR}/ ${RES}
echo
echo "finished at: `date`"
