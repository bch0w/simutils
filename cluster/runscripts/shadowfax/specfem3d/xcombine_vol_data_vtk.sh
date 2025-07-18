QUANTITY=${1}
RES=0

NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
NPROC_START=0
NPROC_END=`expr $NPROC - 1`

DIR_IN=`grep ^LOCAL_PATH DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
DIR_OUT="./OUTPUT_FILES"
mkdir -p ${DIR_OUT}

./bin/xcombine_vol_data_vtk \
		${NPROC_START} ${NPROC_END} ${QUANTITY} ${DIR_IN}/ ${DIR_OUT}/ ${RES}
