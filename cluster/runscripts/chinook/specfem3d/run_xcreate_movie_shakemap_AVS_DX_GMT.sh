#!/bin/bash -e
# Make ParaView movie files if needed
echo "======== Running xcreate_move_shakemap_AVS_DX_GMT ========"
./bin/xcreate_movie_shakemap_AVS_DX_GMT <<- EOF
2
1
${NSTEP}
2
4
EOF
mv ${OUTPUT_FILES}/moviedata* ${RESULTS}
mv ${OUTPUT_FILES}/*inp ${RESULTS}
