#!/bin/bash
# Line 1: Velocity will be shown in movie
# 	1 = create files in OpenDX format
# 	2 = create files in AVS UCD format
# 	3 = create files in GMT xyz Ascii long/lat/Uz format
# 	    any other value = exit
# Line 2: Enter first time step of movie
# Line 3: Enter last time step of movie
# Line 4:  1= norm of velocity  2=velocity x-comp 3=velocity y-comp 4=velocity z-comp
# Line 5: 
#	1 = define file names using frame number
# 	2 = define file names using time step number

NSTEP=`grep ^NSTEP DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

mpiexec -n 16 ./bin/xcreate_movie_shakemap_AVS_DX_GMT << EOF
2
1
$NSTEP
2
4
EOF
