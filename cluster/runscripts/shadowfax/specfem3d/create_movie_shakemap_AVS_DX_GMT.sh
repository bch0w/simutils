#!/bin/sh                                                                        
                                                                                 
# PARAMETERS
# 5 = create files in VTU format with individual files
# 1 = first time step
# -1 = last time step
# 1 = Z  (1=Z, 2=N, 3=E)
NSTEP=`grep ^NSTEP DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`
echo "xcreate_movie_shakemap_AVS_DX_GMT"
echo
echo "`date`"

# Last Entry: 1= norm of velocity  2=velocity x-comp 3=velocity y-comp 4=velocity z-comp
PROMPT="2\n1\n${NSTEP}\n2\n4"
./bin/xcreate_movie_shakemap_AVS_DX_GMT  <<< $'2\n1\n30000\n2\n4'

echo
echo "finished at: `date`"
