#!/bin/sh                                                                        
                                                                                 
# PARAMETERS
# 5 = create files in VTU format with individual files
# 1 = first time step
# -1 = last time step
# 1 = Z  (1=Z, 2=N, 3=E)
echo "xcreate_movie_AVS_DX"
echo
echo "`date`"
./bin/xcreate_movie_AVS_DX  <<< $'5\n1\n-1\n1'
echo
echo "finished at: `date`"
