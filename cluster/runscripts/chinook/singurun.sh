#!/bin/bash -e  

# This is a wrapper script to run container tasks through singularity.
# Arguments to be run are fed in through the command line

# Location of the Docker Image which should have been pulled manually
SIF_IMAGE="/import/c1/ERTHQUAK/bhchow/REPOSITORIES/containers/adjtomo_latest.sif"
CMD=${1}  # cmd line argument

# Flag Descriptions:
# -c/--contain: preserves internal container filesystem
# --bind: binds current working directory on system to a directory /home1 that
#	is inside the container. 
singularity exec -c --bind $(pwd):/home1 ${SIF_IMAGE} bash -c "${CMD}$
