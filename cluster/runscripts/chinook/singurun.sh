#!/bin/bash -e  

# This is a wrapper script to run container tasks through singularity.
# Arguments to be run are fed in through the command line

# Location of files and directories on local filesystem to bind into container
# SIF_IMAGE="/import/c1/ERTHQUAK/bhchow/REPOSITORIES/containers/adjtomo_centos7.sif"
SIF_IMAGE="/import/c1/ERTHQUAK/bhchow/REPOSITORIES/containers/pysep_centos7.sif"
HOME1="/import/c1/ERTHQUAK/${USER}"  # bind as /home1
WORK=$(pwd)

# Runs whatever the User provides as the command line argument after script
# Flag Descriptions:
# -c/--contain: preserves internal container filesystem
# --bind: binds current working directory on system to a directory /home1 that
#	is inside the container. 
singularity exec -c \
    --bind ${WORK}:/work --bind ${HOME1}:/home1 ${SIF_IMAGE} \
    bash -c "cd /work; $*"

