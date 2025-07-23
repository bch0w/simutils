#!/bin/bash
# Re-organize SPECFEM working directory after running `xspecfem3D` to keep
# output directory clean and to set up for subsequent runs so that critical 
# files are not overwritten

set -u  # exit if undefined variables are used
# set -e  # exit if failed command

TAG=${1}
echo ${TAG}
cd OUTPUT_FILES/

# Get rid of unneeded files
rm timestamp??????
rm starttimeloop.txt

# Move waveforms
mkdir -p ${TAG}/
mkdir ${TAG}/SEM
mv *.sem? ${TAG}/SEM

# Move Moviedata
mkdir ${TAG}/MOVDATA
mv moviedata?????? ${TAG}/MOVDATA

mkdir ${TAG}/AVSMOVIE
mv *.inp ${TAG}/AVSMOVIE

# Move log files
mv output_solver.txt ${TAG}
mv output_list_s*.txt ${TAG}
mv gpu_*txt ${TAG}

# Copy source files
cp ../DATA/CMTSOLUTION ${TAG}
cp ../DATA/STATIONS ${TAG}
cp ../DATA/Par_file ${TAG}

