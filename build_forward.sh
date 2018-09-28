#!/bin/bash
# CREATE A NEW SPECFEM RUN BASED ON EVENT ID NUMBER
# SHOULD BE RUN INSIDE THE SPECFEM3DMASTER RUN FOLDER
EVENT_ID=$1
DECOMPOSE=$2

RUNFOLDER=`pwd -P`
TOMO=/scale_wlg_nobackup/filesets/nobackup/nesi00263/bchow
PRIMER=${TOMO}/primer
STORAGE=${TOMO}/storage
CMTSOLUTION=${PRIMER}/cmtsolution_files/${EVENT_ID}CMTSOLUTION
TEMPLATE=${PRIMER}/simutils/run_templates/forward_simulation.sh

# CHECK AND EXITS:
# CHECK IF CMTSOLUTION FILE EXISTS
if ! [ -f ${CMTSOLUTION} ]
then
	echo ${CMTSOLUTION} DOES NOT EXIST
	exit
fi
# CHECK IF RUNFOLDER ALREADY EXISTS
if [ -d ${STORAGE}/${EVENT_ID} ]
then
	echo OUTPUT_FILES ALREADY EXISTS IN STORAGE
	exit
fi

# CHECK IF OUTPUT_FOLDER EXISTS IN RUNFOLDER
if [ -d ${RUNFOLDER}/OUTPUT_FILES ]
then
	echo OUTPUT_FILES ALREADY EXISTS IN RUN FOLDER, MOVING
	source ${PRIMER}/simutils/output_to_storage.sh
	echo
fi

# IF PASS CHECK-STOPS, CREATE AND RUN
if ! [ -d ${RUNFOLDER}/DATA/tomo_files ]
then
	echo tomo_files IS NOT PRESENT IN DATA
	exit
fi
if ! [ -f ${RUNFOLDER}/DATA/STATIONS ]
then
	echo STATIONS IS NOT PRESENT IN DATA
	exit
fi
if ! [ -d ${RUNFOLDER}/MESH/ ]
then
	echo MESH IS NOT PRESENT
	exit
fi

echo
echo | grep "SIMULATION_TYPE" ${RUNFOLDER}/DATA/Par_file
echo | grep "NSTEP" ${RUNFOLDER}/DATA/Par_file
echo | grep "DT  " ${RUNFOLDER}/DATA/Par_file
echo | grep "ATTENUATION" ${RUNFOLDER}/DATA/Par_file
echo | grep "SAVE_SEISMOGRAMS_*" ${RUNFOLDER}/DATA/Par_file
echo

echo SYMLINKING CMTSOLUTION
rm ${RUNFOLDER}/DATA/CMTSOLUTION
ln -s ${CMTSOLUTION} ${RUNFOLDER}/DATA/CMTSOLUTION

echo CREATING FORWARD RUN SCRIPT
cp ${PRIMER}/simutils/run_templates/forward_simulation.sh ${RUNFOLDER}/RUNFORWARD.sh
SED1="sed -i '3s/.*/#SBATCH --job-name="${EVENT_ID}"/' ${RUNFOLDER}/RUNFORWARD.sh"
SED2="sed -i '17s/.*/DECOMPOSE="${DECOMPOSE}"/' ${RUNFOLDER}/RUNFORWARD.sh"
SED3="sed -i '3s/.*/#SBATCH --job-name="${EVENT_ID}"/' ${RUNFOLDER}/RUNFORWARD.sh"
eval ${SED1}
eval ${SED2}
eval ${SED3} 

#sbatch ${RUNFOLDER}/RUNFORWARD.sh
