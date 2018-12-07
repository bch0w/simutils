#!/bin/bash
# TO DO: turn off attenuation, change hessian kernel
# CREATE A NEW SPECFEM RUN BASED ON EVENT ID NUMBER
# SHOULD BE RUN INSIDE THE SPECFEM3DMASTER RUN FOLDER
# EXAMPLE CALL: ./build_adjoint 2018p130600
EVENT_ID=$1
if [ -z "$1" ]
then
	echo "EVENT ID REQUIRED"
	exit
fi

# FLAG FOR change_simulation_type.pl
# -a -- adjoint calculation
# -f -- forward calculation w/ save_forward=.false. (DEFAULT)
# -b -- run both simultaneously
# -F -- forward w/ save_forward=.true.
SIMTYPE="-b"

RUNFOLDER=`pwd -P`
TOMO=/scale_wlg_nobackup/filesets/nobackup/nesi00263/bchow
PRIMER=${TOMO}/primer
STORAGE=${TOMO}/storage
CMTSOLUTION=${PRIMER}/cmtsolution_files/${EVENT_ID}CMTSOLUTION
TEMPLATE=${PRIMER}/simutils/run_templates/forward_simulation.sh

# ECHO CHECK
echo
echo Event ID: ${EVENT_ID}
echo Run Folder: ${RUNFOLDER}
echo

echo CHANGING SIMULATION TYPE
${RUNFOLDER}/utils/change_simulation_type.pl ${SIMTYPE}
echo

RA="RUNADJOINT_${EVENT_ID}.sh"
echo CREATING ADJOINT RUN SCRIPT: ${RA}
cp ${PRIMER}/simutils/run_templates/adjoint_simulation.sh ${RUNFOLDER}/${RA}
SED1="sed -i '3s/.*/#SBATCH --job-name="${EVENT_ID}"_adj/' ${RUNFOLDER}/${RA}"
eval ${SED1}

echo
echo
echo MAKE SURE par_file PARAMETERS ARE SET APPROPRIATELY
echo
echo | grep "SIMULATION_TYPE" ${RUNFOLDER}/DATA/Par_file
echo | grep "NSTEP" ${RUNFOLDER}/DATA/Par_file
echo | grep "DT  " ${RUNFOLDER}/DATA/Par_file
echo | grep "ATTENUATION 	" ${RUNFOLDER}/DATA/Par_file
echo | grep "SAVE_SEISMOGRAMS_*" ${RUNFOLDER}/DATA/Par_file
echo | grep "APPROXIMATE_HESS_KL" ${RUNFOLDER}/DATA/Par_file
echo
echo
