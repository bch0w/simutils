# If OUTPUT_FILES are given in a SPECFEM master folder, check 
# If it is a finished run, move it to the storage folder, if it is not finished
# then don't do anything

STORAGE=/scale_wlg_nobackup/filesets/nobackup/nesi00263/bchow/storage
RUNFOLDER=`pwd -P`
OUTPUT_FILES=${RUNFOLDER}/OUTPUT_FILES

# CHECK FOR OUTPUT FOLDER
if [ -d ${OUTPUT_FILES} ]
then
	if [ -f ${OUTPUT_FILES}/output_solver.txt ]
	then
		if echo | grep -q "End of the simulation" ${OUTPUT_FILES}/output_solver.txt
		then
			# assumes CMTSOLUTION files are created with the same format always
			EVENT_ID="$(grep "event name:" ${OUTPUT_FILES}/CMTSOLUTION | cut -c18-)"
			if [ -d ${STORAGE}/${EVENT_ID} ]
			then
				echo OUTPUT_FILE EXISTS IN STORAGE, SORT IT OUT
				exit
			fi
			mv ${OUTPUT_FILES} ${STORAGE}/${EVENT_ID}
		fi
	fi
fi
