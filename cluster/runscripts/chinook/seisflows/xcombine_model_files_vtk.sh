#!/bin/sh 
                                                                                 
#SBATCH --job-name=combine_model_vtk
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --partition=debug
#SBATCH --time=00:05:00                                                          
#SBATCH --output=combine_all.out
                                                                                 
ulimit -s unlimited                                                              
ulimit -l unlimited    

# USAGE
# xcombine_vol_data_vtk \
# slice_list filename input_topo_dir input_file_dir output_dir high/low-resolution [region]

# xcombine_vol_data_vtk parameters
SLICE_LIST="all"
INPUT_TOPO_DIR="COMBINE_DATABASES/"
INPUT_FILE_DIR="COMBINE_DATABASES/"
RES=0  # 0 for low-res, 1 for hi-res (hi-res is LARGE)
REGION=1

# output parameters
PATH_INPUT="/import/c1/ERTHQUAK/bhchow/work/seisflows/NAKVERSION/output"
PATH_OUT="VISUALS"
MODEL_OUTPUT="VISUALS/MODELS/"
GRADIENT_OUTPUT="VISUALS/GRADIENTS/"
KERNEL_OUTPUT="VISUALS/KERNELS"

# FLAGS
MAKE_KERNELS=false
MAKE_GRADIENTS=true
MAKE_MODELS=false

echo "`date`"

# To remember the calling directory
CWD=`pwd`

# Do this for all models
if ${MAKE_MODELS}; then
    for MODEL in MODEL_04
    do
        MODEL=${PATH_INPUT}/${MODEL}
        # Remove all symlinked data in preparation for next
        cd ${CWD}/${INPUT_FILE_DIR}
        rm proc*_vs?.bin  # e.g., vsv(_kernel)
        rm proc*_vp?.bin  # e.g., vsv(_kernel)
        rm proc*_rho.bin
        rm proc*_eta.bin

        # Symlink in all parameters available in the model
        ln -s ${MODEL}/* .

        # Combine all parameters into vtk files
        cd ${CWD}
        for FILENAME in vpv vph vsv vsh eta rho
        do
            echo ${FILENAME}
            time srun -n 1 ./bin/xcombine_vol_data_vtk \
                ${SLICE_LIST} ${FILENAME} ${INPUT_TOPO_DIR}/ \
                ${INPUT_FILE_DIR}/ ${PATH_OUT}/ ${RES} ${REGION}
        done

        # Rename all vtk files and move them to the correct location
        cd ${CWD}/${PATH_OUT}
        MODNAME=${MODEL##*/}  # Trim everything before the last '/'
        rename .vtk _${MODNAME}.vtk *
        mv *.vtk ${CWD}/${MODEL_OUTPUT}
    done
fi

# Do this for all gradients
if ${MAKE_GRADIENTS}; then
    for GRADIENT in GRADIENT_05
    do
        GRADIENT=${PATH_INPUT}/${GRADIENT}
        # Remove all symlinked data in preparation for next
        cd ${CWD}/${INPUT_FILE_DIR}
        rm proc*_vs?_kernel.bin  
        rm proc*_vp?_kernel.bin 
        rm proc*_rho_kernel.bin
        rm proc*_eta_kernel.bin

        ln -s ${GRADIENT}/* .

        # Do this for all parameters
        cd ${CWD}
        for FILENAME in vpv_kernel vph_kernel vsv_kernel vsh_kernel eta_kernel rho_kernel
        do
            echo ${FILENAME}
            time srun -n 1 ./bin/xcombine_vol_data_vtk \
                ${SLICE_LIST} ${FILENAME} ${INPUT_TOPO_DIR}/ \
                ${INPUT_FILE_DIR}/ ${PATH_OUT}/ ${RES} ${REGION}
        done
        # Rename all vtk files and move them to the correct location
        cd ${CWD}/${PATH_OUT}
        GRADNAME=${GRADIENT##*/}
        rename .vtk _${GRADNAME}.vtk *
        mv *.vtk ${CWD}/${GRADIENT_OUTPUT}
    done
fi

# Do this for all models
if ${MAKE_KERNELS}; then
    for KERNEL in ${PATH_INPUT}/KERNEL_*
    do
        KERNEL=${PATH_INPUT}/${KERNEL}
        # Remove all symlinked data in preparation for next
        cd ${CWD}/${INPUT_FILE_DIR}
        rm proc*_vs?_kernel.bin  
        rm proc*_vp?_kernel.bin 
        rm proc*_rho_kernel.bin
        rm proc*_eta_kernel.bin

        # Symlink in all parameters available in the model
        ln -s ${KERNEL}/* .

        # Combine all parameters into vtk files
        cd ${CWD}
        for FILENAME in vpv_kernel vph_kernel vsv_kernel vsh_kernel eta_kernel rho_kernel
        do
            echo ${FILENAME}
            time srun -n 1 ./bin/xcombine_vol_data_vtk \
                ${SLICE_LIST} ${FILENAME} ${INPUT_TOPO_DIR}/ \
                ${INPUT_FILE_DIR}/ ${PATH_OUT}/ ${RES} ${REGION}
        done

        # Rename all vtk files and move them to the correct location
        cd ${CWD}/${PATH_OUT}
        KERNAME=${KERNEL##*/}  # Trim everything before the last '/'
        rename .vtk _${KERNAME}.vtk *
        mv *.vtk ${CWD}/${KERNEL_OUTPUT}
    done
fi

cd ${CWD}
echo "finished at: `date`"

