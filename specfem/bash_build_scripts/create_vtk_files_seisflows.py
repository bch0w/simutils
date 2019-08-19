#!/usr/bin/env python
"""
This is a finalization script to run a few processes after an inversion is done
"""
import os
import sys
import time
import glob
import subprocess

# sys.argv sets the method
try:
    method = sys.argv[1]
except IndexError:
    method = None

# Set user parameters
sf_dir = "twoevent_6895e2f7"
kernels = ["vs_kernel"]
models = ["vs"]

# Can set these paths if they change
basepath = "/scale_wlg_nobackup/filesets/nobackup/nesi00263/bchow/tomo/"
seisflows = os.path.join(basepath, "seisflows", "checkerboard", sf_dir)
# specfem = os.path.join(basepath, "specfem_87a78d_CrayCCE-19.04")
specfem = os.path.join(basepath, "specfem_6895e2f7_GNU-7.1.0")


output_script = os.path.join(specfem, 'run_temp_xcombine_vol_data_vtk.sh')

# Auto finish pathing
input_sum_dir = "SUM"
output_sum_dir = f"SUM/{sf_dir}"
input_sum = os.path.join(specfem, input_sum_dir)
output_sum = os.path.join(specfem, output_sum_dir)
if not os.path.exists(output_sum):
    os.makedirs(output_sum)

# Specfem is a little stupid and requires relative pathing
sem_insum = f"./{input_sum_dir}/"
sem_outsum = f"./{output_sum_dir}/"
combine_bin = "./bin/xcombine_vol_data_vtk"

# Set paths to kernels and models
glob_kernels = os.path.join(seisflows, "output", "kernels_????")
glob_models = os.path.join(seisflows, "output", "model_????")
pyatoa_vtk_dir = os.path.join(seisflows, "pyatoa.io", "figures", "vtks")

# Remove any old symlinks left in the INPUT_SUM folder
for fid in glob.glob(os.path.join(input_sum, '*.bin')):
    os.unlink(fid)
if os.path.exists(output_script):
    os.remove(output_script)

# Build the submit script and symlink necessary files
if method == "build":
    with open(output_script, 'w') as f:
        f.write("#!/bin/bash\n"
                "#SBATCH --job-name=combine_vol_data_vtk\n"
                "#SBATCH --nodes=1\n"
                "#SBATCH --ntasks=1\n"
                "#SBATCH --cpus-per-task=1\n"
                "#SBATCH --account=gns02750\n"
                "#SBATCH --clusters=maui\n"
                "#SBATCH --partition=nesi_research\n"
                "#SBATCH --hint=nomultithread\n"
                "#SBATCH --output=sum_temp.log\n"
                "#SBATCH --time=0:05:00\n")
        
        f.write('echo "starting at `date`"\n')
        
        # Get the proc values, assuming its the same for model and kernels
        test_model = glob.glob(glob_models)[0]
        procs = glob.glob(os.path.join(test_model, "*vs.bin"))
        procs.sort()
        proc_a = int(os.path.basename(procs[0]).split('_')[0].split('proc')[1])
        proc_z = int(os.path.basename(procs[-1]).split('_')[0].split('proc')[1])
         
        
        # Loop through all the available kernels, created by Specfem
        for kernel_model in glob.glob(glob_kernels):
            model_fid = os.path.basename(kernel_model)
            for kernel_dir in glob.glob(os.path.join(kernel_model, "*")):
                # Get the name of the kernel
                kernel_fid = os.path.basename(kernel_dir)
                print(f"kernel {model_fid} {kernel_fid}")

                # Loop through each of the kernels requested
                for kernel in kernels:
                    # Loop through each of the proc files
                    for src in glob.glob(
                            os.path.join(kernel_dir, f'proc*_{kernel}*')):
                        # Set the symlink target to contain the event id 
                        proc_name = os.path.basename(src).split('.')[0]
                        dst_fid = f"{proc_name}_{model_fid}_{kernel_fid}.bin"
                        
                        dst = os.path.join(input_sum, dst_fid)
                        os.symlink(src, dst)

                    # Write the srun call
                    kernel_id = f"{kernel}_{model_fid}_{kernel_fid}"
                    to_write = (f"srun -n 1 {combine_bin} {proc_a} {proc_z} " + 
                                f"{kernel_id} {sem_insum} {sem_outsum} 0\n")
                    f.write(to_write)

        # Loop through all the available models
        for model_dir in glob.glob(glob_models):
            # Get the name of the model
            model_fid = os.path.basename(model_dir)
            print(f"model {model_fid}")

            # Loop through each of the models requested
            for model in models:
                # Loop through each of the proc files
                for src in glob.glob(
                        os.path.join(model_dir, f'proc*_{model}.bin')):
                    # Set the symlink target to contain the event id 
                    proc_name = os.path.basename(src).split('.')[0]
                    dst_fid = f"{proc_name}_{model_fid}.bin"
                    
                    dst = os.path.join(input_sum, dst_fid)
                    os.symlink(src, dst)

                # Write the srun call
                model_id = f"{model}_{model_fid}"
                to_write = (f"srun -n 1 {combine_bin} {proc_a} {proc_z} " + 
                            f"{model_id} {sem_insum} {sem_outsum} 0\n")
                f.write(to_write)

        # Bash command to move the written files to pyatoa.io
        f.write(f"mv {output_sum}/* {pyatoa_vtk_dir}/\n")

        f.write('echo "finished at `date`"\n')
