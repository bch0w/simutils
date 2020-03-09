#!/usr/bin/env python
"""
This is a finalization script to run a few processes after an inversion is done
It takes the kernels that are saved into the 'output' directory of Seisflows,
and symlinks/renames them into a directory in the Specfem3D run folder.
It then produces the appropriate SBATCH script to run which will generate 
combined VTK files for each of the kernels/models, and after creating them, will
place them into the pyatoa.io VTK directory for easy access
"""
import os
import sys
import time
import glob
import subprocess

# sys.argv sets the method
# method == "build" will set up the symlinks and create the runscript
# method == None removes all symlinks and clears everything
try:
    method = sys.argv[1]
except IndexError:
    method = None

# Set user parameters
kernels = ["vs_kernel"]
models = ["vs"]

# Set paths here
seisflows = "/path/to/seisflows/"
specfem = os.getcwd()

# Auto finish pathing
output_script = os.path.join(specfem, 'run_temp_xcombine_vol_data_vtk.sh')
sf_dir = os.path.basename(seisflows)
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
glob_gradients = os.path.join(seisflows, "output", "gradient_????")
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
                "#SBATCH --account=nesi00263\n"
                "#SBATCH --clusters=maui\n"
                "#SBATCH --partition=nesi_research\n"
                "#SBATCH --hint=nomultithread\n"
                "#SBATCH --output=sum_temp.log\n"
                "#SBATCH --time=0:05:00\n")
        
        f.write('echo "starting at `date`"\n')
        
        # Get the proc values, assuming its the same for model and gradients
        test_model = glob.glob(glob_models)[0]
        procs = glob.glob(os.path.join(test_model, "*vs.bin"))
        procs.sort()
        proc_a = int(os.path.basename(procs[0]).split('_')[0].split('proc')[1])
        proc_z = int(os.path.basename(procs[-1]).split('_')[0].split('proc')[1])
         
        
        # Loop through all the available gradients, created by Specfem
        print("GRADIENTS")
        for gradient_dir in glob.glob(glob_gradients):
            gradient_fid = os.path.basename(gradient_dir)
            print(f"\t {gradient_fid}")

            # Loop through each of the kernels requested
            for kernel in kernels:
                # Loop through each of the proc files
                for src in glob.glob(
                        os.path.join(gradient_dir, f'proc*_{kernel}*')):
                    # Set the symlink target to contain the event id 
                    proc_name = os.path.basename(src).split('.')[0]
                    dst_fid = f"{proc_name}_{gradient_fid}.bin"
                    
                    dst = os.path.join(input_sum, dst_fid)
                    os.symlink(src, dst)

                # Write the srun call
                kernel_id = f"{kernel}_{gradient_fid}"
                to_write = (f"srun -n 1 {combine_bin} {proc_a} {proc_z} " + 
                            f"{kernel_id} {sem_insum} {sem_outsum} 0\n")
                f.write(to_write)

        # Loop through all the available models
        print("MODELS") 
        for model_dir in glob.glob(glob_models):
            # Get the name of the model
            model_fid = os.path.basename(model_dir)
            print(f"\t {model_fid}")

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
        f.write(f"echo MOVED ALL .vtk FILES TO {pyatoa_vtk_dir}\n")
        f.write('echo "finished at `date`"\n')
