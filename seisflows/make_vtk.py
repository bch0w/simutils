"""
An automation script to take the output models from Seisflows, convert the .bin
files into .vtk files, collect into approriate directories, create differences
and plot standard visualizations using Mayavi
"""
import os
import sys
import time
import numpy as np
import subprocess
from glob import glob
from pprint import pprint
from os.path import join as oj
from pyatoa.utils.read import read_specfem_vtk
from pyatoa.utils.write import rcv_vtk_from_specfem, src_vtk_from_specfem


def setup_sum_dir(path_to_specfem):
    """
    Clean directory before progressing
    """
    sum_dir = oj(path_to_specfem, "SUM")    
    if not os.path.exists(sum_dir):
        os.makedirs(sum_dir)

    fids = glob(oj(sum_dir, "*"))
    if fids:
        check = input("files found in SUM directory, clean? [y/(n)]: ")
        if check == "y":
            for f in fids:
                os.remove(f)
        else:
            sys.exit("will not overwrite SUM directory, exiting")

def discover_available(paths):
    """
    Determine which gradients, models, kernels and parameters are available to 
    create vtk files for.
    """
    path_to_output = os.path.join(paths["cwd"], "output")
    assert(os.path.exists(path_to_output)), "output path does not exist"

    def split_tag(t):
        """split the output seisflows tag between the proc_ and .bin"""
        t = os.path.basename(t)
        t = "_".join(t.split("_")[1:])  # e.g. mu_kernel.bin
        return "".join(t.split(".")[:-1])
  
    print("discovering available models/gradients/kernels") 
    all_files = glob(oj(path_to_output, "*_????"))

    labels = set([os.path.basename(_).split("_")[0] for _ in all_files])
    available = {t:{} for t in labels}
    for l in labels:
        # kernel directory must be treated differently
        if l == "kernels":
            kernels = glob(oj(path_to_output, "kernels_*"))  # e.g. kernels_0001
            kernel_tags = [os.path.basename(_) for _ in kernels]

            # !!! Assuming here that all kernels have the same events/params
            events = glob(oj(kernels[0], "*"))
            event_tags = [os.path.basename(_) for _ in events]                
            # summed kernels are already represented by the gradient
            for s in ["sum", "sum_nosmooth"]:
                try:
                    event_tags.remove(s)
                except ValueError:
                    continue    

            # Determine the unique event and parameters used
            bin_files = glob(oj(events[0], "*"))
            params = set([split_tag(_) for _ in bin_files])
           
            available[l]["tags"] = kernel_tags
            available[l]["events"] = event_tags
            available[l]["params"] = list(params)

        else:
            directories = glob(oj(path_to_output, f"{l}_*"))
            tags = [os.path.basename(_) for _ in directories]

            bin_files = glob(oj(directories[0], "*"))
            params = set([split_tag(_) for _ in bin_files])

            # Check that these gradients/models have not been made yet
            # Get a list of existing vtk files in the output directory
            check = [os.path.basename(_).split('.')[0] for _ in 
                                     glob(os.path.join(paths['vtkout'], l,'*'))]

            # Here we assume that if one parameter exists, then all exit
            import ipdb;ipdb.set_trace()
            for t in tags:
                for p in params:
                    if f"{t}_{p}" in check:
                        print(f"\t{t}_{p} already exists in output directory")
                        tags.remove(t)
                        break

            available[l]["tags"] = tags
            available[l]["params"] = list(params)
   
    if not available:
        sys.exit("no available models/gradients/kernels found")
    else:
        for l in available.keys():
            print(f"{l}: {len(available[l]['tags'])}", end=" / ")
        print("")

    return available


def organize_models_gradients(pars, paths):
    """
    Models and gradients are organized by parameter and model number
    """
    c = 0
    for unique_id in pars["tags"]:  # e.g. model_0001
        directory = oj("output", unique_id)
        for par in pars["params"]:
            srcs = glob(oj(directory, f"proc*_{par}*.bin"))
            for src in srcs:
                # e.g. proc000000_vs.bin
                filename = os.path.basename(src)
                # Add a unique tag to the filename, different if kernels
                parts = filename.split("_")
                # Rejoin end of string if e.g. 'vs_kernel' and not 'vs'
                parts[1] = "_".join(parts[1:]) 
                new_filename = "_".join([parts[0], unique_id, parts[1]])
                # Set the source and destination for symlink
                dst = oj(paths["specfem"], "SUM", new_filename)
                os.symlink(os.path.abspath(src), dst)
                c += 1
    return c


def organize_kernels(kernel, pars, paths):
    """
    Kernels are organized by event id, so require a further level over
    gradient and model directories
    """
    c = 0
    # Number of counts relates to time required for xcombine_vol_data_vtk
    events = glob(oj("output", kernel, "*"))
    for event in events:
        event_id = os.path.basename(event)
        for par in pars["params"]:
            srcs = glob(oj(event, f"proc*_{par}*.bin"))
            for src in srcs:
                # e.g. proc000000_vs.bin
                filename = os.path.basename(src)
                # Add a unique tag to the filename, different if kernels
                parts = filename.split("_")
                # Rejoin end of string if e.g. 'vs_kernel' and not 'vs'
                parts[1] = "_".join(parts[1:])
                new_filename = "_".join([parts[0], kernel, event_id, 
                                         parts[1]])
                # Set the source and destination for symlink
                dst = oj(paths["specfem"], "SUM", new_filename)
                os.symlink(os.path.abspath(src), dst)
                c += 1
    return c
    

def run_xcombine_vol_data_vtk(paths):
    """
    Use subprocess to call sbatch of an existing script. A bit hacky but avoids
    having to generate an entire sbatch script in here.
    """
    os.chdir(paths["specfem"])
    stdout = subprocess.run(f"sbatch {paths['xcomb']}", shell=True, 
                            capture_output=True, text=True).stdout

    # Maui specific command to check job status, similar to Seisflows
    for part in stdout.strip().split(" "):
        try:
            # integer ensures no trailing/leading whitespaces
            job_id = str(int(part))
            break
        except ValueError:
            continue
    else:
        sys.exit("subprocess did not return a job number") 

    print(f"waiting for job number {job_id} to finish")
    time.sleep(20)
    while True: 
        time.sleep(10)
        stdout = subprocess.run(
                        f"sacct -nL -o jobid,state -j {job_id}",
                        shell=True, capture_output=True, text=True).stdout

        # Choose the last occurrence of job id to check status from, as 
        # earlier precursor jobs may be completed but main job is still running
        # If job is still launching, else will wait another round
        status = None
        for lines in stdout.strip().split("\n"):
            if lines.split()[0] == job_id:
                status = lines.split()[1]
        if status is None:
            print("waiting for status")
            continue

        if status == "RUNNING":
            print(".", end="")
            continue
        elif status == "PENDING":
            continue
        elif status == "COMPLETED":
            print(f"finished")
            return
        elif status == "FAILED":
            sys.exit(f"failed, exiting")
        else:
            print(f"unknown status {status}")
            continue


def post_organize(paths):
    """
    Remove symlinks and organize files in tagged directories, ensure that the
    directory containing symlinks is empty after organization
    """
    print("organizing output files")

    # Move files
    for src in glob(oj(paths["specfem"], "SUM", "*.vtk")):
        # e.g. model, kernel, gradient
        tag = os.path.basename(src).split("_")[0]
        dst_dir = oj(paths["vtkout"], tag)
        if not os.path.exists(dst_dir):
            os.makedirs(dst_dir)
        dst = oj(dst_dir, os.path.basename(src))

        os.rename(src, dst)

    # Remove symlinks
    for fid in glob(oj(paths["specfem"], "SUM", "*.bin")):
        if os.path.islink(fid):
            os.remove(fid)

    # Just make sure the directory is empty, otherwise something went wrong
    assert(len(glob(oj(paths["specfem"], "SUM", "*"))) == 0)


def create_srcrcv_vtk(paths):
    """
    One-time creation of source and receiver VTK files from the current 
    directory using the CMTSOLUTIONS and STATIONS files
    """
    src_vtk_from_specfem(path_to_data=os.path.join(paths["cwd"], "sfDATA"),
                         path_out=paths["vtkout"])

    rcv_vtk_from_specfem(path_to_data=os.path.join(paths["cwd"], "sfDATA"),
                         path_out=paths["vtkout"])


if __name__ == "__main__":
    # PATHS SET HERE
    paths = {"cwd": os.getcwd(),
             "specfem": oj(os.getcwd(), "scratch", "solver", "mainsolver"),
             "vtkout": oj(os.getcwd(), "scratch", "preprocess", "vtks"),
             "xcomb": ("/scale_wlg_persistent/filesets/project/nesi00263/" 
                       "PyPackages/simutils/seisflows/combvoldatavtk.sh")
             }
    assert os.path.exists(paths["specfem"])

    # Run the workflow
    setup_sum_dir(paths["specfem"])
    available = discover_available(paths)
    pprint(available)
    check = input("Manual remove some files? y/[n]: ")
    if check == "y":
        # File tags are located in `available`
        import ipdb; ipdb.set_trace()

    for tag in available.keys():
        pars = available[tag]
        if tag == "kernels":
            for kernel in pars["tags"]:
                count = organize_kernels(kernel, pars, paths)
                print(f"{kernel}, {len(pars['events'])} events, "
                      f"{count} files symlinked")
                if count:
                    run_xcombine_vol_data_vtk(paths)
                    post_organize(paths)
        else:
            count = organize_models_gradients(pars, paths)
            print(f"{len(pars['tags'])} {tag}s, "
                  f"{len(pars['params'])} parameters, "
                  f"{count} files symlinked")
            if count:
                run_xcombine_vol_data_vtk(paths)
                post_organize(paths)

    create_srcrcv_vtk(paths)




