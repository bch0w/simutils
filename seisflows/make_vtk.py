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
from pyatoa.utils.read import read_specfem_vtk


# Globally accessible paths and parameters
cwd = os.getcwd()
path_to_specfem = os.path.join(cwd, "scratch", "solver", "mainsolver")
path_to_vtks = os.path.join(cwd, "pyatoa.io", "figures", "vtks")
path_to_xcomb = "/home/chowbr/primer/run_xcombine_vol_data_vtk.sh"

# Paths to auxiliary data
coast_fid = "/home/chowbr/primer/auxiliary/coastline/coast.npy"
srcs_fid = os.path.join(path_to_vtks, "srcs.vtk")
rcvs_fid = os.path.join(path_to_vtks, "rcvs.vtk")

# Parameter choices
parameters = ["vs", "vp"]
file_tags = ["model", "gradient"] 
model = "0003"
kernels = ["2*", "3*"]


def setup():
    """
    Clean directory before progressing
    """
    sum_dir = os.path.join(path_to_specfem, "SUM")    
    if not os.path.exists(sum_dir):
        os.makedirs(sum_dir)

    fids = glob(os.path.join(sum_dir, "*"))
    if fids:
        check = input("files found in SUM directory, clean? [y/(n)]: ")
        if check == "y":
            for f in fids:
                os.remove(f)
        else:
            sys.exit(0)


def pre_organize(file_tag="model"):
    """
    Bookkeeping function that symlinks all files with the unique tag names
    so that only one call to xcombine_vol_data_vtk is required

    :type file_tag: str
    :param file_tag: tag from Seisflows, ['model', 'gradient', 'kernel']
    :rtype: int
    :return: number of unique models that need to be turned into .vtk files
    """
    # Address kernels differently
    if file_tag == "kernels":
        count = _organize_kernels()
    else:
        count = _organize_models_gradients(file_tag)

    print(f"{count} {file_tag} unique ids symlinked")
    return count


def _organize_models_gradients(tag):
    """
    Models and gradients are organized by parameter and model number
    """
    count = 0
    # Number of counts relates to time required for xcombine_vol_data_vtk
    directories = glob(os.path.join("output", f"{tag}_{model}"))
    for directory in directories:
        unique_id = os.path.basename(directory)
        for par in parameters:
            count += 1
            srcs = glob(os.path.join(directory, f"proc*_{par}*.bin"))
            for src in srcs:
                # e.g. proc000000_vs.bin
                filename = os.path.basename(src)
                # Add a unique tag to the filename, different if kernels
                parts = filename.split("_")
                # Rejoin end of string if e.g. 'vs_kernel' and not 'vs'
                parts[1] = "_".join(parts[1:])
                new_filename = "_".join([parts[0], unique_id, parts[1]])
                # Set the source and destination for symlink
                dst = os.path.join(path_to_specfem, "SUM", new_filename)
                os.symlink(os.path.abspath(src), dst)

    return count

def _organize_kernels():
    """
    Kernels are organized by event id, so require a further level over
    gradient and model directories
    """
    # Number of counts relates to time required for xcombine_vol_data_vtk
    count = 0
    directories = glob(os.path.join("output", f"kernels_{model}"))
    for directory in directories:
        unique_id = os.path.basename(directory)
        for kernel in kernels:
            events = glob(os.path.join(directory, kernel))
            for event in events:
                event_id = os.path.basename(event)
                for par in parameters:
                    count += 1
                    srcs = glob(os.path.join(event, f"proc*_{par}*.bin"))
                    for src in srcs:
                        # e.g. proc000000_vs.bin
                        filename = os.path.basename(src)
                        # Add a unique tag to the filename, different if kernels
                        parts = filename.split("_")
                        # Rejoin end of string if e.g. 'vs_kernel' and not 'vs'
                        parts[1] = "_".join(parts[1:])
                        new_filename = "_".join([parts[0], unique_id, event_id, 
                                                 parts[1]])
                        # Set the source and destination for symlink
                        dst = os.path.join(path_to_specfem, "SUM", new_filename)
                        os.symlink(os.path.abspath(src), dst)
    
    return count


def run_xcombine_vol_data_vtk(count):
    """
    Use subprocess to call sbatch of an existing script. A bit hacky but avoids
    having to generate an entire sbatch script in here.

    :type count: int
    :param count: number of unique models that should be generated
    :return:
    """
    os.chdir(path_to_specfem)
    subprocess.call(f"sbatch {path_to_xcomb}", shell=True)
    # Wait for all files to be created by sbatch
    wait_for = "n"
    while wait_for != "y":
        wait_for = input("Finished running xcomb? y/[n]: ")

def post_organize():
    """
    Remove symlinks and organize files in tagged directories, ensure that the
    directory containing symlinks is empty after organization
    """
    # Move files
    for src in glob(os.path.join(path_to_specfem, "SUM", "*.vtk")):
        # e.g. model, kernel, gradient
        tag = os.path.basename(src).split("_")[0]
        dst_dir = os.path.join(path_to_vtks, tag)
        if not os.path.exists(dst_dir):
            os.makedirs(dst_dir)
        dst = os.path.join(dst_dir, os.path.basename(src))

        os.rename(src, dst)

    # Remove symlinks
    for fid in glob(os.path.join(path_to_specfem, "SUM", "*.bin")):
        if os.path.islink(fid):
            os.remove(fid)

    # Just make sure the directory is empty, otherwise something went wrong
    assert(len(glob(os.path.join(path_to_specfem, "SUM", "*"))) == 0)


def difference_models(log=True, poissons=True, outdir="model"):
    """
    Run diff_vtk() on all models, only if files don't already exist.
    Need to construct filenames before checking that they exist:
    Create either net model updates with log differences, or poissons ratios
    with poissons

    :type log: bool
    :param log: create net model update using log differences
    :type poissons: bool
    :param poissons: create poissons ratio vtk files
    :type outdir: str
    :param outdir: if given, name of the directory to save to the differenced 
        models to. If not given, saved into the model directory
    """
    outdir = os.path.join(path_to_vtks, outdir)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Create net model update, or log differences
    if log:
        for par in parameters:
            model_init = os.path.join(
                path_to_vtks, "model", f"model_init_{par}.vtk")
            models = glob(os.path.join(path_to_vtks, "model", f"model*{par}*"))
            for model in models:
                if "init" in model:
                    continue
                fid_out = os.path.join(outdir, f"log_{os.path.basename(model)}")
                if not os.path.exists(fid_out):
                    print(f"\t{os.path.basename(fid_out)}")
                    diff_vtk(model_a=model, model_b=model_init, fidout=fid_out,
                             method="log")

    # Calculate Poissons ratio
    if poissons:
        models_vp = glob(os.path.join(path_to_vtks, "model", "model*vp.vtk"))
        for model_vp in models_vp:
            model_vs = model_vp.replace("vp", "vs")
            fid_out = model_vp.replace("vp", "poissons")
            fid_out = os.path.join(outdir, os.path.basename(fid_out))
            if not os.path.exists(fid_out):
                print(f"\t{os.path.basename(fid_out)}")
                diff_vtk(model_a=model_vp, model_b=model_vs, fidout=fid_out,
                         method="poissons")


def diff_vtk(model_a, model_b, fidout, method="subtract"):
    """
    Read each model and scan line by line, difference all necessary values,
    write out a new vtk file that is the combination of A and B

    Different methods are:
    1) subtract: c = a - b
    2) log: c = log(a / b) ~ (a-b)/b, where a=m_i, b=m_0
    3) poisson: c = 0.5 * (a**2 - 2 * b**2) / (a**2 - b**2)
        where a=Vp, b=Vs

    :type model_a: str
    :param model_a: path to model file A
    :type model_b: str
    :param model_b: path to model file B
    :type method: str
    :param method: method for differencing,
        choice = ['subtract', 'log', 'poissons']
    :type fidout: str
    :param fidout: name of the output file to be written
    :rtype differences: list
    :return differences: list of the differences in values between a and b
    """
    # read files
    model_a, header_dict_a = read_specfem_vtk(model_a)
    model_b, header_dict_b = read_specfem_vtk(model_b)

    # check that the files have the same characteristics before parsing
    for key in header_dict_a.keys():
        if key == "scalars":
            continue
        elif header_dict_a[key] != header_dict_b[key]:
            sys.exit("{} not equal".format(key))

    # parse through models together and separate by len of line, skip header
    differences = []
    start = header_dict_a["data_line"]
    for a, b in zip(model_a[start:-1], model_b[start:-1]):
        try:
            a = float(a.strip())
            b = float(b.strip())
            if method == "subtract":
                difference = a - b
            # Take the natural log of the the quotient of a and b, this gives
            # to first order approximation, the percent difference. Yoshi said
            # Albert Tarantola said, "always view models in log space"
            elif method == "log":
                difference = np.log(a / b)
            elif method == "poissons":
                difference = 0.5 * (a**2 - 2 * b**2) / (a**2 - b**2)

            differences.append(difference)
        except ValueError:
            print("value error")

    # Write out the differences to a new file, newline at the end to play nice
    with open(fidout, "w") as f:
        f.writelines(model_a[:header_dict_a["data_line"]])
        for diff in differences:
            if diff == 0:
                f.write("{:13.5f}    \n".format(float(diff)))
            elif abs(diff) > 1:
                f.write("{:13.5f}    \n".format(float(diff)))
            else:
                f.write("{:13.10f}    \n".format(float(diff)))
        f.write("\n")

    return differences


if __name__ == "__main__":
    if not os.path.exists(path_to_specfem):
        sys.exit(-1)
    # Run functions in order
    c = 0
    setup()
    for ftag in file_tags:
        c_ = pre_organize(ftag)
        c += c_
    run_xcombine_vol_data_vtk(c)
    post_organize()
    difference_models(outdir="diffs")


