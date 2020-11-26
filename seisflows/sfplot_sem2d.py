"""
A simple way to plot the output SPECFEM2D binary files that are recovered
during a SeisFlows inversion. 
"""
import os
import numpy as np
from glob import glob

import matplotlib.pyplot as plt
from pyatoa.utils.read import read_fortran_binary

assert(os.path.exists("parameters.yaml")), "Parameter file not found"

# Figure out the number of processors based on the parameter file
with open("parameters.yaml") as f:
    lines = f.readlines()
    for line in lines:
        if "NPROC: " in line:
            nproc = int(line.strip().split(":")[-1])

# Read in the spatial coordinates once
xvals = np.array([])
zvals = np.array([])
for i in range(0, nproc, 1):
    # Read in X values for each processor
    xvals_ = read_fortran_binary(f"output/model_init/proc{i:0>6}_x.bin")
    xvals = np.concatenate((xvals, xvals_))

    # Read in Z values for each processor
    zvals_ = read_fortran_binary(f"output/model_init/proc{i:0>6}_z.bin")
    zvals = np.concatenate((zvals, zvals_))

# Read in model values for each parameter
parameters = ["vp", "vs", "rho"]    
for model_path in glob(f"output/model_*"): 
    model_name = os.path.basename(model_path)
    for par in parameters:
        # Check if the parameter has any .bin files
        if not glob(f"{model_path}/proc*_{par}.bin"):
            print(f"{model_name}_{par} not found, skipping...")
            continue

        # Read in parameter values for each processor
        vals = np.array([])
        for i in range(0, nproc, 1):
            vals_ = read_fortran_binary(f"{model_path}/proc{i:0>6}_{par}.bin")
            vals = np.concatenate((vals, vals_))

        # Plot and save the cwd
        plt.scatter(xvals, zvals, c=vals)
        plt.colorbar()
        plt.title(f"{model_name} {par}")
        plt.xlabel("X")
        plt.ylabel("Z")

        fid_out = f"{model_name}_{par}.png"
        print(f"saving figure {fid_out}")
        plt.savefig(fid_out)
        plt.close()

