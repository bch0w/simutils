"""Generate state files with the correct file naming"""
import os
import sys
from copy import copy
from glob import glob

with open("./model_slices_template.pvsm") as f:
    lines_og = f.readlines()

for fid in glob("./vtk_files/MODEL_*_reg_1_vs?.vtk"):
    filename = os.path.basename(fid)
    fid_out = filename.replace(".vtk", ".pvsm")
    path_out = os.path.join(f"state_files", fid_out)
    if os.path.exists(path_out):
        continue
    modelname = filename.split(".")[0]
    parameter = modelname.split("_")[-1]

    print(fid, modelname, parameter)

    lines = copy(lines_og)
    for i, line in enumerate(lines[:]):
        lines[i] = line.format(PARAMETER=parameter, MODEL_NAME=modelname)

    print(fid_out)
    with open(path_out, "w") as f:
        f.writelines(lines)

