"""
Tired of having to open the SPECFEM3D Par_file every time you want to edit 
a silly parameter? This Python script is designed to allow a User top quickly
access, assess, and edit a parameter in the SPECFEM3D Par_file
"""
import os
import sys


# Determine if were accessing or editing
try:
    name = sys.argv[1]
except IndexError:
    sys.exit(f"\n\tfirst argument must be the parameter name "
             f"(case insensitive)\n")
try:
    newval = sys.argv[2]
except IndexError:
    newval = None

# Figure out where the Par_file is, typically should be run in a SPECFEM workdir
cwd = os.getcwd()
cases = ["Par_file", "DATA/Par_file", "../DATA/Par_file"]
for case in cases:
    path = os.path.relpath(os.path.join(cwd, case))
    if os.path.exists(case):
        fid = case
        break
else:
    import pdb;pdb.set_trace()
    sys.exit(f"Cannot find Par_file, please make sure you are in a work dir")

with open(fid, "r") as f:
    lines = f.readlines()

for i, line in enumerate(lines[:]):
    # Skip comments and newlines
    if not line.strip() or line[0] == "#":
        continue
    try:
        key, val, *_ = [_.strip() for _ in line.strip().split("=")]
        # If comments after the line, remove
        if "#" in val:
            val = val.split("#")[0].strip()
    except ValueError:
        print(f"\n\tUnexpected line: {line}\t")
        continue
    if name.upper() == key:
        # If accessing, print and move on
        if newval is None:
            sys.exit(f"\n\t{line}")
        else:
            lines[i] = line.replace(val, newval)
            print(f"\n\t{key}: {val} -> {newval}\n")
            break
else:
    sys.exit(f"\n\tparameter '{name}' not found\n")

with open(fid, "w") as f:
    f.writelines(lines)

