"""
Small utility function to quickly change SBATCH information in ALL slurm submit
scripts within a directory. New parameters can also be added and old parameters
deleted. Should be run inside the directory containing the SLURM bash scripts
"""
import os
from glob import glob


# Parameters to change in the SBATCH commands. Unknown parameters will be added
change_dict = {"account": "phone_account123",
               "partition": "debug",
               "cpus-per-task": 23,
               "new_par": "new_val",
               }
# Parameters to delete if encountered
delete_list = ["nodes"]
# Files to ignore when globbing through
ignore_list = ["xcombine_vol_data_vtk"]


for fid in glob("*.sh"):
    print(f"considering file: {fid}")
    # Ignore anything in the ignore list, parts of filenames okay
    for ignore in ignore_list:
        if ignore in fid:
            print(f"\tignoring file")
            continue

    # Grab all lines in the text file
    with open(fid) as f:
        lines = f.readlines()
  
    new_lines = []
    found_keys = []

    # Iterate on lines, only consider those which fit the SBATCH command
    for l, line in enumerate(lines):
        if "#SBATCH" not in line:
            new_lines.append(line)
            continue
        
        # ASSUMPTION: Slurm allows two formats for inputting short: 
        # '-p debug' or long: '--partition=debug'. I'm assuming here that 
        # all my scripts are only using the long form
        parts = line.split("--")[1]
        try:
            key, val = parts.strip().split("=")
        except ValueError as e:
            key, val = parts.strip().split(" ")

        # Delete any keys that the user requests removed
        if key in delete_list:
            print(f"\tfound '{key}' in delete list, removing")
            continue

        # Check if key should be changed from current value
        elif key in change_dict.keys():
            found_keys.append(key)
            new_val = change_dict[key]
            print(f"\tchanging '{key}' to '{new_val}'")
            line = lines[l].replace(val, str(new_val))

        new_lines.append(line)

    # Figure out where the last line '#SBATCH' is so we can append after it
    for l, line in enumerate(new_lines):
        if "#SBATCH" in line:
            last_line = l

    # Check if new keys should be added to the file
    i = 1
    for add_key in change_dict.keys():
        if add_key not in found_keys:
            add_val = change_dict[add_key]
            print(f"\tadding new pair: '{add_key}={add_val}'")
            new_line = f"#SBATCH --{add_key}={add_val}\n"
            new_lines.insert(last_line + i, new_line)
            i += 1

    # Finally overwrite the current file with new parameters
    with open(fid, "w") as f:
        f.writelines(new_lines)
        
        

                
            
        
