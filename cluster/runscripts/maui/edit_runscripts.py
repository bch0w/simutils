"""
Small utility function to quickly change the number of nodes and tasks for all
of the run scripts in this directory
"""
import os
from glob import glob

account = "gns03247"
num_tasks = None
num_nodes = None

# num_tasks = int(input("Number of tasks/processors?: "))
# num_nodes = num_tasks // 40

for fid in glob("*.sh"):
    if "xcombine" in fid:
        continue
    with open(fid) as f:
        lines = f.readlines()
    
    # Edit lines in place so we can write back to file
    for l, line in enumerate(lines):
        if "--nodes=" in line and num_nodes is not None:
            lines[l] = f"#SBATCH --nodes={num_nodes}\n"
        elif "--ntasks=" in line and num_tasks is not None:
            lines[l] = f"#SBATCH --ntasks={num_tasks}\n"
        elif "--account=" in line and account is not None:
            lines[i] = f"#SBATCH --account={account}\n"

    with open(fid, "w") as f:
        for line in lines:
            f.write(line)
        
        

                
            
        
