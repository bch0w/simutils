"""
Small utility function to quickly change the number of nodes and tasks for all
of the run scripts in this directory
"""
import os
from glob import glob

num_tasks = int(input("Number of tasks/processors?: "))
num_nodes = num_tasks // 40

for fid in glob("*.sh"):
    with open(fid) as f:
        lines = f.readlines()
    
    # Edit lines in place so we can write back to file
    for l, line in enumerate(lines):
        if "--nodes=" in line:
            lines[l] = f"#SBATCH --nodes={num_nodes}\n"
        elif "--ntasks=" in line:
            lines[l] = f"#SBATCH --ntasks={num_tasks}\n"

    with open(fid, "w") as f:
        for line in lines:
            f.write(line)
        
        

                
            
        