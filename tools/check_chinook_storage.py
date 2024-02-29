"""
Python wrapper for Bash commands to check storage on UAF HPC Chinook by active
User. Intended to be run by a weekly Cronjob with email output to a Users
group. 

Contains built in leeway for fail halfway through and rerunning the process w/o
modification

.. note::

    Rule of thumb 'du -sh' takes 45m/TB (based on ~2TB taking 1h30m)

.. note:: du flags
    
    c: print grand total
    h: human redable
    s: summarize (max depth = 0)
"""
import os
import sys
import subprocess
import time


# - 'test': has a few random Users that do not have much data
# - 'active': active users on ERTHQUAK (last check 2-14-24)
# - 'all': all Users in the ERTHQUAK group
USERS = {"test": ["ykaneko1_ua", "ecasarotti", "gthompson"],
         "active": ["bhchow", "ctape", "ammcpherson",  "jthurin", "agupta7", 
                    "eaonyango", "ytian4"],
         "all": None}
CHOICE = "active"  # select one key in Users
CMD = "du -sh"  # the actual command that's being run


def run_bash_cmd(cmd):
    """Use subprocess to capture the text output of a bash command"""
    try:
        stdout = subprocess.run(cmd, check=True, text=True, 
                                universal_newlines=True, shell=True, 
                                stdout=subprocess.PIPE).stdout
    except subprocess.CalledProcessError as e:
        stdout = f"FAILED '{cmd}'"
    return stdout


# Mark time for log and file name
time_start = time.time()
time_now = time.asctime().upper()
time_str = time.strftime("%Y-%m-%d")

# Input and output paths
base_directory = "/import/c1/ERTHQUAK"  
output_file = f"/import/c1/ERTHQUAK/ERTHQUAK/bhchow/erthquak_group_du.txt"

# Select which Users to run du for 
active_users = USERS[CHOICE]
if active_users is None:
    # Select all available Users in the directory
    active_users = [os.path.basename(_) for _ in 
                    os.listdir(base_directory) if _ != "ERTHQUAK"]

# In the case of re-running this command, this will allow us to skip completed 
if os.path.exists(output_file):
    check_lines = open(output_file, "r").read()
else:
    check_lines = []
    show_storage = run_bash_cmd("show_storage")[394:]  # strip off personal info
    # Only write file header one time
    with open(output_file, "a") as f:
        f.write(f"{'=' * 80}\n")
        f.write(f"AUTOMATED CHINOOK STORAGE OUTPUT LOG (CREATED: {time_now})\n")
        f.write(f"{'=' * 80}\n")
        f.write(f"\n")
        f.write(f"{show_storage}\n")
        f.write(f"{'=' * 80}\n")

# For each active User, run $ du -sh to get total file count, this will take 
# a while. Sort active users list so it's the same each time.
for active_user in sorted(active_users):
    # Allow skipping over active users that have already been run
    if active_user in check_lines:
        continue
    path = os.path.join(base_directory, active_user)
    output = run_bash_cmd(f"{CMD} {path}").strip()
    with open(output_file, "a") as f:
        f.write(f"{output}\n")

# Finish up by writing how long this took
with open(output_file, "a") as f:
    f.write(f"{'=' * 80}\n")
    f.write(f"Completed in {time.time() - time_start:.2f}s\n")


