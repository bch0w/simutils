"""
Python wrapper for Bash commands to check storage on UAF HPC Chinook by active
User. Intended to be run by a weekly Cronjob with email output to a Users
group. 


Contains built in leeway for fail halfway through and rerunning the process w/o
modification

.. note::

    Rule of thumb 'du -sh' takes 45m/TB (based on ~2TB taking 1h30m)
"""
import os
import sys
import subprocess
import time

# Small dry-run version for testing/debugging
TEST = True

def run_bash_cmd(cmd):
    """Use subprocess to capture the text output of a bash command"""
    return subprocess.run(cmd, check=True, text=True, universal_newlines=True, 
                          shell=True, stdout=subprocess.PIPE).stdout

# Mark time for log and file name
time_now = time.asctime().upper()
time_str = time.strftime("%Y-%m-%d")

# Determine where we are checking files and which users
base_directory = "/import/c1/ERTHQUAK"
if TEST:
    active_users = ["ykaneko1_ua"]
else:
    active_users = ["bhchow", "ctape", "ammcpherson", "jthurin", "agupta7", 
                    "eaonyango", "ytian4"]

# Write file somewhere so we can reference/email later
output_file = f"/import/c1/ERTHQUAK/ERTHQUAK/bhchow/erthquak_du_{time_str}.txt"

# In the case of re-running this command, this will allow us to not rerun
if os.path.exists(output_file):
    check_lines = open(output_file, "r").read()
else:
    check_lines = []
    show_storage = run_bash_cmd("show_storage")[77:]  # strip off first sentence
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
    output = run_bash_cmd(f"du -sh {path}")
    with open(output_file, "a") as f:
        f.write(f"{output}\n")


