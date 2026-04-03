"""
SPECFEM Tomo Helper does not allow for attenuation in the output xyz models
so this script just tags on a Qp and Qs value to the end. Header does not 
need to be adjusted. Writes in place
"""
import sys

QP = 350
QS = 350

fid = sys.argv[1]
with open(fid, "r") as f:
    lines = f.readlines()

with open(fid, "w") as f:
    for line in lines[0:4]:
        f.write(line)
    for line in lines[4:]:
        line = f"{line.strip()} {QP:4.1f} {QS:4.1f}\n"
        f.write(line)
        


