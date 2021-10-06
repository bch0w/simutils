"""
Change the POINTS quantity so that different VTK files don't force the same
colorscale on each other
"""
import sys

assert(len(sys.argv) == 3), "input requires 'filename' and 'quantity'"

with open(sys.argv[1]) as f:
    lines = f.readlines()

for i, line in enumerate(lines[:]):
    if "SCALARS" in line:
        scalars, qty, float_ = line.split()
        print(f"{qty} -> {sys.argv[2]}")
        lines[i] = line.replace(qty, sys.argv[2])
        break

with open(sys.argv[1], "w") as f:
    f.writelines(lines)

