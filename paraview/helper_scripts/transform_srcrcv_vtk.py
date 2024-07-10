"""
If I set the origin of the mesh to 0,0,0, then sources and receivers also
require their locations changed. This will read in the file, convert and write
"""
import os
import sys
import numpy as np
from glob import glob

fid_in = sys.argv[1]

# Convert to 0 origin 
x_origin = 171312.0
y_origin = 5286950.0
convert = 1E-3

# Set constant values for each coordinate
set_x = None
set_y = None  # Set fixed Y-value, None not set
set_z = 1

template = "{x:18.6E}{y:18.6E}{z:18.6E}\n"
fids = glob(os.path.join("./", fid_in))
for fid in fids:
    old_file_lines = open(fid, "r").readlines()
    new_fid = os.path.basename(fid).split(".")[0] + "_converted.vtk"
    print(new_fid)
    with open(new_fid, "w") as f:
        for line in old_file_lines:
            try:
                x, y, z = line.strip().split()
                # Shift by origin and convert
                if set_x is not None:
                    x = set_x
                else:
                    x = (float(x) - x_origin) * convert
                if set_y is not None:
                    y = set_y
                else:
                    y = (float(y) - y_origin) * convert
                if set_z is not None:
                    z = set_z
                else:
                    z = float(z) * convert

                f.write(template.format(x=x, y=y, z=z))
            # ValueError thrown for header files
            except ValueError:
                f.write(line)

