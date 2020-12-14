"""
Convert an XYZ point cloud into something renderable by paraview
Point cloud needs to be defined by X, Y and Z values that have the same length
"""
import numpy as np

x, y, z = np.loadtxt("williams_hikurangi_interface_utm60_nonan.xyz").T
npoints = len(x)

with open("interface.vtk", "w") as f:
    # Write header
    f.write("# vtk DataFile Version 2.0\n"
            "loop\n"
            "ASCII\n"
            "DATASET UNSTRUCTURED_GRID\n"
            f"POINTS {npoints} float\n")

    # Write data
    for x_, y_, z_ in zip(x, y, z):
        f.write(f"{x_:15.2f} {y_:15.2f} {z_:15.2f}\n")

    # Write cell data
    f.write(f"CELLS {npoints} {npoints*2}\n")
    for i in range(npoints):
        f.write(f"1 {i}\n")

    # Write cell types
    f.write(f"CELL_TYPES {npoints}\n")
    for i in range(npoints):
        f.write(f"1 \n")

    # Write Z Scalar values
    f.write(f"\nPOINT_DATA {npoints}\n"
            f"SCALARS Z_Value float 1\n"
            f"LOOKUP_TABLE default\n")
    for z_ in z:
        f.write(f"{z_}\n")

