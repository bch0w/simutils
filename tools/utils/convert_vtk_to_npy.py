"""
Simple script to convert a .vtk file into numpy arrays containing the 
coordinates and the corresponding values. Praise be to numpy, slayer of the 
devil Paraview, that wicked program. Fuck me paraview is so annoying.
"""
import sys
import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy

fid = sys.argv[1]
tag = fid.split(".")[0]

# Instantiate the reader
reader = vtk.vtkGenericDataObjectReader()
reader.SetFileName(fid)
reader.Update()

# Grab coordinates and values from the read in file
coords = vtk_to_numpy(reader.GetOutput().GetPoints().GetData())
values = vtk_to_numpy(reader.GetOutput().GetPointData().GetScalars())

# Append values as a column to the coordinate file
arr_out = np.c_[coords, values]

# Save
np.save(f"{tag}", arr_out)
