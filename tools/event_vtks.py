"""
Script to make the .vtk files for each event from CMTSOLUTION files
"""
import os
import glob
import numpy as np
from pyatoa.utils.operations.source_receiver import lonlat_utm

vtk_header = ("# vtk DataFile Version 2.0\n"
              "Source and Receiver VTK file from Pyatoa\n"
              "ASCII\n"
              "DATASET POLYDATA\n"
              "POINTS\t1 float\n")

specfem_data = "./"
cmtsolutions = glob.glob(os.path.join(specfem_data, 'CMTSOLUTION_*'))
cmtsolutions.sort()
for i, cmtsolution_fid in enumerate(cmtsolutions):
    cmtsolution = np.genfromtxt(cmtsolution_fid, dtype='str', skip_header=1,
                                delimiter=':')
    event_id = cmtsolution[0][1].strip()
    event_lat = float(cmtsolution[3][1].strip())
    event_lon = float(cmtsolution[4][1].strip())
    event_depth = float(cmtsolution[5][1].strip()) * -1E3
    event_x, event_y = lonlat_utm(event_lon, event_lat, -60)

    with open(f"{event_id}.vtk", "w") as f:
        f.write(vtk_header)
        f.write(f"{event_x:18.6E}{event_y:18.6E}{event_depth:18.6E}\n\n")
        
        

