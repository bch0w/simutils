"""
A simple script to visualize the station and receiver locations on a map
with event id and station name annotation
"""
import os
import glob
import numpy as np
from pyatoa.utils.operations.calculations import myround
from pyatoa.utils.operations.source_receiver import lonlat_utm

specfem_data = "./"

cmtsolutions = glob.glob(os.path.join(specfem_data, 'CMTSOLUTION_*'))
for i, cmtsolution_fid in enumerate(cmtsolutions):
    cmtsolution = np.genfromtxt(cmtsolution_fid, dtype='str', skip_header=1, 
                                delimiter=':')
    event_id = cmtsolution[0][1].strip()
    event_lat = float(cmtsolution[3][1].strip())
    event_lon = float(cmtsolution[4][1].strip())
    z = float(cmtsolution[5][1].strip())
    z *= -1E3

    # convert to UTM -60
    x, y = lonlat_utm(event_lon, event_lat, -60)
    print(f"{x:18.6E}{y:18.6E}{z:18.6E}")
