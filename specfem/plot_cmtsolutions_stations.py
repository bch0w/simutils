"""
A simple script to visualize the station and receiver locations on a map
with event id and station name annotation
"""
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from pyatoa.utils.operations.source_receiver import lonlat_utm

mpl.rcParams['lines.markersize'] = 4
fontsize=6

specfem_data = '/home/chowbr/tomo/seisflows/checkerboard/eightevent_1603deb/scratch/solver/2012p783751/DATA'
# specfem_data = '/scale_wlg_nobackup/filesets/nobackup/nesi00263/bchow/tomo/seisflows/specfem3d_1603deb/DATA_seisflows'
coastline_npy = '/scale_wlg_persistent/filesets/home/chowbr/primer/auxiliary/coastline/nz_resf_coast_mod_utm60H_xyz.npy'


f = plt.figure(dpi=100)

# Set the corners of the map
lat_min = -43
lat_max = -37
lon_min = 173
lon_max = 179
x_min, y_min = lonlat_utm(lon_min, lat_min, -60)
x_max, y_max = lonlat_utm(lon_max, lat_max, -60)

plt.xlim([x_min, x_max])
plt.ylim([y_min, y_max])

# Gather all the CMTSOLUTIONS
cmtsolutions = glob.glob(os.path.join(specfem_data, 'CMTSOLUTION_*'))
for i, cmtsolution_fid in enumerate(cmtsolutions):
    cmtsolution = np.genfromtxt(cmtsolution_fid, dtype='str', skip_header=1, 
                                delimiter=':')
    event_id = cmtsolution[0][1].strip()
    event_lat = float(cmtsolution[3][1].strip())
    event_lon = float(cmtsolution[4][1].strip())

    # convert to UTM -60
    x, y = lonlat_utm(event_lon, event_lat, -60)

    # Plot as a circle, annotate name
    plt.scatter(x, y, c='None', marker='o', edgecolor='k')
    anno_x = 0.8 * (x_max - x_min) + x_min
    anno_y = np.linspace(y_min, 0.5 * (y_max -y_min) + y_min, len(cmtsolutions) + 1)
    plt.annotate(xy=(x,y), s=i, fontsize=10)
    plt.annotate(xy=(x,y), xytext=(anno_x, anno_y[i]),  s=f"{i}: {event_id}", 
                 fontsize=fontsize)
    
# Gather all stations
station_file = os.path.join(specfem_data, 'STATIONS')
stations = np.genfromtxt(station_file, dtype='str')

for station in stations:
    station_ids = ".".join([station[1], station[0]])
    
    # convert to UTM -60
    station_lat = float(station[2])
    station_lon = float(station[3])
    x, y = lonlat_utm(station_lon, station_lat, -60)
    
    plt.scatter(x, y, c='None', marker='v', edgecolor='k')
    plt.annotate(xy=(x,y), s=station_ids, fontsize=fontsize)

# Plot any auxiliary data
if coastline_npy:
    coastline = np.load(coastline_npy)
    plt.scatter(coastline[:,0], coastline[:,1], c='k', marker='.', 
                s=0.5, alpha=0.25, zorder=2)

plt.savefig('source_receivers.png')
    
    


