"""
A simple script to visualize the station and receiver locations on a map
with event id and station name annotation
"""
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from pyatoa.utils.operations.calculations import myround
from pyatoa.utils.operations.source_receiver import lonlat_utm

mpl.rcParams['lines.markersize'] = 5
mpl.rcParams['xtick.labelsize'] = 'x-small'
mpl.rcParams['ytick.labelsize'] = 'x-small'

fontsize = 8.5

specfem_data = "./"
coastline_npy = "/scale_wlg_persistent/filesets/home/chowbr/primer/auxiliary/coastline/nz_resf_coast_mod_utm60H_xyz.npy"

plot_cmtsolutions = True
plot_stations = False

f, ax = plt.subplots(1, dpi=200)

# Set the corners of the map
lat_min = -43
lat_max = -37
lon_min = 173
lon_max = 179
x_min, y_min = lonlat_utm(lon_min, lat_min, -60)
x_max, y_max = lonlat_utm(lon_max, lat_max, -60)

# Set major and minor ticks for grid
major_tick_spacing_km = 50*1E3
minor_tick_spacing_km = 10*1E3
major_xticks = np.arange(myround(x_min, 10000), x_max, major_tick_spacing_km)
minor_xticks = np.arange(myround(x_min, 10000), x_max, minor_tick_spacing_km)
major_yticks = np.arange(myround(y_min, 10000), y_max, major_tick_spacing_km)
minor_yticks = np.arange(myround(y_min, 10000), y_max, minor_tick_spacing_km)
ax.set_xticks(major_xticks)
ax.set_xticks(minor_xticks, minor=True)
ax.set_yticks(major_yticks)
ax.set_yticks(minor_yticks, minor=True)
ax.grid(which='both')

# Set plot boundaries
plt.xlim([x_min, x_max])
plt.ylim([y_min, y_max])

# Gather all the CMTSOLUTIONS
if plot_cmtsolutions:
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
        anno_y = np.linspace(y_min, 0.5 * (y_max -y_min) + y_min, 
                                                          len(cmtsolutions) + 1)
        plt.annotate(xy=(x,y), s=i, fontsize=10)
        plt.annotate(xy=(x,y), xytext=(anno_x, anno_y[i]),  
                     s=f"{i}: {event_id}", fontsize=fontsize)

    
if plot_stations:
    # Gather all stations
    station_file = os.path.join(specfem_data, 'STATIONS')
    stations = np.genfromtxt(station_file, dtype='str')

    for i, station in enumerate(stations):
        station_ids = ".".join([str(i),station[1], station[0]])
        
        # convert to UTM -60
        station_lat = float(station[2])
        station_lon = float(station[3])
        x, y = lonlat_utm(station_lon, station_lat, -60)
        
        plt.scatter(x, y, c='None', marker='v', edgecolor='k')
        plt.annotate(xy=(x,y), s=station_ids, fontsize=fontsize)

# Plot any auxiliary data
if coastline_npy:
    if os.path.exists(coastline_npy):
        coastline = np.load(coastline_npy)
        plt.scatter(coastline[:,0], coastline[:,1], c='k', marker='.', 
                    s=0.5, alpha=0.25, zorder=2)

plt.show()
plt.savefig('source_receivers.png')
    
    


