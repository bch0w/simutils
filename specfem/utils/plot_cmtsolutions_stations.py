"""
A simple script to visualize the station and receiver locations on a map
with event id and station name annotation
"""
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['lines.markersize'] = 5
mpl.rcParams['xtick.labelsize'] = 'x-small'
mpl.rcParams['ytick.labelsize'] = 'x-small'


def myround(x, base=5, choice='near'):
    """
    Round value x to nearest base, round 'up','down' or to 'near'est base

    :type x: float
    :param x: value to be rounded
    :type base: int
    :param base: nearest integer to be rounded to
    :type choice: str
    :param choice: method of rounding, 'up', 'down' or 'near'
    :rtype roundout: int
    :return: rounded value
    """
    if choice == 'near':
        roundout = int(base * round(float(x)/base))
    elif choice == 'down':
        roundout = int(base * np.floor(float(x)/base))
    elif choice == 'up':
        roundout = int(base * np.ceil(float(x)/base))

    return roundout


def utm_zone_from_lat_lon(lat, lon):
    """
    Calculate the UTM zone longitude value using quick maffs.
    Get the sign of the UTM zone based on the latitude value.

    :type lat: float
    :param lat: latitude coordinate in degrees
    :type lon: float
    :param lon: longitude coordinate in degrees
    :rtype: int
    :return: UTM zone number
    """
    try:
        sign = lat / abs(lat)  # silly way to figure out if lat is +/-
    except ZeroDivisionError as e:
        raise Exception("latitude is 0, UTM zone is ambigious") from e
    return sign * np.ceil((lon + 180) / 6)


def lonlat_utm(lon_or_x, lat_or_y, utm_zone=None, inverse=False):
    """
    Convert latitude and longitude coordinates to UTM projection using PyProj

    :type lon_or_x: float or int
    :param lon_or_x: longitude value in WGS84 or X in UTM-'zone' projection
    :type lat_or_y: float or int
    :param lat_or_y: latude value in WGS84 or Y in UTM-'zone' projection
    :type utm_zone: int
    :param utm_zone: UTM zone for conversion from WGS84
    :type inverse: bool
    :param inverse: if inverse == False, latlon => UTM, vice versa.
    :rtype: tuple (float, float)
    :return: (x in UTM or longitude in WGS84, y in UTM or latitude in WGS84)
    """
    from pyproj import Proj

    # If converting latlon to utm and no utm zone given, calculate utm zone
    if utm_zone is None and not inverse:
        utm_zone = utm_zone_from_lat_lon(lat_or_y, lon_or_x)
    elif utm_zone is None and inverse:
        raise TypeError(
            "lonlat_utm() missing 1 required positional argument: 'utm_zone'"
        )
    # Determine if the projection is north or south
    if utm_zone < 0:
        direction = "south"
    else:
        direction = "north"

    # Proj doesn't accept negative zones
    utm_zone = abs(utm_zone)

    projstr = (f"+proj=utm +zone={utm_zone} +{direction} +ellps=WGS84"
               " +datum=WGS84 +units=m +no_defs")
    projection = Proj(projstr)

    x_or_lon, y_or_lat = projection(lon_or_x, lat_or_y, inverse=inverse)

    return x_or_lon, y_or_lat


# =============================================================================
specfem_data = "./"
lat_min = 40.5
lat_max = 45.5
lon_min = 126.
lon_max = 134.
utm_zone = 52
plot_cmtsolutions = True
plot_stations = True
fontsize = 4
# =============================================================================

f, ax = plt.subplots(1, dpi=200)

# Set the corners of the map
x_min, y_min = lonlat_utm(lon_min, lat_min, utm_zone)
x_max, y_max = lonlat_utm(lon_max, lat_max, utm_zone)

# Set major and minor ticks for grid
if False:
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
    # ax.grid(which='both')

    # Set plot boundaries
    plt.xlim([x_min, x_max])
    plt.ylim([y_min, y_max])

if plot_stations:
    # Gather all stations
    station_file = os.path.join(specfem_data, "STATIONS")
    stations = np.genfromtxt(station_file, dtype='str')

    for i, station in enumerate(stations):
        # station_id = ".".join([station[1], station[0]])
        station_id = station[0]
        
        # convert to UTM 
        station_lat = float(station[2])
        station_lon = float(station[3])
        x, y = lonlat_utm(station_lon, station_lat, utm_zone)
        
        plt.scatter(x, y, c='coral', marker='v', edgecolor='k')
        plt.annotate(xy=(x,y), text=station_id, fontsize=fontsize, 
                     color="coral")

# Gather all the CMTSOLUTIONS
if plot_cmtsolutions:
    cmtsolutions = glob.glob(os.path.join(specfem_data, "CMTSOLUTION*"))
    cmtsolutions.sort()
    for i, cmtsolution_fid in enumerate(cmtsolutions):
        cmtsolution = np.genfromtxt(cmtsolution_fid, dtype='str', skip_header=1, 
                                    delimiter=':')
        event_id = cmtsolution[0][1].strip()
        event_lat = float(cmtsolution[3][1].strip())
        event_lon = float(cmtsolution[4][1].strip())

        # get magnitude from moment tensor
        moment = 0
        for m in range(6, 12):
            moment += float(cmtsolution[m][1].strip()) ** 2
        moment = 1/np.sqrt(2) * np.sqrt(moment) 
        event_mw = 2/3 * np.log10(abs(moment)) - 10.7

        # convert to UTM -60
        x, y = lonlat_utm(event_lon, event_lat, utm_zone)

        # Plot as a circle, annotate name
        plt.scatter(x, y, c='teal', marker='o', edgecolor='k')
        anno_x = 0.75 * (x_max - x_min) + x_min
        anno_y = np.linspace(y_min, 0.6 * (y_max -y_min) + y_min, 
                                                          len(cmtsolutions) + 1)
        anno_y = anno_y[::-1]
        plt.annotate(xy=(x,y), text=i, fontsize=10, color="teal")
        plt.annotate(xy=(x,y), xytext=(anno_x, anno_y[i]),  
                     text=f"{i}: M{event_mw:.1f}; {event_id}", fontsize=fontsize, 
                     color="teal")

plt.savefig('source_receivers.png')
plt.show()
    
    


