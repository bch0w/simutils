"""
A few shared functions that are useful for creating meshes
"""
import numpy as np


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
        roundout = int(base * round(float(x) / base))
    elif choice == 'down':
        roundout = int(base * np.floor(float(x) / base))
    elif choice == 'up':
        roundout = int(base * np.ceil(float(x) / base))

    return roundout


def lonlat_utm(lon_or_x, lat_or_y, utm_zone=-60, inverse=False):
    """
    convert latitude and longitude coordinates to UTM projection
    From Pyatoa

    :type lon_or_x: float or int
    :param lon_or_x: longitude value in WGS84 or X in UTM-'zone' projection
    :type lat_or_y: float or int
    :param lat_or_y: latude value in WGS84 or Y in UTM-'zone' projection
    :type utm_zone: int
    :param utm_zone: UTM zone for conversion from WGS84
    :type inverse: bool
    :param inverse: if inverse == False, latlon => UTM, vice versa.
    :rtype x_or_lon: float
    :return x_or_lon: x coordinate in UTM or longitude in WGS84
    :rtype y_or_lat: float
    :return y_or_lat: y coordinate in UTM or latitude in WGS84
    """
    from pyproj import Proj

    # Determine if the projection is north or south
    if utm_zone < 0:
        direction = "south"
    else:
        direction = "north"
    # Proj doesn't accept negative zones
    utm_zone = abs(utm_zone)

    # Proj requires a string to tell it how to convert the coordinates
    projstr = (f"+proj=utm +zone={utm_zone}, +{direction} +ellps=WGS84"
               " +datum=WGS84 +units=m +no_defs")

    # Initiate a Proj object and convert the coordinates
    my_proj = Proj(projstr)
    x_or_lon, y_or_lat = my_proj(lon_or_x, lat_or_y, inverse=inverse)

    return x_or_lon, y_or_lat
