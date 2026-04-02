"""
Create a uniform grid of points for a given domain and output a Specfem3D
STATIONS file with unique station naming 
"""
import numpy as np
from pyproj import Proj

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
    p = Proj(proj="utm", zone=utm_zone, ellps="WGS84", preserve_units=False)
    x_or_lon, y_or_lat = p(lon_or_x, lat_or_y, inverse=inverse)

    return x_or_lon, y_or_lat


def uniform_grid_latlon(lat_min, lat_max, lon_min, lon_max, dlat, dlon):
    """
    Return a uniform grid of points
    """
    lats = np.arange(lat_min, lat_max, dlat)
    lons = np.arange(lon_min, lon_max, dlon)
    lat_grid, lon_grid = np.meshgrid(lats, lons)

    return lat_grid, lon_grid


def write_to_stations(lat_grid, lon_grid, network="NN", sta_tag="S{:0>3}",
                      fid="./STATIONS"):
    """
    Given two 2D arrays, write a STATIONS file with unique station naming
    """
    template = "{s:>6}{n:>6}{lat:12.4f}{lon:12.4f}{d:7.1f}{e:7.1f}\n"
    i = 1
    with open(fid, "w") as f:
        for lat, lon in zip(lat_grid, lon_grid):
            for lat_, lon_ in zip(lat, lon):
                f.write(template.format(s=sta_tag.format(i), n=network,
                                        lat=lat_, lon=lon_, d=0., e=0.
                                        ))
                i += 1
    print(f"{i} stations written")


if __name__ == "__main__":
    kwargs = {"lat_min": 4483695.801447477,
              "lat_max": 4683695.801447477, 
              "lon_min": 418796.1202446909,
              "lon_max": 568796.1202446909,
              "dlat": 0.05,
              "dlon": 0.05}
    
    utm_zone = 52  # if not None, conerts coordinates to or from lat lon
    if utm_zone is not None:
        kwargs["lon_min"], kwargs["lat_min"] = lonlat_utm(kwargs["lon_min"],
                                                          kwargs["lat_min"],
                                                          utm_zone,
                                                          inverse=True)

        kwargs["lon_max"], kwargs["lat_max"] = lonlat_utm(kwargs["lon_max"],
                                                          kwargs["lat_max"],
                                                          utm_zone,
                                                          inverse=True)

    lat_grid, lon_grid = uniform_grid_latlon(**kwargs)
    write_to_stations(lat_grid, lon_grid)




