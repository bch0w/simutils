"""
Create a uniform grid of points for a given domain and output a Specfem3D
STATIONS file with unique station naming 

Note that station numbering will be 6 characters long,
"""
import numpy as np
from pyproj import Proj

# Expected format for STATIONS file
TEMPLATE = "{s:>6}{n:>6}{lat:12.4f}{lon:12.4f}{d:7.1f}{e:7.1f}\n"


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


def write_to_stations(lat_grid, lon_grid, network="NN", sta_tag="{:0>5}",
                      fid="./STATIONS"):
    """
    Given two 2D arrays, write a STATIONS file with unique station naming
    """
    i = 1
    with open(fid, "w") as f:
        for lat, lon in zip(lat_grid, lon_grid):
            for lat_, lon_ in zip(lat, lon):
                f.write(template.format(s=sta_tag.format(i), n=network,
                                        lat=lat_, lon=lon_, d=0., e=0.
                                        ))
                i += 1
    print(f"{i} stations written")


def main():
    choice = "NK2026_PAPER"

    if choice == "NK2026_PAPER":
        lon_min = 420_000
        lon_max = 565_000
        lat_min = 4_490_000
        lat_max = 4_680_000
        nlatlon = (100, 100)
        dlatlon = None
        utm_zone = 52
        network = "XX"
        fid_out = ("/home/bhchow/REPOS/spectral/research/watc/simblast/"
                   "SPECFEM_DATA/STATIONS/STATIONS_PAPER_NK_GRID")
    else:
        lat_min = 0
        lat_max = 0
        lon_min = 0
        lon_max = 0
        nlatlon = None
        dlatlon = (0, 0)
        utm_zone = None
        network = "XX"
        fid_out = "STATIONS"

    # Convert from UTM to Lon/Lat 
    if utm_zone:
        lon_min, lat_min = lonlat_utm(lon_min, lat_min, utm_zone, True)
        lon_max, lat_max = lonlat_utm(lon_max, lat_max, utm_zone, True)

    #  Create grid
    if nlatlon:
        nlat, nlon = nlatlon
        lats = np.linspace(lat_min, lat_max, nlat)
        lons = np.linspace(lon_min, lon_max, nlon)
    elif dlatlon:
        dlat, dlon = dlatlon
        lats = np.arange(lat_min, lat_max, dlat)
        lons = np.arange(lon_min, lon_max, dlon)
    print(f"writing {len(lats)*len(lons)} stations")

    with open(fid_out, "w") as f:
        for i, lat in enumerate(lats):
            for j, lon in enumerate(lons):
                station_name = f"{i:0>2}{j:0>2}"  # e.g., 0000 is the ULHC
                f.write(TEMPLATE.format(s=station_name, n=network,
                                        lat=lat, lon=lon, d=0., e=0.))


if __name__ == "__main__":
    main()




