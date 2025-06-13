"""
Create a uniform grid of points for a given domain and output a Specfem3D
STATIONS file with unique station naming 
"""
import numpy as np


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

if __name__ == "__main__":
    kwargs = {"lat_min": 40.0,
              "lat_max": 46.0,
              "lon_min": 125.0,
              "lon_max": 135.0,
              "dlat": 0.5,
              "dlon": 0.5}

    lat_grid, lon_grid = uniform_grid_latlon(**kwargs)
    write_to_stations(lat_grid, lon_grid)




