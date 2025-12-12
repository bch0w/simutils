"""
Create a uniform grid of points for a given domain and output an AxiSEM3D 
STATIONS file with unique station naming 
"""
import numpy as np
# import pyshtools as pysh


def uniform_grid_latlon(lat_min, lat_max, lon_min, lon_max, dlat, dlon):
    """
    Return a uniform grid of points
    """
    lats = np.arange(lat_min, lat_max, dlat)
    lons = np.arange(lon_min, lon_max, dlon)
    lat_grid, lon_grid = np.meshgrid(lats, lons)

    return lat_grid, lon_grid


def write_to_stations(lat_grid, lon_grid, network="NN", sta_tag="S{:0>3}",
                      elevation=None, fid="./STATIONS", _print=False):
    """
    Given two 2D arrays, write a STATIONS file with unique station naming
    """
    if _print:
        clm = pysh.datasets.Venus.VenusTopo719() / 1E3
    else:
        clm = None

    template = "{s:>6}{n:>6}{lat:12.4f}{lon:12.4f}{d:7.1f}{e:9.2f}\n"
    i = 1
    with open(fid, "w") as f:
        for lat, lon in zip(lat_grid, lon_grid):
            for lat_, lon_ in zip(lat, lon):
                write_str = template.format(s=sta_tag.format(i), n=network,
                                            lat=lat_, lon=lon_, d=0., e=0.)
                if _print:
                    ref_radius = 6051.9
                    # Print for posterity and to show the relative elevation for 
                    # each station
                    elv = (clm.expand(lat=lat_, lon=lon_) - ref_radius) * 1E3
                    prntscrn = template.format(s=sta_tag.format(i), n=network,
                                               lat=lat_, lon=lon_, d=0, 
                                               e=elv)[:-1]
                    print(prntscrn)

                
                f.write(write_str)
                i += 1
    print(f"{i} stations written")


def write_grid_to_stations(lat_min, lat_max, lon_min, lon_max, dlat, dlon, 
                           network="NN", elevation=None, fid="./STATIONS",):
    """
    Given two 2D arrays, write a STATIONS file with unique station naming
    """
    template = "{s:>6}{n:>6}{lat:12.4f}{lon:12.4f}{d:7.1f}{e:9.2f}\n"

    if elevation is None:
        clm = pysh.datasets.Venus.VenusTopo719() / 1E3
    else:
        clm = None

    lats = np.arange(lat_min, lat_max, dlat)
    lons = np.arange(lon_min, lon_max, dlon)
    c = 0
    with open(fid, "w") as f:
        for i, lat_ in enumerate(lats):
            for j, lon_ in enumerate(lons):
                s = f"{i:0>2}{j:0>2}"

                if clm:
                    ref_radius = 6051.9
                    # Print for posterity and to show the relative elevation for 
                    # each station
                    elv = (clm.expand(lat=lat_, lon=lon_) - ref_radius) * 1E3
                else:
                    elv = 0

                write_str = template.format(s=s, n=network, lat=lat_, lon=lon_, 
                                            d=0., e=elv)
                
                f.write(write_str)
                c += 1
    print(f"{c} stations written")

if __name__ == "__main__":
    network = "VN"
    kwargs = {"lat_min": -30.,
              "lat_max": 60.,
              "lon_min": -140.,
              "lon_max": -20.0,
              "dlat": 1,
              "dlon": 1}
    if True:
        write_grid_to_stations(network=network, fid="./STATIONS_GRID", 
                               elevation=False, **kwargs)

    if False:
        lat_grid, lon_grid = uniform_grid_latlon(**kwargs)
        write_to_stations(lat_grid, lon_grid, network, _print=True)




