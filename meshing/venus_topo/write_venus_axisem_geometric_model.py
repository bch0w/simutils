"""
Read in VenusTopo and write out a geometric model for the surface topography
Following AxiSEM user manual and building off Python scripts provided:
https://github.com/benjaminfernando/axisem3d-3dmodels/blob/main/template_input/\
    python%20scripts/moho.py
"""
import os
import numpy as np
import matplotlib.pyplot as plt 
import pyshtools as pysh
from netCDF4 import Dataset
from scipy.ndimage import gaussian_filter, zoom
from matplotlib.ticker import FixedFormatter


def build_venus_topo(nlat=180, nlon=360, downsample=4, save=True, plot=False):
    """
    Read in VenusTopo file and export coordinates in the correct format
    """
    clm = pysh.datasets.Venus.VenusTopo719()  # units in m

    # downsample
    nlat = int(nlat / downsample)
    nlon = int(nlon / downsample)
    print(f"nlat={nlat}; nlon={nlon}")

    # Define coordinates, do not wrap around
    lats = np.linspace(-89.999, 89.999, nlat)
    lons = np.linspace(-179.999, 179.999, nlon)

    # Used for saving and checking for saved files
    fid = f"venus_topo_lat{nlat}_lon{nlon}.npy"

    # Read or create radii data
    if os.path.exists(fid):
        radii = np.load(fid)
    else:
        # Generate the depth arrays for each line of constant longitude to get
        # a 2D array
        radii = np.zeros((nlon, nlat))
        for i, lon in enumerate(lons):
            fixed_lon = np.ones(len(lats)) * lon
            print(f"Lon={lon} {i}/{len(lons)}")
            radii[i] = clm.expand(lat=lats, lon=fixed_lon)

        # Save because calculating takes a while
        if save:
            np.save(radii, fid)

    # Simple imshow to confirm resolution
    if plot:
        plt.imshow(np.flipup(radii.T))
        plt.title(f"nlat={nlat}, nlon={nlon}")
        plt.show()

    return lats, lons, radii


def write_data_to_netcdf_file(lats, lons, depth_m=None, radii_m=None, 
                              ref_val_m=6052E3, smooth_sigma=None, 
                              output="venustopo.nc",):
    """
    Writes lat/lon/depth to NetCDF used by AxiSEM3D for topo implementation.
    Modified from AxiSEM3D example `moho.py`

    :type lat: np.array
    :param lat: list of lat values, needs to start at -90 to work with AxiSEM
    :type lat: np.array
    :param lat: list of lon values, needs to start at -180 to work with AxiSEM
    :type radii_m: np.array
    :param radii_m: 2D array of radii corresponding to the above lat/lons
    :type ref_val_m: float
    :type ref_val_m: radius of the layer to undulate, could be surface, or moho
    :type smooth_sigma: float
    :param smooth_sigma: optional std. for gaussian smoothing, if not given then
        no smoothing
    :type output: str
    :param output: filename for file to be written
    """
    assert(depth_m or radii_m)

    # Smooth the entire grid with a Gaussian however at the poles you will have
    # offset values due to the tails of the Gaussian
    if smooth_sigma is not None: 
        depth_m = gaussian_filter(depth_m, smooth_sigma)

    try:    
        ncfile = Dataset(output, "w", format="NETCDF4_CLASSIC")
    except PermissionError:
        Dataset(output).close()
        ncfile = Dataset(output, "w", format="NETCDF4_CLASSIC")
    
    # Let NetCDF know the dimensions of the output file
    ncfile.createDimension("latitude", len(lats))  # CHECK THIS
    ncfile.createDimension("longitude", len(lons))
    ncfile.createDimension("relative_topo_radius", 0)  # we will add to this

    # NetCDF Variables
    latitude = ncfile.createVariable("latitude", "d", ("latitude",)) 
    longitude = ncfile.createVariable("longitude", "d", ("longitude",))
    rel_topo = ncfile.createVariable("relative_topo_radius", "d",
                                     ("latitude","longitude")
                                     )

    # Add data to the dimensions
    latitude[:] = lats
    longitude[:] = lons
    
    # Ensure that the actual depth values are set at the correct radius
    if radii_m:
        radius_points = ref_val_m - radii_m
    else:
        radius_points = ref_val_m + depth_m
    rel_topo[:] = np.ravel(radius_points)
    
    ncfile.close()


def plot_ncfile(fid):
    """
    Read in a target ncfile and plot the values in there for confirmation    
    """
    ncfile  = Dataset(fid)
    


if __name__ == "__main__":
    lats, lons, radii = build_venus_topo()
    write_data_to_netcdf_file(lats, lons, radii)
