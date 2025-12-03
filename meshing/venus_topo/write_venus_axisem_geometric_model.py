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
from scipy.ndimage import gaussian_filter


CMAP = "cividis"


def build_venus_topo(lats, lons, fid, plot=False):
    """
    Read in VenusTopo file and export coordinates in the correct format
    """
    clm = pysh.datasets.Venus.VenusTopo719()  # units in m

    # Read or create radii data
    if os.path.exists(fid):
        print(f"'{fid}' exists, reading")
        radii = np.load(fid)
    else:
        # Generate the depth arrays for each line of constant longitude to get
        # a 2D array
        radii = np.zeros((nlat, nlon))
        for i, lat in enumerate(lats):
            fixed_lat = np.ones(len(lons)) * lat
            print(f"Lat={lat} {i}/{len(lats)}")
            radii[i] = clm.expand(lat=fixed_lat, lon=lons)

        # Save because calculating takes a while
        np.save(fid, radii)

    # Simple imshow to confirm resolution
    if plot:
        # Flip because our ordering starts from south pole, convert m -> km
        radii_plot = np.flipud(radii) * 1E-3
        im = plt.imshow(radii_plot, cmap=CMAP)
        plt.xlabel("nLat")
        plt.ylabel("nLon")
        plt.title(f"VenusTopo (nlat={nlat}, nlon={nlon})")
        plt.colorbar(im, fraction=0.024, pad=0.005, label="Radius [km]")
        plt.show()

    return radii


def write_data_to_netcdf_file(lats, lons, fid, depth_m=None, radii_m=None, 
                              ref_val_m=6052E3, smooth_sigma=None):
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
    assert((depth_m is not None) or (radii_m is not None))

    # Smooth the entire grid with a Gaussian however at the poles you will have
    # offset values due to the tails of the Gaussian
    if smooth_sigma is not None: 
        depth_m = gaussian_filter(depth_m, smooth_sigma)

    if os.path.exists(fid):
        print(f"'{fid}' already exists, skipping")
        return

    try:    
        ncfile = Dataset(fid, "w", format="NETCDF4_CLASSIC")
    except PermissionError:
        Dataset(fid).close()
        ncfile = Dataset(fid, "w", format="NETCDF4_CLASSIC")
    
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
    if radii_m is not None:
        radius_points = radii_m - ref_val_m 
    else:
        radius_points = ref_val_m + depth_m

    rel_topo[:] =  radius_points
    
    ncfile.close()
    print("nc file written successfully")


def read_ncfile(fid):
    """
    Read in the ncfile and return its outputs
    """
    model  = Dataset(fid)
    latitude = model.variables["latitude"]
    longitude = model.variables["longitude"]
    rel_topo  = model.variables["relative_topo_radius"]

    return longitude, latitude, rel_topo


def plot_ncfile(fid):
    """
    Read in a target ncfile and plot the values in there for confirmation.
    Assuming parameter keys are the same as the written file
    """
    x, y, z = read_ncfile(fid)
    X, Y = np.meshgrid(x, y)

    pc = plt.pcolor(X, Y, z[:] * 1E-3)  # m -> km
    plt.colorbar(pc, fraction=0.024, pad=0.005, label="Relative Topo [km]")

    plt.title("ncFile Venus Topo")
    plt.xlabel("Latitude")
    plt.ylabel("Longitude")
    plt.gca().set_aspect("equal")
    plt.show()
    

if __name__ == "__main__":
    # PARAMETERS
    nlon = 360  
    nlat = 180
    downsample = 0.25
    tag =  "venus_topo"
    plot = True
    save = True

    # --- 
    # downsample
    nlat = int(nlat / downsample)
    nlon = int(nlon / downsample)
    print(f"nlat={nlat}; nlon={nlon}")
    
    # build filename
    fid = f"{tag}_lat{nlat}_lon{nlon}"
    npfid = f"{fid}.npy"
    ncfid = f"{fid}.nc"
    print(tag)

    # define coordinates, do not wrap around
    lats = np.linspace(-89.999, 89.999, nlat)
    lons = np.linspace(-179.999, 179.999, nlon)
    # ---
    radii = build_venus_topo(lats, lons, npfid, plot=False)
    output = write_data_to_netcdf_file(lats, lons, ncfid, radii_m=radii)
    
    # last check to see that the file we wrote looks good
    if plot:
        plot_ncfile(ncfid)
