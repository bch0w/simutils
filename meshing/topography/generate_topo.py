"""
A script to generate topography data files for use in Meshfem3D
Designed around the topography/bathymetry files of SRTM30P by Sandwell et al.

https://topex.ucsd.edu/WWW_html/srtm30_plus.html (6/30/25 last access)

.. note::
    Underlying topography files should be in .nc (NetCDF) format, these will be
    converted to numpy arrays internally
"""
import os
import json
import sys
import traceback
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from pyproj import Proj
from scipy.io import netcdf_file
from scipy.interpolate import griddata


def set_parameters(fid):
    """
    Load parameters from configuration file
    """
    try:
        with open(fid, "r") as f:
            parameters = json.load(f)
    except Exception as e:
        print("\n\nERROR READING JSON FILE") 
        traceback.print_exc()
        print("ERROR READING JSON FILE\n\n")
        sys.exit(-1)

    return parameters


def read_nc(fid):
    """
    Read a NetCDF file (.nc) and return the underlying data. It is assumed these
    nc files are from srtm30p and that the internal variables are: x, y, z
    given in lon, lat, meters (positive z up)

    The file is given as a grid so we also want to convert it into a single
    3-column array so it's easier to work with

    :type fid: str
    :param fid: .nc file to read
    :rtype: np.array
    :return: Nx3 array with columns representing x, y, z
    """
    with netcdf_file(fid) as data:
        try:
            x = data.variables["x"][:]
            y = data.variables["y"][:]
            z = data.variables["z"][:]

            # convert x, y, and z to MxN arrays
            x, y = np.meshgrid(x, y)

            # convert x, y, z to Nx1 arrays
            x = x.flatten()
            y = y.flatten()
            z = z.flatten()

            # convert x, y, z to a single Nx3 array
            return np.vstack((x, y, z)).T
        except KeyError:
            print("Unexpected internal structure of .nc file")
            sys.exit(-1)


def convert_coordinates(data, utm_projection=-60):
    """
    If the topography file is in lat/lon coordinates, convert to UTM projection
    Returns irregular grid since the conversion does not retain uniform spacing

    :type data: np.array
    :param data: data to convert, should be in format [lat, lon,
    """
    # Load in the all the data and distribute
    lon = data[:, 0]
    lat = data[:, 1]
    x_irreg, y_irreg = lonlat_utm(lon, lat, utm_projection)

    # Replace lon/lat with x,y
    data_out = np.copy(data)
    data_out[:, 0] = x_irreg
    data_out[:, 1] = y_irreg

    return data_out


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


def cut_topography(data, x_or_lon_min, x_or_lon_max, y_or_lat_min,
                   y_or_lat_max):
    """
    Take a topography data array and cut it to the correct size and shape of
    the mesh based on mesh dimensions and element spacing.

    :type data: np.array
    :param data: data to be cut, should  be in the format [x, y, z]
    :type x_or_lon_min: float
    :param x_or_lon_min: minimum x or longitude value
    :type x_or_lon_max: float
    :param x_or_lon_max: maximum x or longitude value
    :type y_or_lat_min: float
    :param y_or_lat_min: minimum y or latitude value
    :type y_or_lat_max: float
    :param y_or_lat_max: maximum y or latitude value
    :return:
    """
    # Determine where the data falls outside the bounds
    x_too_small = np.where(data[:, 0] < x_or_lon_min)[0]
    x_too_large = np.where(data[:, 0] > x_or_lon_max)[0]
    y_too_small = np.where(data[:, 1] < y_or_lat_min)[0]
    y_too_large = np.where(data[:, 1] > y_or_lat_max)[0]

    # Get rid of duplicates
    to_remove = np.unique(np.concatenate(
        (x_too_small, x_too_large, y_too_small, y_too_large), 0)
    )

    # delete in place
    data = np.delete(data, to_remove, 0)

    return data


def interpolate_points(data, x_min, x_max, y_min, y_max, spacing_m, plot=False):
    """
    Interpolates data to a regular grid, based on mesh dimensions

    :type data: np.array
    :param data: data to be initerpolated, [x, y, z]
    :type x_min: float
    :param x_min: minimum x value
    :type x_max: float
    :param x_max: maximum x value
    :type y_min: float
    :param y_min: minimum y value
    :type y_min: float
    :param y_max: maximum y value
    :type spacing_m: int
    :param spacing_m: spacing of the new interpolated mesh
    :return:
    """
    # Parse data
    points = data[:, :2]
    values = data[:, 2]

    # Create the regular grid to be interpolated onto
    x_reg = np.arange(x_min, x_max, spacing_m)
    y_reg = np.arange(y_min, y_max, spacing_m)

    x_grid, y_grid = np.meshgrid(x_reg, y_reg)
    x_out = x_grid.flatten()
    y_out = y_grid.flatten()

    # Interpolate the topography data
    interp_vals = griddata(points=points, values=values, xi=(x_grid, y_grid))

    if plot:
        plt.imshow(interp_vals)
        plt.gca().invert_yaxis()
        plt.title("NAlaska UTM3 SRTM30P Topography")
        plt.xlabel("X [km]")
        plt.ylabel("Y [km]")
        plt.show()

    # Create the ndarray by creating column vectors and mushing em together
    z_out = interp_vals.flatten()

    # Arrange the data in the desired format
    data_out = np.column_stack((x_out, y_out, z_out))

    return data_out


def write_data(data, method, tag, path="./"):
    """
    Saves the topo as an npy and an ascii, with a descriptive file name

    :type data: np.array
    :param data: Nx3 array representing x, y, z for topography
    :type method: str
    :param method: different write methods for 'trelis' and 'meshfem'
    :type tag: str
    :param tag: tag for the file name, extension will be appended automatically
    """
    if not os.path.exists(path):
        os.makedirs(path)

    # Meshfem only requires a single value of Z-values
    if method == "meshfem":
        fid = os.path.join(path, f"{tag}_topography.dat")
        np.savetxt(os.path.join(path, fid), data[:, 2], fmt='%d')

        topo_fid = os.path.join(path, "interface_topo.dat")
        if os.path.exists(topo_fid):
            os.remove(topo_fid)
        os.symlink(fid, topo_fid)

    # Trelis requires XYZ files. If a flat moho is required, it must be the same
    # sampling rate as the topo
    elif method == "trelis":
        np.savetxt(os.path.join(path, tag + '.xyz'), data, fmt='%.1f')


def print_meshfem_stats(data, utm_projection, log="./log.txt"):
    """
    Meshfem3D requires an interfaces.dat file which describes the format of
    the topography file to the mesher. Information required:
    npts_x, npts_y, lon_min, lat_min, d_lon, d_lat

    :type data: np.array
    :param data: Nx3 array containing x, y, z data
    :type utm_projection: int
    :param utm_projection: utm projection for conversion back into lat lon
    :return:
    """
    x, y = data[:, :2].T
    x = np.unique(x)
    y = np.unique(y)
    assert((len(x) * len(y)) == len(data)), "data lengths don't line up, weird"

    lon_min, lat_min = lonlat_utm(x.min(), y.min(), utm_projection, True)
    lon_max, lat_max = lonlat_utm(x.max(), y.max(), utm_projection, True)

    # abs() because we only care about the magnitude
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    deltax = abs(abs(x.max()) - abs(x.min()))
    deltay = abs(abs(y.max()) - abs(y.min()))
    deltalat = abs(abs(lat_max) - abs(lat_min))
    deltalon = abs(abs(lon_max) - abs(lon_min))
    dlat = deltalat / len(y)
    dlon = deltalon / len(x)

    with open(log, "w") as f:
        print(f"\nTOPOGRAPHY INFO\n{'='*50}\n"
            f"TOPO_MIN =  {data[:,2].min():.2f}\n"
            f"TOPO_MAX =  {data[:,2].max():.2f}\n"
            f"TOPO_MEAN = {data[:,2].mean():.2f}\n"
            f"TOPO_MED =  {np.median(data[:,2]):.2f}\n",
            file=f
            )

        print(f"\nMESHFEM INTERFACE STATS\n{'='*50}\n"
            f"NPTS_LON = {len(x)}\n"
            f"NPTS_LAT = {len(y)}\n"
            f"DELTA_X  = {deltax}\n"
            f"DELTA_Y  = {deltay}\n"
            f"DELTA_LAT  = {deltalat}\n"
            f"DELTA_LON  = {deltalon}\n"
            f"dX  = {dx}\n"
            f"dY  = {dy}\n"
            f"dLON    = {dlon}\n"
            f"dLAT    = {dlat}\n"
            f"LON_MIN = {lon_min}\n"
            f"LAT_MIN = {lat_min}\n"
            f"LON_MAX = {lon_max}\n"
            f"LAT_MAX = {lat_max}\n",
            file=f
            )

        print(f"Place one of the following in your 'interfaces.dat' file:\n\n",
              file=f)
        
        print(f" .false. {len(x)} {len(y)} {lon_min:.2f}d0 {lat_min:.2f}d0 "
              f"{dlon:.6f}d0 {dlat:.6f}d0", file=f)
        print(f" .true. {len(x)} {len(y)} {x.min():.0f}d0 {y.min():.0f}d0 "
              f"{dx:.2f}d0 {dy:.2f}d0", file=f)
                

def overwrite_interfaces_line(interfaces_fid, data):
    """
    Sort of brute force approach to modify the interfaces file expecting that
    the topography file name is `interface_topo.dat`
    """
    assert os.path.exists(interfaces_fid)

    x, y = data[:, :2].T
    xvals = np.unique(x)
    yvals = np.unique(y)
    dx = xvals[1] - xvals[0]
    dy = yvals[1] - yvals[0]

    print("replacing topography line in interfaces file")

    # Replacement topo line matching the created file. We use the UTM 
    # version because it probably has cleaner numbers
    topo_line = (
        f" .true. {len(xvals)} {len(yvals)} "
        f"{xvals.min():.0f}d0 {yvals.min():.0f}d0 "
        f"{dx:.2f}d0 {dy:.2f}d0\n"
        )

    lines = open(interfaces_fid, "r").readlines()
    for l, line in enumerate(lines):
        if "interface_topo.dat" in line:
            lines[l-1] = topo_line
            break
    with(open(interfaces_fid, "w")) as f:
        f.writelines(lines)


def main(tag, method, srtm_files, x_min, x_max, y_min, y_max, coords, spacing_m,
         border_m=10E3, buffer_m=None, utm_projection=-60, plot=False, 
         path="./"
         ):
    """
    Create topography files for a mesher using underlying SRTM30P data files

    :type tag: str
    :param tag: tag for file saving
    :type method: str
    :param method: file writing method, either 'trelis' or 'meshfem'
        - trelis: writes all points in x, y, z
        - meshfem: only writes z points
    :type srtm_files: list of str
    :param srtm_files: path locations of the SRTM30P files to be used to
        generate the topography files
    :type x_min: float
    :param x_min: minimum x value
    :type x_max: float
    :param x_max: maximum x value
    :type y_min: float
    :param y_min: minimum y value
    :type y_min: float
    :param y_max: maximum y value
    :type coords: str
    :param coords: coordinate system to use for the mesh, either 'latlon' or
        'xyz'. If 'latlon', the coordinates are given in latitude and longitude
        and will be converted to UTM coordinates. If 'xyz', the coordinates are
        given in UTM coordinates (easting, northing) and will not be converted.
    :type spacing_m: int
    :param spacing_m: spacing of the new interpolated mesh
    :type border_m: float
    :param border_m: border to add to the initial cut of the topography file
        (in meters) so that interpolation won't hit the edge of the file
    :type utm_projection: int
    :param utm_projection: projection label, defaults to -60 / 60S for 
        New Zealand
    """
    if coords == "latlon":
        print(f"converting latlon to UTM {utm_projection} coordinates")
        x_min, y_min = lonlat_utm(x_min, y_min, utm_projection)
        x_max, y_max = lonlat_utm(x_max, y_max, utm_projection)

    # Add buffer to the given bounds for flexibility
    if buffer_m is not None:
        x_min -= buffer_m
        y_min -= buffer_m
        x_max += buffer_m
        y_max += buffer_m

    print(f"x_min, x_max: {x_min*1e-3:.2f}, {x_max*1e-3:.2f} km")
    print(f"y_min, y_max: {y_min*1e-3:.2f}, {y_max*1e-3:.2f} km")
    print(f"delta_x ={(x_max - x_min) * 1e-3:.2f}")
    print(f"delta_y ={(y_max - y_min) * 1e-3:.2f}")

    topography = None
    for i, fid in enumerate(srtm_files):
        print(f"reading file {i+1}/{len(srtm_files)}: {os.path.basename(fid)}")
        # Convert from latitude/longitude to UTM 
        print(f"\tconverting latlon to UTM {utm_projection}")
        topo_utm = convert_coordinates(data=read_nc(fid),
                                       utm_projection=utm_projection
                                       )
        print(f"\ttopography file has size {np.shape(topo_utm)}")
        print(f"\txmin,xmax = {topo_utm[:,0].min()}, {topo_utm[:,0].max()}")
        print(f"\tymin,ymax = {topo_utm[:,1].min()}, {topo_utm[:,1].max()}")

        # Cut the topography file to the desired mesh boundaries with some buffr
        print("\tcutting topography file to desired dimensions")
        topo_cut = cut_topography(data=topo_utm,
                                  x_or_lon_min=x_min - border_m,
                                  x_or_lon_max=x_max + border_m,
                                  y_or_lat_min=y_min - border_m,
                                  y_or_lat_max=y_max + border_m
                                  )
        print(f"\ttopography file cut to size {np.shape(topo_cut)}")
        if topography is None:
            topography = topo_cut
        else:
            topography = np.concatenate((topography, topo_cut))

    print(f"final topo file has size {np.shape(topography)} ")

    # Interpolate the topography file to a regular grid
    print("interpolating to regular grid")
    topo_interp = interpolate_points(data=topography,
                                     x_min=x_min, x_max=x_max,
                                     y_min=y_min, y_max=y_max,
                                     spacing_m=spacing_m,
                                     plot=plot
                                     )
    print(f"\tinterpolated file has size {np.shape(topo_interp)}")
    
    # Check if output topography file has NaNs
    if np.any(np.isnan(topo_interp)):
        print(f"ERROR: NaNs found in output topography, check that the input "
              f"nc files cover the bounding box you are cutting and use "
              f"`border_m` to cut a larger bounding box than needed")   
        # vals = topo_interp[np.isnan(topo_interp[:,-1])][:,:2]
        # print(vals)

        sys.exit()

    # Write the topo and moho files
    print(f"writing to {tag}")
    write_data(data=topo_interp, method=method, tag=tag, path=path)

    if method == "meshfem":
        log = os.path.join(path, f"{tag}_log.txt")
        print_meshfem_stats(topo_interp, utm_projection, log=log)

        # Overwrite the interfaces.dat file with the new topography
        interfaces_fid = os.path.join(path, "interfaces.dat")
        overwrite_interfaces_line(interfaces_fid, topo_interp)
        

if __name__ == "__main__":
    try:
        fid = sys.argv[1]
    except IndexError:
        print(f"choice must be in 'NZ', 'AK', 'NK'")
        sys.exit()  

    # Read the configuration file
    pars = set_parameters(fid)

    # Find the topography data files that are to be used
    srtm_files = glob(pars["topography_path"])
    if not srtm_files:
        sys.exit("No input .nc files found")

    tag = os.path.basename(fid).split(".")[0]

    main(tag=tag, method=pars["method"], srtm_files=srtm_files, 
         x_min=pars["x_min"],  x_max=pars["x_max"], 
         y_min=pars["y_min"], y_max=pars["y_max"], 
         coords=pars["coords"], spacing_m=pars["spacing_m"], 
         border_m=pars["border_m"], utm_projection=pars["utm_projection"], 
         plot=pars["plot"], path= pars["path"]
         )
    
