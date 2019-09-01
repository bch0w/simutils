
"""
There are some values required by the internal mesher that control the sizes
of spectral elements in the horizontal and vertical directions. This collection
of functions should help determine those values in a standardized manner
"""
import numpy as np
from scipy.interpolate import griddata
from mesh_utils import lonlat_utm, myround


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


def cut_topography(data, x_or_lon_min=150000., x_or_lon_max=750000.,
                   y_or_lat_min=5000000., y_or_lat_max=6000000.):
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


def interpolate_points(data, x_min=174000., x_max=680000., y_min=5200000.,
                       y_max=5950000., spacing_m=1000):
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
    x_irreg = data[:, 0]
    y_irreg = data[:, 1]
    topography = data[:, 2]

    # Create the regular grid to be interpolated onto
    x_reg = np.arange(x_min, x_max, spacing_m)
    y_reg = np.arange(y_min, y_max, spacing_m)

    x_grid, y_grid = np.meshgrid(x_reg, y_reg)
    x_out = x_grid.flatten()
    y_out = y_grid.flatten()

    # Interpolate the topography data
    interp_vals = griddata(
        points=(x_irreg, y_irreg), values=topography, xi=(x_grid, y_grid)
    )

    # Create the ndarray by creating column vectors and mushing em together
    z_out = interp_vals.flatten()

    # Arrange the data in the desired format
    data_out = np.column_stack((x_out, y_out, z_out))

    return data_out


def write_data(data, method, moho=-33000., topo_fid="topo", moho_fid="moho"):
    """
    Saves the topo as an npy and an ascii, with a descriptive file name

    :param data:
    :return:
    """
    # Meshfem only requires a single value of Z-values
    if method == "meshfem":
        np.savetxt(topo_fid + '.dat', data[:, 2], fmt='%d')
    # Trelis requires XYZ files. If a flat moho is required, it must be the same
    # sampling rate as the topo
    elif method == "trelis":
        np.savetxt(topo_fid + '.xyz', data, fmt='%.1f')
        if moho:
            data_moho = np.copy(data)
            data_moho[:, 2] = moho
            np.savetxt(moho_fid + '.xyz', data_moho, fmt='%.1f')
    else:
        print("Invalid method")


def create_topo_nz_north_utm():
    """
    Create topo and moho files for a custom NZ North mesh

    NZ_NORTH_68
    I set the mesh bounds to be LLC = -43,173 to URC = -37,179
    this translates to UTM -60 LLC = -173e3, 5231e3 to URC = 677e3, 5903e3

    But I want a simulation buffer to prevent any mesh edge effects such as
    spurious reflections, so I put a ~50km buffer to each mesh side, and then
    round to a nice round number in UTM -60. This corresponds to

    x_min = 125E3
    x_max = 725E3
    y_min = 5150E3
    y_max = 5950E3
    moho = -33E3
    spacing_m = 1E3

    """
    # Set the topography parameters here
    x_min = 125E3
    x_max = 725E3
    y_min = 5150E3
    y_max = 5950E3
    moho = -100E3
    spacing_m = 4E3
    utm_projection=-60

    # Load the topography file to be interpolated
    topography_file = (
        "/seis/prj/fwi/bchow/data/mapping/topography/srtm30_e140s10.npy"
    )
    topo_latlon = np.load(topography_file)

    # Create the output file id based on the parameters
    x_length = int((x_max - x_min) * 1E-3)
    y_length = int((y_max - y_min) * 1E-3)
    space_km = int(spacing_m * 1E-3)
    moho_write = int(abs(moho * 1E-3))

    topo_fid = f"topo_utm60_x{x_length}_y{y_length}_{space_km}km"
    moho_fid = f"moho{moho_write}_utm60_x{x_length}_y{y_length}_{space_km}km"

    # Convert from latitude/longitude to UTM -60/ UTM 60south
    print(f"converting latlon to UTM {utm_projection}")
    topo_utm60s = convert_coordinates(data=topo_latlon,
                                      utm_projection=utm_projection
                                      )

    # Cut the topography file to the desired mesh boundaries with some buffer
    buffer = 5E3
    print("cutting topography file to desired dimensions")
    topo_cut = cut_topography(data=topo_utm60s,
                              x_or_lon_min=x_min - buffer,
                              x_or_lon_max=x_max + buffer,
                              y_or_lat_min=y_min - buffer,
                              y_or_lat_max=y_max + buffer
                              )

    # Interpolate the topography file to a regular grid
    print("interpolating to regular grid")
    topo_interp = interpolate_points(data=topo_cut, x_min=x_min, x_max=x_max,
                                     y_min=y_min, y_max=y_max,
                                     spacing_m=spacing_m
                                     )

    # Write the topo and moho files
    print("writing")
    write_data(data=topo_interp, method="trelis", moho=moho, topo_fid=topo_fid,
               moho_fid=moho_fid)


if __name__ == "__main__":
    create_topo_nz_north_utm()
