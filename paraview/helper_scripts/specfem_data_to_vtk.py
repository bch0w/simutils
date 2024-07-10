"""
Some functions removed from Pyatoa.utils.read and write, might be useful 
somewhere else later if making a small package to deal with converting vtk files
"""

def read_specfem_vtk(path_to_vtk):
    """

    Read the unstructured grid VTlK files that are output by Specfem for model,
    gradient, and kernel visualizations. Returns a header as a dictionary, and
    the lines of the data file. Useful for manipulating VTK files in place.

    :type path_to_vtk: str
    :param path_to_vtk: full path to the .vtk file to read
    :rtype: tuple (list, dict)
    :return: first index of tuple is a list read line by line from readlines(),
        second index of tuple is a dictionary with all relevant header 
        information from an unstructured_grid vtk file
    """
    with open(path_to_vtk, "r") as f:
        lines = f.readlines()

    # determine important line numbers, headers are specific to specfem3d output
    for i, line in enumerate(lines):
        if "POINTS" in line:
            points_n = int(line.split()[1])
            points_line = i + 1
        elif "CELL" in line and "CELL_TYPE" not in line:
            cells_n = int(line.strip().split()[1])
            cells_size = int(line.strip().split()[2])
            cells_line = i
        elif "CELL_TYPES" in line:
            cell_types_n = int(line.strip().split()[1])
            cell_types_line = i
        elif "POINT_DATA" in line:
            point_data_n = int(line.split()[1])
            point_data_line = i
        elif "SCALARS" in line:
            scalars = line.strip().split()[1]
            data_line = i + 2

    # easier returns in a dictionary
    header_dict = {"points_n": points_n, "points_line": points_line,
                   "cells_n": cells_n, "cells_size": cells_size,
                   "cells_line": cells_line, "cell_types_n": cell_types_n,
                   "cell_types_line": cell_types_line,
                   "points_data_n": point_data_n,
                   "points_data_line": point_data_line, "scalars": scalars,
                   "data_line": data_line
                   }

    return lines, header_dict


def rcv_vtk_from_specfem(path_to_data, path_out="./", utm_zone=-60, z=3E3):
    """
    Creates source and receiver VTK files based on the STATIONS and
    CMTSOLUTIONS from a Specfem3D DATA directory.

    :type path_to_data: str
    :param path_to_data: path to specfem3D/DATA directory
    :type path_out: str
    :param path_out: path to save the fiels to
    :type utm_zone: int
    :param utm_zone: utm zone for converting lat lon coordinates
    :type z: float
    :param z: elevation to put stations at
    """
    from pyatoa.utils.srcrcv import lonlat_utm

    # Templates for filling
    vtk_line = "{x:18.6E}{y:18.6E}{z:18.6E}\n"
    vtk_header = ("# vtk DataFile Version 2.0\n" 
                  "Source and Receiver VTK file from Pyatoa\n"
                  "ASCII\n"
                  "DATASET POLYDATA\n"
                  "POINTS\t{} float\n")

    stations = np.loadtxt(os.path.join(path_to_data, "STATIONS"),
                          usecols=[2, 3], dtype=str)
    lats = stations[:, 0]
    lons = stations[:, 1]

    with open(os.path.join(path_out, "rcvs.vtk"), "w") as f:
        f.write(vtk_header.format(len(stations)))
        for lat, lon in zip(lats, lons):
            rx, ry = lonlat_utm(lon_or_x=lon, lat_or_y=lat, utm_zone=utm_zone,
                                inverse=False)
            f.write(vtk_line.format(x=rx, y=ry, z=z))
        f.write("\n")


def src_vtk_from_specfem(path_to_data, path_out="./", utm_zone=-60, cx=None,
                         cy=None, cz=False):
    """
    Creates source and receiver VTK files based on the STATIONS and
    CMTSOLUTIONS from a Specfem3D DATA directory.

    :type path_to_data: str
    :param path_to_data: path to specfem3D/DATA directory
    :type path_out: str
    :param path_out: path to save the fiels to
    :type utm_zone: int
    :param utm_zone: utm zone for converting lat lon coordinates
    :type cx: float
    :param cx: Constant X-value for creating an Y-slice, should be in units of
        meters, in the UTM coordinate system
    :type cy: float
    :param cy: Constant Y-value for creating an X-slice, should be in units of
        meters, in the UTM coordinate system
    :type cz: float
    :param cz: Constant Z-value for creating a Z-slice, should be in units of
        meters with positve z-axis so negative Z-values are deeper
    """
    from pyatoa.utils.srcrcv import lonlat_utm

    def quick_read_cmtsolution(path):
        """utility function to read cmtsolution file into dictionary object"""
        dict_out = {}
        cmt = np.loadtxt(path, skiprows=1, dtype="str", delimiter=":")
        for arr in cmt:
            # Replace spaces with underscores in the key
            key, value = arr[0].replace(" ", "_"), arr[1].strip()
            # Most values will be float except event name
            try:
                value = float(value)
            except ValueError:
                pass
            dict_out[key] = value
        return dict_out

    # Templates for filling
    vtk_line = "{x:18.6E}{y:18.6E}{z:18.6E}\n"
    vtk_header = ("# vtk DataFile Version 2.0\n" 
                  "Source and Receiver VTK file from Pyatoa\n"
                  "ASCII\n"
                  "DATASET POLYDATA\n"
                  "POINTS\t{} float\n")

    # Gather all the sources
    sources = glob.glob(os.path.join(path_to_data, "CMTSOLUTION*"))

    # Open files that need to be written
    f_std = open(os.path.join(path_out, "srcs.vtk"), "w")
    f_xslice = f_yslice = f_zslice = None
    # Constant X-value means a slice parallel to the Y-Axis, a bit confusing
    if cx:
        f_yslice = open(os.path.join(path_out, f"srcs_yslice_{cx}.vtk"), "w")
    if cy:
        f_xslice = open(os.path.join(path_out, f"srcs_xslice_{cy}.vtk"), "w")
    if cz:
        f_zslice = open(os.path.join(path_out, f"srcs_zslice_{abs(cz)}.vtk"),
                        "w")

    # Write in the headers, use a try-except to write all even if None
    for f in [f_std, f_xslice, f_yslice, f_zslice]:
        try:
            f.write(vtk_header.format(len(sources)))
        except AttributeError:
            continue

    # Create VTK file for Sources, assuming CMTSOLUTION format
    for source in sources:
        src = quick_read_cmtsolution(source)
        sx, sy = lonlat_utm(lon_or_x=src["longitude"], lat_or_y=src["latitude"],
                            utm_zone=utm_zone, inverse=False)
        sz = src["depth"] * -1E3
        # Write data to all files using try-except
        f_std.write(vtk_line.format(x=sx, y=sy, z=sz))
        if cy:
            f_xslice.write(vtk_line.format(x=sx, y=cy, z=sz))
        if cx:
            f_yslice.write(vtk_line.format(x=cx, y=sy, z=sz))
        if cz:
            f_zslice.write(vtk_line.format(x=sx, y=sy, z=cz))

    # Close all the files
    for f in [f_std, f_xslice, f_yslice, f_zslice]:
        try:
            f.write("\n")
            f.close()
        except AttributeError:
            continue



