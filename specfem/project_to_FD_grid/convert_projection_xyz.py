"""
Small script to read in and plot slices of a GLL model that has been projected
to a regular FD grid using the SPECFEM3D Cartesian auxiliary function 
xproject_and_combine_vol_data_on_regular_grid
"""
import os
import sys
import matplotlib.pyplot as plt
import numpy as np


def read_fortran_binary(filename):
    """
    Reads Fortran-style unformatted binary data into numpy array.

    .. note::
        The FORTRAN runtime system embeds the record boundaries in the data by
        inserting an INTEGER*4 byte count at the beginning and end of each
        unformatted sequential record during an unformatted sequential WRITE.
        see: https://docs.oracle.com/cd/E19957-01/805-4939/6j4m0vnc4/index.html

    :type filename: str
    :param filename: full path to the Fortran unformatted binary file to read
    :rtype: np.array
    :return: numpy array with data with data read in as type Float32
    """
    nbytes = os.path.getsize(filename)
    with open(filename, "rb") as file:
        # read size of record
        file.seek(0)
        n = np.fromfile(file, dtype="int32", count=1)[0]

        if n == nbytes - 8:
            file.seek(4)
            data = np.fromfile(file, dtype="float32")
            return data[:-1]
        else:
            file.seek(0)
            data = np.fromfile(file, dtype="float32")
            return data


def create_grid(fd_proj_grid): 
    """
    Read the single array projected file and convert into a NumPy array
    """
    # Get the defining dimensions of the FD projection grid
    # o? = origin, h? = sampling rate, n? = npoints
    ox, oy, oz, hx, hy, hz, nx, ny, nz = np.loadtxt(fd_proj_grid).ravel()

    # n? values are used for indexing, need to be ints
    nx = nx.astype(int)
    ny = ny.astype(int)
    nz = nz.astype(int)
    
    # Restructure data array for plotting
    x = np.linspace(ox, ox + (hx * (nx - 1)), nx)
    y = np.linspace(oy, oy + (hy * (ny - 1)), ny)
    z = np.linspace(oz, oz + (hz * (nz - 1)), nz)

    arr = []
    for z_ in z:
        for y_ in y:
            for x_ in x:
                arr.append([x_, y_, z_])

    xyz = np.array(arr)

    return xyz


if __name__ == "__main__":
    choice = "smooth" 
    fd_proj_grid = "fd_proj_grid.txt"

    if choice == "smooth":
        # fid_template = "OUTPUT_FILES/projected_smooth_3x1km/{}_projected.bin"
        # fid_template = "OUTPUT_FILES/projected_smooth_4x2km/{}_projected.bin"
        fid_template = "OUTPUT_FILES/projected_smooth_500x500m/{}_projected.bin"
    elif choice == "raw":
        fid_template = "OUTPUT_FILES/projected_raw/{}_projected.bin"

    # Read Vp first to initialize the array
    arr = read_fortran_binary(fid_template.format("vp"))

    # Then read the remaining parameters and stack alongside
    for fid in ["vs", "rho", "qmu", "qkappa"]:
        arr = np.vstack((arr, read_fortran_binary(fid_template.format(fid))))

    xyz = create_grid(fd_proj_grid)

    # Combine grid with data
    data = np.hstack((xyz, arr.T))

    np.savetxt(f"{choice}_shallow_tomo_file.xyz", data, fmt="%6.3f")

