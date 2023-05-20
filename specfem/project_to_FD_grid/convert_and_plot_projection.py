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


def read_and_convert_projection(filename, fd_proj_grid): 
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
    
    # Read in the actual data array 
    data = read_fortran_binary(filename)

    # Restructure data array for plotting
    # x_ = np.linspace(ox, ox + (hx * nx), nx)
    # y_ = np.linspace(oy, oy + (hy * ny), ny)
    # z_ = np.linspace(oz, oz + (hz * nz), nz)

    # x, y, z = np.meshgrid(x_, y_, z_, indexing="ij")

    # Restructure data array for plotting
    x = np.linspace(ox, ox + (hx * nx), nx)
    y = np.linspace(oy, oy + (hy * ny), ny)
    z = np.linspace(oz, oz + (hz * nz), nz)

    i = 0
    arr = []
    for z_ in z:
        for y_ in y:
            for x_ in x:
                arr.append([x_, y_, z_, data[i]])
                i += 1

    arr = np.array(arr)

    return arr


def plot3d_projection(x, y, z, data):
    """Plot the 3D visualization of the data"""
    f = plt.figure()
    ax = f.add_subplot(111, projection="3d")

    kw = {
        'vmin': data.min(),
        'vmax': data.max(),
        'levels': np.linspace(data.min(), data.max(), 10),
    }

    _ = ax.contourf(
                x[:, :, 0], y[:, :, 0], data[:, :, 0],
                    zdir='z', offset=0, **kw
                    )
    _ = ax.contourf(
                y[0, :, :], data[0, :, :], z[0, :, :],
                    zdir='y', offset=0, **kw
                    )
    C = ax.contourf(
                data[:, -1, :], y[:, -1, :], z[:, -1, :],
                    zdir='x', offset=x.max(), **kw
                    )

    xmin, xmax = x.min(), x.max()
    ymin, ymax = y.min(), y.max()
    zmin, zmax = z.min(), z.max()
    ax.set(xlim=[xmin, xmax], ylim=[ymin, ymax], zlim=[zmin, zmax])

    # Set labels and zticks
    ax.set(
        xlabel='X [km]',
        ylabel='Y [km]',
        zlabel='Z [m]',
    )

    # Set zoom and angle view
    ax.view_init(40, -30, 0)
    ax.set_box_aspect(None, zoom=0.9)


    f.colorbar(C, ax=ax, fraction=0.02, pad=0.1, label='Vs [m/s]')

    plt.show()


def _find_nearest(array, value):
    """find the nearest value in a given array"""
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


def plot2d_xsection(arr, xval=350000.):
    """2D depth slices at a given Z value""" 
    xval_true = _find_nearest(arr[:,0], xval)
    print(f"xval: {xval} -> {xval_true}")

    dslice = arr[np.where(arr[:,0] == xval_true)]
    yzd = np.delete(dslice, 0, axis=1)

    x_ = yzd[:, 0]
    y_ = yzd[:, 1]
    z_ = yzd[:, 2]
    print(len(x_))
    sc = plt.scatter(x_, y_, c=z_, s=1, cmap="rainbow")
    plt.colorbar(sc)
    plt.grid()

    # test = np.reshape(z_, (100, 400))
    # plt.imshow(test)
    # plt.gca().invert_yaxis()
    
    plt.show()


def plot2d_depth_slice(arr, depth_z=-20000.):
    """2D depth slices at a given Z value""" 
    depth_z_true = _find_nearest(arr[:,2], depth_z)
    print(f"depth_z: {depth_z} -> {depth_z_true}")

    dslice = arr[np.where(arr[:,2] == depth_z_true)]
    xyd = np.delete(dslice, 2, axis=1)

    x_ = xyd[:, 0]
    y_ = xyd[:, 1]
    z_ = xyd[:, 2]
    sc = plt.scatter(x_, y_, c=z_)
    plt.colorbar(sc)
    

    # z_ = np.reshape(xyd[:, 2], (len(x), len(y)))
    # plt.contourf(x_, y_, z_)

    plt.show()
    

if __name__ == "__main__":
    filename = "OUTPUT_FILES/vs_projected.bin"
    fd_proj_grid = "fd_proj_grid.txt"
   
    arr = read_and_convert_projection(filename, fd_proj_grid) 
    
    plot2d_xsection(arr)
    # plot2d_depth_slice(arr)
    # plot3d_projection(x, y, z, data)
