"""
Functionality to plot regular XYZ grid files which define 3D models. These are 
expected to be created as SPECFEM external tomography models, or from other 
utility functions. 

.. note::
    Originally located in Pyatoa visuals directory. Moved here for more rapid
    development. Will be moved back to Pyatoa once it is more complete.
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt


class XYZViewer:
    """
    A class to read, manipulate and plot structured grid files (xyz) that are
    in the format of the Specfem external tomography file.
    """
    def __init__(self, variables=None):
        """
        Initiate the class with some empty variables

        :type variables: list
        :param variables: variable names in the order that the columns are
            defined in the file  
        """
        self.variables = variables or ["x", "y", "z", "vp", "vs", 
                                       "rho", "qp", "qs"]
        self.zero_origin = False
        self.xgrid = None
        self.ygrid = None
        self._coast = None

    def read(self, fid, fmt="specfem"):
        """
        A wrapper for np.loadtxt to read in xyz files.

        .. note::
            For now this function assumes that it is reading an external
            tomography file that is formatted how SPECFEM expects it.

        :type fid: str
        :param fid: file id to read from
        """
        if fmt == "specfem":
            values = np.loadtxt(fid, dtype=float, skiprows=4).T
        elif fmt == "semslicer":
            values = np.loadtxt(fid, dtype=float, delimiter=",").T
        else:
            raise ValueError("fmt must be 'specfem' or 'semslicer'")

        # It is possible that attenuation is not included, if so change the
        # default variables
        assert(len(values) == len(self.variables)), \
                (f"Number of columns {len(values)} does not match number of "
                 f"excepted variables {len(self.variables)}")

        # Set each of the columns as an internal variable
        for i, variable in enumerate(self.variables):
            setattr(self, variable, values[i])

        self.check()

    def read_semslicer(self, fid):
        """
        A wrapper for np.loadtxt to read in xyz files that are outputted by the
        semslicer fortran script

        :type fid: str
        :param fid: file id to read from
        """
        values = np.loadtxt(fid, dtype=float, skiprows=4).T

        # It is possible that attenuation is not included
        if len(values) != len(self.variables):
            variables = self.variables[:6]
            assert(len(values) == len(variables)), \
                f"Number of columns ({len(values)}) does not match known format"

        # Set each of the columns as an internal variable
        for i, variable in enumerate(self.variables):
            setattr(self, variable, values[i])

        self.check()

    def compare(self, fid, choice=""):
        """
        Read in a corresponding xyz file and difference values
        :param fid:
        :return:
        """
        pass

    def show(self, **kwargs):
        """
        A wrapper for pyplot.show() to avoid importing pyplot just for show()
        """
        plt.show(**kwargs)

    def savefig(self, fid):
        """
        A wrapper for pyplot.savefig()
        """
        if not os.path.exists(os.path.dirname(os.path.abspath(fid))):
            os.mkdir(os.path.dirname(fid))

        plt.savefig(fid)

    def close(self, *args, **kwargs):
        """
        A wrapper for pyplot.close('all')
        """
        plt.close(*args, **kwargs)

    def grid(self):
        """
        Define a Numpy mesh grid that can be used for contour plotting, using
        internal variables. Will be used to reshape data arrays
        """
        if self.zero_origin:
            self.xgrid, self.ygrid = np.meshgrid(
                np.arange(0, self.x_max - self.x_min + self.dx, self.dx),
                np.arange(0, self.y_max - self.y_min + self.dy, self.dy)
            )
        else:
            self.xgrid, self.ygrid = np.meshgrid(
                np.arange(self.x_min, self.x_max + self.dx, self.dx),
                np.arange(self.y_min, self.y_max + self.dy, self.dy)
            )
        assert (np.shape(self.xgrid) == (self.ny, self.nx)), "Error in gridding"

        # Create Z-grid for cross-sections (currently NO zero origin option)
        self.xxgrid, self.zxgrid = np.meshgrid(
            np.arange(self.x_min, self.x_max + self.dx, self.dx),
            np.arange(self.z_min, self.z_max + self.dz, self.dz)
        )

        self.yygrid, self.zygrid = np.meshgrid(
            np.arange(self.y_min, self.y_max + self.dy, self.dy),
            np.arange(self.z_min, self.z_max + self.dz, self.dz)
        )

    def coast(self, fid=None):
        """
        Plot the coastline based on a two column text file which defines 
        X and Y coordinates of the coastline which should fall into the 
        coordinate system 

        :type fid: str or None
        :param fid: text file containing coastline data. If not given, assumes
            this function has been run before
        """
        if self._coast is None:
            assert(fid is not None), f"2-column ascii required for coast"
            self._coast = np.loadtxt(fid)

        if self.zero_origin and self._coast[0][0] != 0:
            for i in [0, 1]:
                self._coast[:, i] = self._coast[:, i] - self._coast[:, i].min()

        plt.scatter(self._coast[:, 0], self._coast[:, 1], c="k", s=1, zorder=5)

    def check(self):
        """
        Check the min/max values, grid spacing, npts etc. assign internally
        """
        for variable in self.variables:
            values = getattr(self, variable)
            setattr(self, f"{variable}_min", values.min())
            setattr(self, f"{variable}_max", values.max())
            if variable in ["x", "y", "z"]:
                unique_values = np.unique(values)
                spacing = abs(unique_values[1] - unique_values[0])
                setattr(self, f"unique_{variable}", unique_values)
                setattr(self, f"n{variable}", len(unique_values))
                setattr(self, f"d{variable}", spacing)

        setattr(self, "npts", len(values))
        self.grid()

    def print_vals(self):
        """
        Print out values created from the check() function
        """
        print("\n")
        for variable in self.variables:
            val_min = getattr(self, f"{variable}_min")
            val_max = getattr(self, f"{variable}_max")

            print(f"{val_min} <= {variable} <= {val_max}; ")
            
            if variable in ["x", "y", "z"]:
                unique_values = getattr(self, f"unique_{variable}")
                dval = getattr(self, f"d{variable}")
                nval = getattr(self, f"n{variable}")

                print(f"d{variable}={dval}, n{variable}={nval}")
                print(unique_values)
            print("\n")

    def decimate(self, factor):
        """
        Crudely decimate the values, if, e.g. only a low-resolution plot needs
        to be made, speeds up plotting calls  

        :type factor: int
        :param factor: factor to decimate npts by
        """
        # Two methods for decimating, one if factor is a factor of npts
        if not self.npts % factor:
            for variable in self.variables:
                values = getattr(self, variable)
                values = values.reshape(-1, factor).max(1)
                setattr(self, variable, values)
        else:
            for variable in self.variables:
                values = getattr(self, variable)
                values = np.maximum.reduceat(values,
                                             np.arange(0, self.npts, factor)
                                             )
                setattr(self, variable, values)

        self.check()

    def depth_slice(self, depth_km, variable, **kwargs):
        """
        Plot a slice across the XY plane at a given Z-axis value `depth_km` 
        Kwargs are passed to pyplot.countourf

        :type depth_km: float
        :param depth_km: depth value to set depth slice at. Must match a depth
            value in the XYZ model
        :type variable: str
        :param variable: 
        :return:
        """
        assert(depth_km in self.unique_z), f"{depth_km} is not a valid Z value"
        assert(variable in self.variables), \
            f"{variable} not in {self.variables}"

        idx = np.where(self.z == depth_km)[0]
        value = getattr(self, variable)[idx]

        value = value.reshape(np.shape(self.xgrid))
        assert(idx.any()), "No values found for Z={depth_km}km"

        f, ax = plt.subplots()
        plt.contourf(self.xgrid, self.ygrid, value, **kwargs)
        cbar = plt.colorbar(fraction=0.046)
        cbar.ax.set_ylabel(f"{variable.title()}", rotation=270, labelpad=15)

        ax.set_xlim(ax.get_xlim())
        ax.set_ylim(ax.get_ylim())
        ax.set_aspect("equal")

        plt.xlabel("X [m]")
        plt.ylabel("Y [m]")
        plt.title(f"{variable.title()} at Z={depth_km} km")

        return f, ax

    def x_cross_section(self, xval, variable, vmin=None, vmax=None, **kwargs):
        """
        Plot a cross section across the YZ plane given a unique X-axis value
        """
        values = getattr(self, "unique_x"), f"{xval} is not a valid X value"
        assert(variable in self.variables), \
            f"{variable} not in {self.variables}"

        idx = np.where(self.x == xval)[0]
        value = getattr(self, variable)[idx]

        value = value.reshape(np.shape(self.yygrid))
        assert(idx.any()), "No values found for X={xval}km"

        f, ax = plt.subplots()
        plt.contourf(self.yygrid, self.zygrid, value, **kwargs)
        cbar = plt.colorbar(fraction=0.046)
        cbar.ax.set_ylabel(f"{variable.title()}", rotation=270, labelpad=15)

        ax.set_xlim(ax.get_xlim())
        ax.set_ylim(ax.get_ylim())

        plt.xlabel("Y [m]")
        plt.ylabel("Z [m]")
        plt.title(f"{variable.title()} at X={xval} km")

        return f, ax

    def y_cross_section(self, yval, variable, vmin=None, vmax=None, **kwargs):
        """
        Plot a cross section across the XZ plane given a unique Y-axis value
        """
        values = getattr(self, "unique_y"), f"{yval} is not a valid Y value"
        assert(variable in self.variables), \
            f"{variable} not in {self.variables}"

        idx = np.where(self.y == yval)[0]
        value = getattr(self, variable)[idx]

        value = value.reshape(np.shape(self.xxgrid))
        assert(idx.any()), "No values found for Y={yval}km"

        f, ax = plt.subplots()
        plt.contourf(self.xxgrid, self.zxgrid, value, **kwargs)
        cbar = plt.colorbar(fraction=0.046)
        cbar.ax.set_ylabel(f"{variable.title()}", rotation=270, labelpad=15)

        ax.set_xlim(ax.get_xlim())
        ax.set_ylim(ax.get_ylim())

        plt.xlabel("X [m]")
        plt.ylabel("Z [m]")
        plt.title(f"{variable.title()} at Y={yval} km")

        return f, ax

    def volume(self, variable):
        """
        Plot the 3D volume
        :return:
        """
        raise NotImplementedError

        assert(variable in self.variables), f"{variable} not found"
        f = plt.figure()
        ax = plt.axes(projection="3d")
        ax.plot_surface(self.x, self.y, self.z, color=getattr(self, variable))


if __name__ == "__main__":
    # Example run script to plot chosen depth slices. Call it with e.g.,:
    # $ python xyz_viewer.py tomography_model_8km.xyz

    #               vvvvv Set parameters here vvvvv

    variables = ["x", "y", "z", "vp", "vs", "rho"]  # columns of .xyz file

    # Choose what slices you want to plot
    # If you set them to None, it will plot ALL available slices/depths
    # Otherwise provide a list of matching values. Empty lists will skip
    depths = [0]  # XY plane depth slice with chosen Z-axis depth value
    x_sections = [118139.5]  # YZ plane cross-section w/ chosen X-axis value
    y_sections = [6817222.6]  # XZ plane cross-section w/ chosen Y-axis value
    decimate = None  # integer values if you want to reduce array size

    # Choose how they should be colored
    plot_parameters = ["vs"]  # which 'variables' to plot, e.g., vp, vs, rho
    colormap = "RdYlBu"  # Must match Matplotlib Colormap 
    colorlvl = 256  # Discretization of the colorscale (integer)
    vmin = None  # Set colorscale min value 
    vmax = None  # Set colorscale max value

    # Additional plotting choices
    zero_origin = False  # Set the origin at (0, 0, 0) if not already
    coastline_file = None  # Optional two-column XY file for coastline

    #               ^^^^^ Set parameters here ^^^^^


    # Instantiate XYZViewer class with correct variables
    xyz = XYZViewer(variables=variables)
    file_id = sys.argv[1]
    base_fid = os.path.splitext(file_id)[0]  # for figure naming

    print(f"Reading/plotting XYZ file: {file_id}")
    xyz.read(file_id, fmt="specfem")
    xyz.print_vals()

    # Set internal parameters
    xyz.zero_origin = zero_origin
    
    # Set default parameters if NoneTypes given
    if plot_parameters is None:
        plot_parameters = variables[3:]
    if depths is None:
        depths = xyz.unique_z
    if x_sections is None:
        x_sections = xyz.unique_x
    if y_sections is None:
        y_sections = xyz.unique_y

    if decimate:
        print(f"\tdecimating arrays by {decimate}")
        xyz.decimate(decimate)

    for depth in depths:
        print(f"plotting depth: {depth}km")
        for variable in plot_parameters:
            print(f"\tplotting variable: {variable}")
            xyz.depth_slice(depth, variable, cmap=colormap, levels=colorlvl,
                            vmin=vmin, vmax=vmax)

            if coastline_file:
                xyz.coast(coastline_file)

            savefid = f"{base_fid}_{variable}_z_{depth}m.png"
            xyz.savefig(savefid)
            xyz.close()

    for x_section in x_sections:
        print(f"plotting X-axis cross-section (YZ plane)")
        for variable in plot_parameters:
            print(f"\tplotting variable: {variable}")
            xyz.x_cross_section(x_section, variable, cmap=colormap, 
                               levels=colorlvl, vmin=vmin, vmax=vmax)

            savefid = f"{base_fid}_{variable}_x_{x_section}m.png"
            xyz.savefig(savefid)
            xyz.close()

    for y_section in y_sections:
        print(f"plotting Y-axis cross-section (XZ plane)")
        for variable in plot_parameters:
            print(f"\tplotting variable: {variable}")
            xyz.y_cross_section(y_section, variable, cmap=colormap, 
                              levels=colorlvl, vmin=vmin, vmax=vmax)

            savefid = f"{base_fid}_{variable}_y_{y_section}m.png"
            xyz.savefig(savefid)
            xyz.close()
