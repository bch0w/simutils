"""
Script to read, convert, and format the NZ Wide 2.X velocity model of 
Eberhart-Phillips et al. (2021) into a .xyz file format readable by 
SPECFEM3D Cartesian.
Partially based off Carl Tape's Matlab scripts (GEOTOOLS)
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.ndimage.filters import gaussian_filter
from pyproj import Proj


class Data(dict):
    """
    Dictionary that can be accessed using variables
    """
    def __setattr__(self, key, value):
        self[key] = value

    def __getattr__(self, key):
        return self[key]


def myround(x, base=5, choice='near'):
    """
    Round value x to nearest base, round 'up','down' or to 'near'est base

    :type x: float
    :param x: value to be rounded
    :type base: int
    :param base: nearest integer to be rounded to
    :type choice: str
    :param choice: method of rounding, 'up', 'down' or 'near'
    :rtype roundout: int
    :return: rounded value
    """
    if choice == 'near':
        roundout = int(base * round(float(x)/base))
    elif choice == 'down':
        roundout = int(base * np.floor(float(x)/base))
    elif choice == 'up':
        roundout = int(base * np.ceil(float(x)/base))

    return roundout


class Converter():
    """
    Collection of functions that are used to convert from NZwide velocity model
    to the SPECFEM3D Cartesian .xyz grid
    """
    def __init__(self, fid, fid_qp=None, fid_qs=None):
        """
        Start by reading in the data files
        """
        data = self.read(fid, fid_qp, fid_qs)
        self.data = data

    def debug(self):
        """
        For internal debug use
        """
        from IPython import embed
        import ipdb
        ipdb.set_trace()
        embed(colors="neutral")


    def read(self, fid, fid_qp=None, fid_qs=None):
        """
        Read the .txt file provided by Eberhart-Phillips which contains the 
        velocity file. Optional attach Qp and Qs data points to the velocity model

        Header on velocity file is assumed to be:
            vp, vp/vs, vs, rho, sf_vp, sf_vpvs, x, y, depth, lat, lon
        """
        print("reading data")
        # Load in the main velocity model file
        data = np.loadtxt(fid, skiprows=2).T
        vp, vpvs, vs, rho, sf_vp, sf_vpvs, x, y, z, lat, lon = data

        dict_out = Data(vp=vp, vpvs=vpvs, vs=vs, rho=rho, sf_vp=sf_vp, 
                        sf_vpvs=sf_vpvs, x=x, y=y, z=z, lat=lat, lon=lon)

        # Pretty print bounds of each value
        spacing = max([len(_) for _ in dict_out.keys()])
        for key, val in dict_out.items():
            print(f"\t{key:<{spacing}}: {val.min()} to {val.max()}")
        print(f"\t{'depths':<{spacing}}: {np.unique(dict_out['z'])}")

        # Optional addition of attenuation data (Qp, Qs)
        if fid_qp:
            qp = np.loadtxt(fid_qp, skiprows=2).T[0]
            assert(len(qp) == len(vp))
            dict_out["qp"] = qp
        if fid_qs:
            qs = np.loadtxt(fid_qs, skiprows=2).T[0]
            assert(len(qs) == len(vp))
            dict_out["qs"] = qs

        return dict_out

    def convert_data(self):
        """
        Convert value formats to correct units and directions.
        """
        print("converting data units")
        self.data.vp *= 1E3  # km/s -> m/s
        self.data.vs *= 1E3  # km/s -> m/s
        self.data.x *= 1E3   # km -> m
        self.data.y *= 1E3   # km -> m
        self.data.z *= -1E3   # km -> m and LH coord -> RH coord

    def cut(self, lat_min=None, lat_max=None, lon_min=None, lon_max=None, 
            z_min=None, z_max=None):
        """
        Cut the original mesh dimensions down in size to reduce computation
        Assuming that the data has been converted already
        
        south_island: -47.5,-40 and 165,176
        """
        print("cutting coordinates to desired domain")
        idx = np.where((self.data.lat > lat_min) & (self.data.lat < lat_max) & 
                       (self.data.lon > lon_min) & (self.data.lon < lon_max))[0]
        print(f"{len(self.data.lat)} -> {len(idx)} data points after cut")
        data = {key: val[idx] for key, val in self.data.items()} 
        self.data = Data(data)

    def convert_coords(self, utm_zone=-60):
        """
        Convert form WGS84 Lat Lon to UTM zone for easier .xyz grid spacing
        Uses PyProj, taken from Pyatoa
        """
        print("converting to UTM zone")

        # Format UTM zone because the projection doesn't want negative values
        south = bool(utm_zone < 0)

        # Define the projection from coordinate system A to B
        projection = Proj(proj="utm", zone=abs(utm_zone), south=south, 
                          ellps="WGS84", preserve_units=False)
        x, y = projection(self.data.lon, self.data.lat, inverse=False)

        # Overwriting the standard x and y set in the original model
        self.data.x = x
        self.data.y = y

    def regularize(self, dx=1E3, dy=1E3, method="linear", sigma=0.7):
        """
        Convert the irregular node spacing of NZwide into a regular grid
        """
        print(f"regularizing onto uniform grid")
        for z in sorted(np.unique(self.data.z))[::-1]:
            idx = np.where(self.data.z == z)[0]

            # Rounding the values for clean UTM bounds
            x_min = myround(self.data.x.min(), 100, "down")
            x_max = myround(self.data.x.max(), 100, "up")
            y_min = myround(self.data.y.min(), 100, "down")
            y_max = myround(self.data.y.max(), 100, "up")

            xi = np.linspace(x_min, x_max, int(dx))
            yi = np.linspace(y_min, y_max, int(dy))
            xi, yi = np.meshgrid(xi, yi)

            # Interpolate unstructured grid using scipy
            for key in ["vs", "vp", "rho", "qp", "qs"]:
                zi = griddata((self.data.x, self.data.y), self.data[key], 
                              (xi, yi), method=method)

                # Smooth the interpolated array
                print(f"smoothing {key} with gaussian sigma={sigma}")
                zi = gaussian_filter(zi, sigma=sigma)

                self.plot_contour(xi, yi, zi, self.data[key], key)
                a=1/0

    def nearest_node(self, x0, y0):
        """
        Smoothing will depend on the nearest node. More smoothing for isolated
        nodes. This function simply finds the nearest node value.
        """



    def plot_coast(self):
        """
        Plot NZ Coastline
        """
        coast = ("/Users/Chow/Documents/academic/vuw/data/carto/coastline/"
                 "nz_coast_full.txt")
        x, y = np.loadtxt(coast).T
        plt.scatter(x/1e3, y/1e3, s=.01, c="k", zorder=20)

    def plot_contour(self, xi, yi, zi, z, key):
        """
        Plot the contour of the new gridded data
        """
        plt.scatter(self.data.x, self.data.y, c=self.data[key], s=.5, 
                    marker="o", zorder=10)
        plt.contourf(xi, yi, zi, zorder=5, levels=64)
        self.plot_coast()

        plt.title(f"{z} m depth, {key}")
        plt.xlabel("X [m]")
        plt.ylabel("Y [m]")
        plt.colorbar(label=key)
        plt.show()

    def plot_xyz(self, coord="xy", choice="vs"):
        """
        Plot the nodal points from the original velocity model to look at 
        data coverage
        """
        unique_z = np.unique(self.data.z)
        idx = np.where(self.data.z == unique_z.max())[0]

        if coord == "xy":
            x = self.data.x[idx]
            y = self.data.y[idx]
        elif coord == "latlon":
            x = self.data.lon[idx]
            x[x<0] += 360  # wrap lat lon so that it's not -180 to 180
            y = self.data.lat[idx]

        val = self.data[choice][idx]
         
        plt.scatter(x/1e3, y/1e3, marker="o", c="None", s=2, edgecolor="k", 
                    linewidth=1)
        if coord == "xy":
            plt.xlabel("X [km]")
            plt.ylabel("Y [km]")
        elif coord == "latlon":
            plt.xlabel("Longitude")
            plt.ylabel("Latitude")

        plt.title(f"NZWide2.2 ({choice}) at Z={unique_z.max()/1E3}km")
        plt.colorbar(label=choice)
        plt.gca().set_aspect(1)
        self.plot_coast()
        self.set_axes(plt.gca())

        plt.show()

    def set_axes(self, ax):
        """
        Give axis a default look
        """
        ax.title.set_fontsize(15)
        ax.xaxis.label.set_fontsize(14)
        ax.yaxis.label.set_fontsize(14)
        ax.tick_params(axis="both", which="both", width=1.5,
                       direction="in", labelsize=12, length=5)
        for axis in ["top", "bottom", "left", "right"]:
            ax.spines[axis].set_linewidth(2)



if __name__ == "__main__":
    # File inputs
    base_dir = "/Users/Chow/Documents/academic/vuw/data/nzwide"
    vel_fid =  os.path.join(base_dir, "vlnzw2p2dnxyzltln.tbl.txt")
    qp_fid = os.path.join(base_dir, "Qpnzw2p2xyzltln.tbl.txt")
    qs_fid = os.path.join(base_dir, "Qsnzw2p2xyzltln.tbl.txt")

    # Parameter set
    interp_method = "cubic"
    lat_min = -47.75
    lat_max = -39.75
    lon_min = 164.75
    lon_max = 176.25
    z_min = None
    z_max = 400000
    utm_zone = -59
    sigma = 0.7

    # Run 
    conv = Converter(vel_fid, qp_fid, qs_fid)
    conv.convert_data()
    conv.convert_coords(utm_zone=utm_zone)
    conv.regularize(method=interp_method, sigma=sigma)

    conv.plot_xyz(coord="xy", choice="sf_vp")
    # conv.cut(lat_min, lat_max, lon_min, lon_max, z_min, z_max)
