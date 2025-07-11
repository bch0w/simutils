"""
SPECFEM3D can output movie files as .xyz files which list LAT LON VAL
We can plot these as frames of a movie and collect them later as a .gif
Can also add text, coastlines etc. easily with Matplotlib
"""
import sys
import os
from glob import glob
from subprocess import run
import numpy as np
from scipy import interpolate
import matplotlib as mpl
import matplotlib.pyplot as plt
from PIL import Image
from pyproj import Proj


# COORDINATE CONVERSION CONSTANTS
UTM_ZONE = 52
XMIN = 0
XMAX = 100E3
YMIN = 0
YMAX = 100E3



def read(fid):
    """
    Simple read and parse of the XYZ data files using numpy
    :return: lon, lat, z
    """
    data = np.loadtxt(fid)
    assert(data.shape[1] == 3), ".xyz file is in the wrong format"
    x, y, z = data.T
    return x, y, z


def find(path="./", ext=".xyz"):
    """
    Find all the availalbe data files with the given extension in the given
    path and return an ordered list of all the frames and let the user know
    """
    wildcard = os.path.join(path, "*" + ext)
    files = sorted(glob(wildcard))
    print(f"{len(files)} files with matching extension '{ext}' found")
    print(f"{files[0]} ... {files[-1]}")
    return files


def time_step(fid, dt=1):
    """
    Parse the standard Specfem movie name into a time stamp
    Expected format is e.g., 'gmt_movie_012600.xyz'
    """
    parts = os.path.basename(fid).split("_")
    step = parts[2].split(".")[0]
    try:
        step = int(step)
    except ValueError:
        print("file name does not adhere to expected format")
    return float(step * dt)


def convert_coords(lon, lat, utm_zone=None):
    """
    Convert from the native Lat/Lon coordinate system into UTM 60S coordinates
    This only needs to be done for the first file because the remaining files
    are assumed to follow the same coordinate points.
    """
    if utm_zone is None:
        utm_zone = UTM_ZONE

    if utm_zone < 0:
        south=True
    elif utm_zone > 0:
        south = False

    projection = Proj(proj="utm", zone=abs(utm_zone), south=south, 
                      ellps="WGS84", preserve_units=False)
    x, y = projection(lon, lat, inverse=False)

    # Zero out the origin
    x -= XMIN
    y -= YMIN

    # Convert units of m to km
    x /= 1E3
    y /= 1E3

    return x, y


def plot(ax, x, y, z, min_val=None, max_val=None, show=True, save=False, 
         **kwargs):
    """
    Plot the xyz file in the same fashion each time. Use tricontourf because
    the data is not in a uniform grid so we use a triangular interpolation
    """
    # Mask out values before plotting
    if min_val:
        tri = mpl.tri.Triangulation(x, y)
        masked_vals = np.less(z, min_val)
        mask = np.all(np.where(masked_vals[tri.triangles], True, False), axis=1)
        tri.set_mask(mask)
        # Set custom levels to keep the colorbar segmented the same across figs
        kwargs["levels"] = np.linspace(0, max_val, kwargs["levels"])
        try:
            tcf = ax.tricontourf(tri, z, vmin=0, vmax=max_val, extend="max", 
                                 **kwargs)
        except ValueError:
            # ValueError occurs when we mask out all the values (e.g. T < 0s)
            tcf = ax.tricontourf(x, y, np.zeros(len(z)), vmin=0, vmax=max_val,
                                 extend="max", **kwargs) 
    # Or just plot the values straight up
    else:
        tcf = ax.tricontourf(x, y, z, **kwargs) 

    # Colorbar
    cbar = plt.colorbar(tcf, label=cbar_title, shrink=0.5,aspect=8, 
                        pad=.03, ticks=[0, max_val/2, max_val])
    cbar.ax.yaxis.set_offset_position("left")
    # cbar.ax.set_yticklabels(["0", "2", ">4"])
    # Mark where the min value threshold is set
    cbar.ax.plot([0, 1], [min_val, min_val], "cyan")

    cbar.update_ticks()

    # Accoutrement
    for axis in ["top", "bottom", "left", "right"]:
        ax.spines[axis].set_linewidth(2)


def srcrcv(source=None, receiver=None, convert=False):
    """
    Plot the source and receiver as simple markers
    """
    if source:
        if convert:
            x, y = convert_coords(source[0], source[1])
        else:
            x, y = source
        plt.scatter(x, y, marker="*", s=100, color="k", linewidth=1, 
                    edgecolors="k")
    if receiver:
        if convert:
            x, y = convert_coords(receiver[0], receiver[1])
        else:
            x, y = receiver
        plt.scatter(x, y, marker="v", s=100, color="w", linewidth=1.5, 
                    edgecolors="k")


def gif(path, duration, fid_out="output.gif"):
    """
    Generate a .gif file from all the resulting .png files
    PIL runs into a '[Errno 24] Too many open files' for large numbers of files
    So use ImageMagick convert wrapped with subprocess instead for those
    """
    files = sorted(glob(os.path.join(path, "*.png")))
    if len(files) < 256:
        img, *imgs = [Image.open(f) for f in files] 
        img.save(fp=fid_out, format="GIF", append_images=imgs, save_all=True,
                 duration=duration, loop=1)
    else:
        os.chdir(path)
        call = f"convert -delay 0 -loop 1 *.png {fid_out}"
        run(call.split(" "))


if __name__ == "__main__":
    try:
        input_path = sys.argv[1]
    except IndexError:
        input_path = "./"

    # =========================================================================
    # ACTIONS
    make_pngs = 1
    make_gif = 1
    trial_run = 0
    # =========================================================================
    # PARAMETER SET HERE
    output_path = os.path.join(input_path, "frames")
    gif_fid = os.path.join(input_path, "sim_mov.gif")
    file_ext = ".xyz"
    min_val = 9e-7
    dt = 0.005
    normalized = True
    gif_duration_ms = 200  # milliseconds
    convert = False

    max_val = 1
    source = [25E3, 25E3]
    receiver = None
    text = ""
    # =========================================================================
    # Controls on colorbar
    if normalized:
        kwargs = {"cmap": "gist_ncar_r", # "gist_heat_r",
                  "norm": plt.Normalize(0, max_val), 
                  "levels": 101,
                  }
        cbar_title = "norm of velocity [m/s]"
    else:
        # For single component 
        kwargs = {"cmap": "seismic", 
                  "norm": plt.Normalize(-1 * max_val, max_val), 
                  "levels": 100,
                  }
        cbar_title = "z comp. velocity [m/s]"

    files = []
    # Test files to sample random data points to get an idea of relative amps
    if trial_run:
        files = find(input_path, file_ext)
        n = int(len(files))
        files = [files[0], files[40], files[int(n/3)], files[int(2*n/3)], 
                 files[-1]]

    files = files[:1000]
    # =========================================================================

    # Prep the file system
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    if not files:
        files = find(input_path, file_ext)
        
    # Make the .png files
    xy = None
    if convert:
        x_utm, y_utm = True, True
    else:
        x_utm, y_utm = None, None
    if make_pngs:
        ts_max = time_step(files[-1], dt=dt)
        for i, file_ in enumerate(files):
            # Set up the plot
            ts = time_step(file_, dt=dt)
            x, y, z = read(file_)
            if x_utm is not None:
                if i == 0:
                    x_utm, y_utm = convert_coords(x, y)
                x, y = x_utm, y_utm

            f, ax = plt.subplots(1)
            ax.set_aspect(1)
            print(f"{int(ts):0>3}/{int(ts_max)}")

            # Plot that ish
            plot(ax, x, y, z, min_val=min_val, max_val=max_val,  **kwargs)
            srcrcv(source, receiver, convert=convert)

            plt.title(f"t={ts:6.2f} s")
            plt.xlim([x.min(), x.max()])
            plt.ylim([y.min(), y.max()])
            # plt.xlim([173, 178.5])
            # plt.ylim([-42.5, -37])

            # Clean up the end
            fid_out = os.path.join(output_path, 
                                   os.path.basename(file_) + ".png")
            f.tight_layout()
            plt.savefig(fid_out)
            plt.close()

    # Create the output .gif
    if make_gif:
        gif(output_path, gif_duration_ms)

