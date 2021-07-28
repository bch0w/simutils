"""
SPECFEM3D can output movie files as .xyz files which list LAT LON VAL
We can plot these as frames of a movie and collect them later as a .gif
Can also add text, coastlines etc. easily with Matplotlib
"""
import os
from glob import glob
from subprocess import run
import numpy as np
from scipy import interpolate
import matplotlib as mpl
import matplotlib.pyplot as plt
from PIL import Image


def read(fid):
    """
    Simple read and parse of the XYZ data files using numpy
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
    cbar = plt.colorbar(tcf, label=kwargs["cbar_title"], shrink=0.5,aspect=8, 
                        pad=.03, ticks=[0, max_val/2, max_val])
    cbar.ax.yaxis.set_offset_position("left")
    # cbar.ax.set_yticklabels(["0", "2", ">4"])
    cbar.update_ticks()

    # Accoutrement
    plt.xlabel("Longitude")
    plt.ylabel("Latitude", rotation=90)
    for axis in ["top", "bottom", "left", "right"]:
        ax.spines[axis].set_linewidth(2)


def srcrcv(source=None, receiver=None):
    """
    Plot the source and receiver as simple markers
    """
    if source:
        plt.scatter(source[0], source[1], marker="*", s=70, color="k", 
                    linewidth=1, edgecolors="k")
    if receiver:
        plt.scatter(receiver[0], receiver[1], marker="v", s=50, color="w", 
                    linewidth=1.5, edgecolors="k")


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


def nznorth_extras(f, ax, text=None, faults=None):
    """
    Plot extra features on NZNorth such as labels, coastline, bathymetry,
    whatever
    """
    # Plot coastline for NZ ATOM North
    coast_fid = ("/Users/Chow/Documents/academic/vuw/data/carto/coastline/"
                 "extras/coast_latlon.txt")
    lat, lon = np.loadtxt(coast_fid).T
    plt.scatter(lon, lat, c="k", s=.005, marker=".")

    if text:
        plt.text(x=178.4, y=-42.4, s=text, color="k", fontsize=10, 
                 verticalalignment="bottom", horizontalalignment="right")

    # Plot active faults
    # if not faults:
    #     faults_fid = ("/Users/Chow/Documents/academic/vuw/data/carto/"
    #               "fault_coordinates/forest_mesh_gns_active_faults.txt")
    #     faults = np.loadtxt(faults_fid).T

    # f, lat, lon = faults
    # for f_ in np.unique(f):
    #     idx = np.where(f == f_)
    #     plt.plot(lon[idx], lat[idx], c="k", linewidth=.5)


if __name__ == "__main__":
    # =========================================================================
    # ACTIONS
    make_pngs = 1
    make_gif = 0
    # =========================================================================
    # PARAMETER SET HERE
    choice = "2018p130600"
    input_path = f"./inputs/{choice}/norm_vel"
    output_path = "./simulation"
    gif_fid = "sim_mov.gif"
    file_ext = ".xyz"
    min_val = 9e-7
    dt = .0125
    normalized = True
    gif_duration_ms = 200  # milliseconds
    trial = 0
    if choice == "2018p130600":
        max_val = 4E-5
        source = [176.300, -39.949]  # 2018p130600
        receiver = [177.674, -39.022]  # NZ.KNZ
        text = ("2018p130600\n"
                "2018-02-18T07:43:48Z\n"
                "M5.2; 21km depth")
    elif choice == "2016p105478":
        max_val = 2.E-5
        source = [173.075, -42.068]  # 2016p105478
        receiver = [178.257, -38.072]  # NZ.PUZ
        text = ("2016p105478\n"
                "2016-02-09T00:39:00Z\n"
                "M5.7, 48km depth")
    elif choice == "2016p881118":
        max_val = 8e-4
        source = [177.229, -40.674]  # 2016p881118
        receiver = [177.528, -38.334]  # NZ.MWZ
        text = ("2016p881118\n"
                "2016-11-22T00:19:42Z\n"
                "M5.5, 29km depth")

    # =========================================================================
    # Controls on colorbar
    if normalized:
        kwargs = {"cmap": "gist_heat_r",
                  "norm": plt.Normalize(0, max_val), 
                  "levels": 101,
                  "cbar_title": "norm of velocity [m/s]"
                  }
    else:
        # For single component 
        kwargs = {"cmap": "seismic", 
                  "norm": plt.Normalize(-1 * max_val, max_val), 
                  "levels": 100,
                  "cbar_title": "z velocity [m/s]"
                  }

    files = []
    # Test files to sample random data points to get an idea of relative amps
    if trial:
        files = ["gmt_movie_000080.xyz",
                 "gmt_movie_005600.xyz",
                 "gmt_movie_010000.xyz",
                 "gmt_movie_019680.xyz",]
        files = [os.path.join(input_path, _) for _ in files]
    # =========================================================================

    # Prep the file system
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    if not files:
        files = find(input_path, file_ext)
        
    # Make the .png files
    faults = None
    if make_pngs:
        ts_max = time_step(files[-1], dt=dt)
        for file_ in files:
            # Set up the plot
            ts = time_step(file_, dt=dt)
            x, y, z = read(file_)


            f, ax = plt.subplots(1)
            ax.set_aspect(1)
            print(f"{int(ts):0>3}/{int(ts_max)}")

            # Plot that ish
            plot(ax, x, y, z, min_val=min_val, max_val=max_val,  **kwargs)
            srcrcv(source, receiver)
            faults = nznorth_extras(f, ax, text, faults)

            plt.title(f"t={ts:6.2f} s")
            plt.xlim([173, 178.5])
            plt.ylim([-42.5, -37])

            # Clean up the end
            fid_out = os.path.join(output_path, 
                                   os.path.basename(file_) + ".png")
            f.tight_layout()
            plt.savefig(fid_out)
            plt.close()

    # Create the output .gif
    if make_gif:
        gif(output_path, gif_duration_ms)

