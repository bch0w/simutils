"""
SPECFEM3D can output movie files as .xyz files which list LAT LON VAL
We can plot these as frames of a movie and collect them later as a .gif
Can also add text, coastlines etc. easily with Matplotlib
"""
import os
from glob import glob
import numpy as np
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


def plot(ax, x, y, z, show=True, save=False, **kwargs):
    """
    Plot the xyz file in the same fashion each time
    """
    tcf = ax.tricontourf(x, y, z, **kwargs) 
    cbar = plt.colorbar(tcf, format="%.2E", label="Disp. [m]")

    plt.xlabel("Longitude")
    plt.ylabel("Latitude", rotation=90)


def gif(path, duration, fid_out="output.gif"):
    """
    Generate a .gif file from all the resulting .png files
    """
    files = sorted(glob(os.path.join(path, "*.png")))
    img, *imgs = [Image.open(f) for f in files] 
    img.save(fp=fid_out, format="GIF", append_images=imgs, save_all=True,
             duration=duration, loop=1)


def nznorth_extras(f, ax):
    """
    Plot extra features on NZNorth such as labels, coastline, bathymetry,
    whatever
    """
    coast_fid = ("/Users/Chow/Documents/academic/vuw/data/carto/coastline/"
                 "extras/coast_latlon.txt")


if __name__ == "__main__":
    # =========================================================================
    # PARAMETER SET HERE
    input_path = "./z_disp"
    output_path = "./output"
    gif_fid = "nznorth_2018p13600_zdisp.gif"
    file_ext = ".xyz"
    kwargs = {"cmap": "seismic", 
              "norm": plt.Normalize(-8e-4, 8e-4), 
              "levels": 100
              }
    gif_duration_ms = 300  # milliseconds
    dt = .0125
    files = []
    files = ["./z_disp/gmt_movie_012600.xyz"]
    # =========================================================================
    # ACTIONS
    make_pngs = True
    make_gif = 0
    # =========================================================================

    # Prep the file system
    if not os.path.exists(output_path):
        os.mkdirs(output_path)
    if not files:
        files = find(input_path, file_ext)

    # Make the .png files
    if make_pngs:
        for file_ in files:
            # Set up the plot
            x, y, z = read(file_)
            f, ax = plt.subplots(1)
            ts = time_step(file_, dt=dt)

            # Plot that ish
            plot(ax, x, y, z, **kwargs)
            nznorth_extras(f, ax)
            plt.title(f"t={ts:6.2f} s")

            # Clean up the end
            fid_out = os.path.join(output_path, 
                                   os.path.basename(file_) + ".png")
            plt.savefig(fid_out)
            plt.close()

    # Create the output .gif
    if make_gif:
        gif(output_path, gif_duration_ms)

