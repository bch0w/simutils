"""
Plot Movie frames outputted by SPECFEM3D_GLOBE with ASCII format
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from glob import glob


# Set parameters here
dt = 3.25E-2
path = "./"
path_out = "./figures"

if not os.path.exists(path_out):
    os.mkdir(path_out)

# Load spatial coordinates
x, y = np.loadtxt(os.path.join(path, "ascii_movie.xy")).T

# Get the min and max amplitude value for scaling
minval, maxval = np.inf, 0
for fid in sorted(glob(os.path.join(path, "ascii_movie_*.d"))):
    vals = np.loadtxt(fid)
    if vals.min() < minval:
        minval = vals.min()
    if vals.max() > maxval:
        maxval = vals.max()


# Pick every 100th frame to plot things quicker
files = sorted(glob(os.path.join(path, "ascii_movie_*.d")))
for i in range(0, len(files), 10):
    fid = files[i]
    # Get the timestep from the filename
    fid_strip = os.path.splitext(os.path.basename(fid))[0]
    nframe = int(fid_strip.split("_")[-1])
    timestep = nframe * dt

    # Get the movie amplitudes from data
    d = np.loadtxt(fid)

    # Plot the movie frame
    plt.scatter(x, y, c=d, cmap="jet", vmin=minval * 2, vmax=maxval / 3)
    plt.title(f"{timestep:.2f}s")
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")

    plt.savefig(os.path.join(path_out, f"movie_{nframe:0>5}.png"))

a=1/0

# Load each timestep one at a time
for fid in sorted(glob(os.path.join(path, "ascii_movie_*.d"))):
    # Get the timestep from the filename
    fid_strip = os.path.splitext(os.path.basename(fid))[0]
    nframe = int(fid_strip.split("_")[-1])
    timestep = nframe * dt

    # Get the movie amplitudes from data
    d = np.loadtxt(fid)

    # Plot the movie frame
    plt.scatter(x, y, c=d, cmap="jet", vmin=minval, vmax=maxval)
    plt.title(f"{timestep:.2f}s")
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")

    plt.savefig(os.path.join(path_out, f"movie_{nframe:0>5}.png"))



