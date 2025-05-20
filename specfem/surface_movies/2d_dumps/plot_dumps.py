"""
Plot wavefield dumps, this is designed for the Venus work but should be
generalizable to all SPECFEM2D dumps
"""
import os
from glob import glob
import numpy as np
import matplotlib.pyplot as plt


# PARAMETERS
output = "./figures"
cmap = "Spectral"
levels = 125
plot_choice = "x"  # x, z, m[agnitude]
show_colorbar = True
dt = 0.25 # used for labeling time in title
vmax = 2E-5
vthresh = 1E-7

if not os.path.exists(output):
    os.makedirs(output)

x, z = np.loadtxt("wavefield_grid_for_dumps.txt").T
x *= 1E-3  # m -> km
z *= 1E-3  # m -> km

# Now go and plot the data
# for fid in sorted(glob("wavefield*.txt")):
for fid in ["wavefield0008000_01.txt"]:
    if "dumps" in fid:
        continue
    x_data, z_data = np.loadtxt(fid).T

    if plot_choice == "x":
        c = x_data
    elif plot_choice == "z":
        c = z_data
    elif plot_choice == "m":
        c = x_data + z_data

    # Don't plot amplitudes below a certain threshold
    c[np.abs(c) < vthresh] = 0

    if vmax is None:
        vmax = max([abs(c.min()), abs(c.max())])

    f = plt.figure(figsize=(5, 5))
    p = plt.tricontourf(x, z, c, levels=levels, cmap=cmap, vmin=-1*vmax, 
                        vmax=vmax)
    if show_colorbar:
        cbar = plt.colorbar(p, shrink=0.8, 
                            label=f"{plot_choice.upper()} Disp. (m)",
                            format="%.0E")
    ax = plt.gca()
    for name in ["top", "right"]:
        ax.spines[name].set_visible(False)
    for name in ["bottom", "left"]:
        ax.spines[name].set_linewidth(1)
    ax.set_aspect("equal")

    # Figure out what time it is
    t = int(os.path.basename(fid).split("_")[0][10:])
    plt.title(f"{t*dt:.2f} s")
    plt.xlabel("X [km]")
    plt.ylabel("Z [km]")
    
    plt.tight_layout()
    plt.savefig(os.path.join(output, 
                             os.path.basename(fid).replace("txt", "png")),
                dpi=200)
    plt.close("all")
    # plt.show()
