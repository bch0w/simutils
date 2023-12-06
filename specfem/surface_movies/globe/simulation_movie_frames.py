"""
SPECFEM3D can output movie files as .xyz files which list LAT LON VAL
We can plot these as frames of a movie and collect them later as a .gif
Can also add text, coastlines etc. easily with Matplotlib
"""
import os
import cartopy.crs as ccrs
import cartopy.feature as cf
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor
from glob import glob
from subprocess import run
from PIL import Image


def read_xy(fid):
    """
    Simple read of the spatial coordinates
    """
    data = np.loadtxt(fid)
    x, y = data.T
    return x, y


def read_data(fid):
    """
    Simple read and parse of the XYZ data files using numpy
    :return: lon, lat, z
    """
    data = np.loadtxt(fid)
    return data


def find(path="./", ext=".d"):
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


def plot_contour(ax, x, y, z, min_val=None, max_val=None, mask=(0, 0), 
                 cbar_title="", **kwargs):
    """
    Plot the xyz file in the same fashion each time. Use tricontourf because
    the data is not in a uniform grid so we use a triangular interpolation
    """
    # Mask out values before plotting
    tri = mpl.tri.Triangulation(x, y)

    # Determine where the amplitude array is within the mask range
    masked_vals = (z >= mask[0]) & (z <= mask[1])

    # Mask out the x, y values corresponding to the mask range
    mask = np.all(np.where(masked_vals[tri.triangles], True, False), axis=1)
    tri.set_mask(mask)
    # Set custom levels to keep the colorbar segmented the same across figs
    kwargs["levels"] = np.linspace(min_val, max_val, kwargs["levels"])
    try:
        tcf = ax.tricontourf(tri, z, vmin=min_val, vmax=max_val, 
                             extend="both", **kwargs)
    # Edge case when we have masked out all values (e.g., prior to origin time)
    except ValueError:
        tcf = ax.tricontourf(x, y, np.zeros(len(z)), vmin=min_val, vmax=max_val,
                             extend="both", **kwargs)

    # Colorbar
    cbar = plt.colorbar(tcf, label=cbar_title, shrink=0.5, aspect=8, 
                        pad=.02, ticks=[min_val, 0, max_val])
    cbar.ax.yaxis.set_offset_position("left")

    # Mark off where the mask is thresholding on the colorbar
    # for val in mask:
    #     cbar.ax.plot([0, 1], [val, val], color="k", lw=1)


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


# !!! Set parameters here, pretty hacky but it avoids passing a bunch of args
# !!! through concurrent.futurdef plot(fid, x, y, input_path="./data", output_path="./output", 
def plot(fid, x, y, input_path="./data", output_path="./simulation", 
         extent=(-168., -139., 63.5, 72),  source=(-148.4868, 68.0748), dpi=300, 
         axis_lw=2, min_val=-1e-5, max_val=1e-5, mask=(-1e-7, 1e-7), 
         dt=3.25E-2,  **kwargs):
    """
    Main plotting function for a single frame. This can be parallelized as it 
    doesn't require information from other frames

    EXTENT: 
    FORCE001: -147.8616, 64.8736
    """
    # Set up the plot
    plt.figure(dpi=dpi)
    ts = time_step(fid, dt=dt)
    z = read_data(fid)

    # Initialize Cartopy map with coastline and borders
    min_lon, max_lon, min_lat, max_lat = extent
    central_longitude = (min_lon + max_lon) / 2
    central_latitude = (min_lat + max_lat) / 2
    ref_proj = ccrs.PlateCarree()
    projection = ccrs.Stereographic(central_longitude=central_longitude,
                                    central_latitude=central_latitude)
    
    ax = plt.axes(projection=projection)
    ax.set_extent(extent, crs=ref_proj)
    ax.coastlines(lw=axis_lw)
    ax.add_feature(cf.BORDERS)

    # Plot the contour plot of the surface velocity field
    plot_contour(ax, x, y, z, min_val=min_val, max_val=max_val, 
                 mask=mask, transform=ref_proj, cmap="seismic", 
                 norm= plt.Normalize(min_val, max_val), levels=256,
                 cbar_title="Vertical Velocity [m/s]")

    # Plot source and receivers
    plt.scatter(source[0], source[1], transform=ref_proj, marker="*", 
                s=80, color="y", linewidth=0.75,  edgecolors="k",
                zorder=10)
    
    stax, stay = np.loadtxt(os.path.join(input_path, "STATIONS"), 
                            usecols=(3, 2)).T
    plt.scatter(stax, stay, transform=ref_proj, marker="v", 
                s=10, color="None", linewidth=0.75,  edgecolors="k",
                zorder=10)

    # Set some plot attributes
    plt.text(0.05, 0.9, f"t={ts:6.2f} s", transform=ax.transAxes,
                fontsize=12, color="k")
    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False,
                        y_inline=False, linewidth=.25, alpha=0.25, 
                        color="k")
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {"rotation": 0}
    ax.spines["geo"].set_linewidth(axis_lw)

    # Clean up the end
    fid_out = os.path.join(output_path, os.path.basename(fid) + ".png")
    print(fid_out)
    plt.savefig(fid_out)
    plt.close()


if __name__ == "__main__":
    # =========================================================================
    make_pngs = bool(input("make .pngs? [y/n]") == "y")
    make_gif = bool(input("make .gif? [y/n]") == "y")
    test_run = 0
    parallel = True

    input_path = "./data"
    output_path = "./simulation"
    file_ext = ".d"

    gif_fid = "sim_mov.gif"
    gif_duration_ms = 200  # milliseconds
    # =========================================================================
    # Prep the file system
    if make_pngs:
        files = []
        if not os.path.exists(output_path):
            os.mkdir(output_path)
        if not files:
            files = find(input_path, file_ext)

        if test_run:
            n = int(len(files)//2)
            files = files[n-5:n+5]
            
        # Get some geographic information
        x, y = read_xy(os.path.join(input_path, "ascii_movie.xy"))
        if parallel:
            with ProcessPoolExecutor() as executor:
                futures = [executor.submit(plot, file_, x, y) 
                           for file_ in files]
            for future in futures:
                future.result()
        else:
            for i, file_ in enumerate(files):
                plot(file_, x, y)

    # Create the output .gif
    if make_gif:
        gif(output_path, gif_duration_ms, gif_fid)

