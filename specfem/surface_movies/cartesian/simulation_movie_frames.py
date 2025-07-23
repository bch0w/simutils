"""
SPECFEM3D can output movie files as .xyz files which list LAT LON VAL
We can plot these as frames of a movie and collect them later as a .gif
Can also add text, coastlines etc. easily with Matplotlib
"""
import argparse
import os
import glob
import numpy as np
import matplotlib.pyplot as plt

from concurrent.futures import ProcessPoolExecutor, wait
from PIL import Image
from pyproj import Proj


def parse_args():
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser()

    # PATHS
    parser.add_argument("-D", "--data", default="DATA", type=str,
                        nargs="?", help="Path to SPECFEM DATA/ directory")
    parser.add_argument("-m", "--mov_files", type=str, nargs="?",
                        help="Path to moviedata*.xyz text files")
    parser.add_argument("-o", "--output", default="./frames", type=str,
                        nargs="?", help="Path to store output files")

    # PARAMETERS
    parser.add_argument("--min_val", default=None, type=float, nargs="?",
                        help="force minimum value on cmap")
    parser.add_argument("--max_val", default=None, type=float, nargs="?",
                        help="force maximum value on cmap")
    parser.add_argument("-c", "--cmap", default=None, type=str, nargs="?",
                        help="Colormap to plot image with")
    parser.add_argument("-u", "--utm", default=None, type=int, nargs="?",
                        help="UTM zone to convert coordinates if given")

    # FLAGS
    parser.add_argument("-f", "--frames", nargs="+", type=int,
                        help="start and end index of frames you want to make")
    parser.add_argument("-p", "--parallel", action="store_true",
                        help="run plotting in parallel with subprocess")
    parser.add_argument("-d", "--dry_run", type=int, default=None, 
                        help="Test out with a user-selected number of frames")
    parser.add_argument("--plot_source", type=bool, default=False,
                        help="Plot the source location")
    parser.add_argument("--plot_stations", type=bool, default=False,
                        help="Plot the station locations")

    return parser.parse_args()


def parse_data(path="./DATA"):
    """
    Get some key information from the DATA directory
    """
    # Par_file
    with open(os.path.join(path, "Par_file"), "r") as f:
        for line in f.readlines():
            if line.startswith("NSTEP"):
                nstep = float(line.strip().split("=")[-1])
            if line.startswith("DT"):
                dt = float(line.strip().split("=")[-1])
                break

    # CMTSOLUTION
    with open(os.path.join(path, "CMTSOLUTION"), "r") as f:
        for line in f.readlines():
            if line.startswith("latitude"):
                src_lat = float(line.strip().split(":")[-1])
            if line.startswith("longitude"):
                src_lon = float(line.strip().split(":")[-1])
                break
    source = {"lat": src_lat, "lon": src_lon}

    # STATIONS
    stations = {}
    with open(os.path.join(path, "STATIONS"), "r") as f:
        for line in f.readlines():
            sta, net, lat, lon, *_ = line.strip().split()
            stations[f"{net}.{sta}"] = {"lat": float(lat), "lon": float(lon)}

    return nstep, dt, source, stations


def convert_coords(lon, lat, utm_zone=None):
    """
    Convert from the native Lat/Lon coordinate system into UTM 60S coordinates
    This only needs to be done for the first file because the remaining files
    are assumed to follow the same coordinate points.
    """
    if utm_zone < 0:
        south=True
    elif utm_zone > 0:
        south = False

    projection = Proj(proj="utm", zone=abs(utm_zone), south=south, 
                      ellps="WGS84", preserve_units=False)
    x, y = projection(lon, lat, inverse=False)

    # Zero out the origin
    x -= x.min()
    y -= y.min()

    # Convert units of m to km
    x /= 1E3
    y /= 1E3

    return x, y


def plot_frame(fid, x, y, dt, min_val=None, max_val=None, source=None, 
               stations=None, utm=False, save="./", dpi=200, norm=None, 
               cmap="viridis", levels=100):
    """
    Plot the xyz file in the same fashion each time. Use tricontourf because
    the data is not in a uniform grid so we use a triangular interpolation
    """
    fid_out = os.path.join(save, os.path.basename(fid) + ".png")
    if os.path.exists(fid_out):
        return


    # Get some values from the actual data file
    z = np.loadtxt(fid, usecols=-1) 
    f, ax = plt.subplots(1)

    # Mask out values before plotting
    tcf = ax.tricontourf(x, y, z, norm=norm, cmap=cmap, levels=levels,
                         vmin=min_val, vmax=max_val) 

    # !!! FIX THIS, NEED TO FIND OUT HOW TO GET MIN ABS VALUE FOR MASKING
    # if mask_pct:
    #     raise NotImplementedError("fix me bro")
    #     tri = mpl.tri.Triangulation(x, y)
    #     mask_val = z.max() * mask_pct
    #     masked_vals = np.where((z > mask_val) & ())
    #     mask = np.all(np.where(masked_vals[tri.triangles], True, False), axis=1)
    #     tri.set_mask(mask)
    #     # Set custom levels to keep the colorbar segmented the same across figs
    #     kwargs["levels"] = np.linspace(0, max_val, kwargs["levels"])
    #     try:
    #         tcf = ax.tricontourf(tri, z, vmin=0, vmax=max_val, extend="max", 
    #                              **kwargs)
    #     except ValueError:
    #         # ValueError occurs when we mask out all the values (e.g. T < 0s)
    #         tcf = ax.tricontourf(x, y, np.zeros(len(z)), vmin=0, vmax=max_val,
    #                              extend="max", **kwargs) 
    # # Or just plot the values straight up
    # else:
    #     tcf = ax.tricontourf(x, y, z, **kwargs) 

    # Source and stations
    if source:        
        plt.scatter(source["lon"], source["lat"], marker="*", s=150, color="y", 
                    linewidth=1, edgecolors="k"
                    )
    if stations:
        for key, val in stations.items():
            stax, stay = val["lon"], val["lat"]
            plt.text(stax, stay, key, fontsize=4, alpha=0.5, ha="center")
            plt.scatter(stax, stay, marker="v", s=20, color="w", linewidth=0.5, 
                        edgecolors="k", alpha=0.5)
            
    # Colorbar
    if min_val > 0:
        ticks = np.linspace(min_val, max_val, 5)
    else:
        ticks = [min_val, 0, max_val]
    cbar = plt.colorbar(tcf, shrink=0.5, aspect=8, 
                        pad=.03, ticks=ticks, format="%4.2E")
    cbar.ax.yaxis.set_offset_position("left")
    ticklabs = cbar.ax.get_yticklabels()
    cbar.ax.set_yticklabels(ticklabs, fontsize=6)

    # Mark where the min value threshold is set
    # if mask_pct:
    #     cbar.ax.plot([0, 1], [min_val, min_val], "cyan")

    cbar.update_ticks()

    # Accoutrement
    timestep = int(os.path.basename(fid).split("_")[2].split(".")[0]) * dt
    ax.set_title(f"t={timestep:.2f}s")
    if utm:
        plt.xlabel("X [m]")
        plt.ylabel("Y [m]")
    else:
        plt.xlabel("Longitude")
        plt.ylabel("Latitude")

    for axis in ["top", "bottom", "left", "right"]:
        ax.spines[axis].set_linewidth(2)
    ax.set_aspect(1)
    plt.xlim([x.min(), x.max()])
    plt.ylim([y.min(), y.max()]) 
    f.tight_layout()

    plt.savefig(fid_out, dpi=dpi)

    plt.close()


# def gif(path, duration, fid_out="output.gif"):
#     """
#     Generate a .gif file from all the resulting .png files
#     PIL runs into a '[Errno 24] Too many open files' for large numbers of files
#     So use ImageMagick convert wrapped with subprocess instead for those
#     """
#     files = sorted(glob(os.path.join(path, "*.png")))
#     if len(files) < 256:
#         img, *imgs = [Image.open(f) for f in files] 
#         img.save(fp=fid_out, format="GIF", append_images=imgs, save_all=True,
#                  duration=duration, loop=1)
#     else:
#         os.chdir(path)
#         call = f"convert -delay 0 -loop 1 *.png {fid_out}"
#         run(call.split(" "))


if __name__ == "__main__":
    args = parse_args()

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    print("parsing data from files")
    # Get information from SPECFEM
    nstep, dt, source, stations = parse_data(args.data)

    files = sorted(glob.glob(os.path.join(args.mov_files, "*")))

    # Determine min and max val out of all the movie frames
    print("getting min/max amplitudes from movie files")
    minval = np.inf
    maxval = 0
    for i in np.arange(0, len(files), 1, dtype=int):
        arr = np.loadtxt(files[i], usecols=-1)
        if arr.min() < minval:
            minval =  arr.min()
        if arr.max() > maxval:
            maxval = arr.max()
    print(f"minimum value = {minval}")
    print(f"maximum value = {maxval}")

    if args.dry_run:
        idxs = np.linspace(0, len(files)-1, args.dry_run, dtype=int)
        files = np.array(files)[idxs]
    elif args.frames:
        start, end = args.frames
        files = files[start:end]
    
    print("determining amplitude bounds for colorbars")
    # Parse through the movie files and find the min and max values
    # Check every 10 seconds to avoid reading every file
    if not args.dry_run:
        df = int(10 * 1 // dt)  # 10 seconds in samples
    else:
        df = 1

    # Figure out if we have a diverging colormap (like velocity) 
    if minval < 0:
        levels = 101  # 1 for zero
        cmap = args.cmap or "seismic"
        _cmap_max = np.max([np.abs(minval), np.abs(maxval)])
        norm = plt.Normalize(-1 * _cmap_max, _cmap_max)  # ensure 0 at center 
    # or sequential (like norm of velocity)
    else:
        levels = 100
        cmap = args.cmap or "viridis"
        norm = plt.Normalize(0, maxval)

    # Assuming all the files have the same grid so we only need to read geo
    # information one time
    lon, lat, _ = np.loadtxt(files[0]).T 

    # Convert to UTM coordinates if necessary
    if args.utm:
        print(f"converting coordinates to UTM {args.utm}")
        x, y = convert_coords(lon, lat, args.utm)

        srcx, srcy = convert_coords(source["lon"], source["lon"], args.utm)
        source = {"lat": srcy, "lon": srcx}

        stations_xy = {}
        for sta, val in stations.items():
            stax, stay = convert_coords(val["lon"], val["lat"], args.utm)
            stations_xy[sta] = {"lat": stay, "lon": stax}
        stations = stations_xy
    else:
        x, y = lon, lat

    # Read through and plot each of the moviefiles
    print("plotting movie frames")
    if args.parallel:
        with ProcessPoolExecutor(max_workers=os.cpu_count()-1) as executor:
            futures = [executor.submit(plot_frame, fid, x, y, dt, 
                                       minval, maxval, source, stations, 
                                       args.utm, args.output, 200,
                                       norm, cmap, levels) 
                                       for fid in files]
        wait(futures)
        futures[0].result()
    else:
        for i, fid in enumerate(files):
            plot_frame(fid=fid, x=x, y=y, dt=dt, min_val=minval, max_val=maxval, 
                       source=source, stations=stations, 
                       levels=levels, norm=norm, cmap=cmap, utm=args.utm,
                       save=args.output)
   
