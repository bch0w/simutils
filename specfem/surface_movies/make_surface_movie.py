"""
Make movie using SPECFEM moviedata files, should be run from directly inside
the SPECFEM working directory
"""
import os
import cartopy.crs as ccrs
import cartopy.feature as cf
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from glob import glob


def run_xcreate_movie():
    """
    Run
    """
    call = "./bin/screate_movie_GMT_global <<< $'{options}'"
    assert(os.path.exists(call))
    options = [
        "1",    # first frame of movie
        "-1",   # last frame of movie
        "1",    # 1=Z, 2=N, 3=E
        "F",  # ASCII (F) or binary (T) 
        "T",  # mute source area (T) or not (F)
        "T",  # use moving average for normalization (T) or not (F)
        "1",    # absolute value for normalization
    ]
    options = "\n".join(options)
    os.system(call.format(options=options))

def read_data(fid):
    x, y, z = np.loadtxt(fid).T
    return x, y, z

def setup_map(x, y, dpi=100, axis_lw=0.75):
    """
    Setup the basemap to plot data onto
    """
    f = plt.figure(dpi=dpi)

    # [min_lon, max_lon, min_lat, max_lat]
    extent = [x.min(), x.max(), y.min(), y.max()]
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

    return f, ax

def tplot_contour(x, y, z, ax=None, min_max=(None, None), mask=(0, 0), 
                  cbar_title="", **kwargs):
    """
    Plot the xyz file in the same fashion each time. Use tricontourf because
    the data is not in a uniform grid so we use a triangular interpolation
    """
    # Mask out values before plotting
    tri = mpl.tri.Triangulation(x, y)

    import pdb;pdb.set_trace()

    # Determine where the amplitude array is within the mask range
    masked_vals = (z >= mask[0]) & (z <= mask[1])

    # Mask out the x, y values corresponding to the mask range
    mask = np.all(np.where(masked_vals[tri.triangles], True, False), axis=1)
    tri.set_mask(mask)

    # Set custom levels to keep the colorbar segmented the same across figs
    min_val, max_val = min_max
    if min_val:
        kwargs["levels"] = np.linspace(min_val, max_val, kwargs["levels"])
        try:
            tcf = ax.tricontourf(tri, z, vmin=min_val, vmax=max_val, 
                                 extend="both", **kwargs)
        # Edge case when we have masked out all values (e.g., prior origin time)
        except ValueError:
            tcf = ax.tricontourf(x, y, np.zeros(len(z)), vmin=min_val, 
                                 vmax=max_val, extend="both", **kwargs)
    else:
        tcf = ax.tricontourf(tri, z, extend="both", **kwargs)

    # Colorbar
    cbar = plt.colorbar(tcf, label=cbar_title, shrink=0.5, aspect=8, 
                        pad=.02, ticks=[min_val, 0, max_val])
    cbar.ax.yaxis.set_offset_position("left")

    # Mark off where the mask is thresholding on the colorbar
    # for val in mask:
    #     cbar.ax.plot([0, 1], [val, val], color="k", lw=1)

def plot_contour(x, y, z, ax):
    plt.contourf(x, y, z, transform=ccrs.PlateCarree())


def plot_frame(fid):
    """
    Plot a single frame given a data file
    """
    x, y, z = read_data(fid)
    f, ax = setup_map(x, y)
    plot_contour(x, y, z, ax)
    plt.show()
    # plt.savefig(f"{fid}.png")
    # plt.close(fig=f)


if __name__ == "__main__":
    plot_frame("gmt_movie_005500.Z.xyz")
