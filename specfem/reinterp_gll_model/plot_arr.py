"""
Small script to read in and plot slices of a GLL model that has been projected
to a regular FD grid using the SPECFEM3D Cartesian auxiliary function 
xproject_and_combine_vol_data_on_regular_grid
"""
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
    
    
CHOICE = "vs"
TAG = "raw"
CMAP = "Spectral"


def _find_nearest(array, value):
    """find the nearest value in a given array"""
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


def plot_single(arr, choice, val, section):
    """"""
    choiceval = ["x", "y", "z"][choice]

    trueval = _find_nearest(arr[:, choice], val)
    print(f"{choiceval}: {val} -> {trueval}")

    # Raw data
    dslice = arr[np.where(arr[:, choice] == trueval)]
    xyd = np.delete(dslice, choice, axis=1)

    # Create figures
    f = plt.figure(1)

    # Left subplot is raw data
    x_ = np.unique(xyd[:, 0]) * 1E-3
    y_ = np.unique(xyd[:, 1]) * 1E-3
    z_ = np.reshape(xyd[:, 2], (len(y_), len(x_))) * 1E-3

    sc = plt.contourf(x_, y_, z_, levels=128, cmap=CMAP)

    plt.colorbar(sc, label="Vs [km/s]", pad=0.008)

    # Finish off the plot
    if choice == 0:
        plt.xlabel("UTM-60 Northing [km]")
        plt.ylabel("Depth [km]")
    elif choice == 1:
        plt.xlabel("UTM-60 Easting [km]")
        plt.ylabel("Depth [km]")
    elif choice == 2:
        plt.xlabel("UTM-60 Easting [km]")
        plt.ylabel("UTM-60 Northing [km]")

    plt.title(f"{section.title()} ({choiceval.upper()}~={val*1e-3:.2f}km)")

    plt.savefig(f"{section}_{CHOICE}_{choiceval}_{int(val)}.png")

    plt.close()


def plot2d(arr1, arr2, choice, val):
    """2D depth slices at a given Z value""" 
    choiceval = ["x", "y", "z"][choice]

    trueval = _find_nearest(arr1[:, choice], val)
    print(f"{choiceval}: {val} -> {trueval}")

    # Raw data
    dslice = arr1[np.where(arr1[:, choice] == trueval)]
    raw_xyd = np.delete(dslice, choice, axis=1)

    # Smooth data
    dslice = arr2[np.where(arr2[:, choice] == trueval)]
    smt_xyd = np.delete(dslice, choice, axis=1)

    # Find min and max values for colorscales
    min_val = min([raw_xyd[:,2].min(), smt_xyd[:,2].min()])
    max_val = max([raw_xyd[:,2].max(), smt_xyd[:,2].max()])

    # Create figures
    f = plt.figure(1, figsize=(12.8,4.8))

    # Left subplot is raw data
    s1 = plt.subplot(121)
    x_ = np.unique(raw_xyd[:, 0])
    y_ = np.unique(raw_xyd[:, 1])
    z_ = np.reshape(raw_xyd[:, 2], (len(y_), len(x_)))
    sc1 = plt.contourf(x_, y_, z_, levels=128, cmap=CMAP, 
                       vmin=min_val, vmax=max_val)
    # plt.colorbar(sc1, label="Vs [m/s]")
    plt.title("original model")

    # Right subplot is smooth
    s1 = plt.subplot(122)
    x_ = np.unique(smt_xyd[:, 0])
    y_ = np.unique(smt_xyd[:, 1])
    z_ = np.reshape(smt_xyd[:, 2], (len(y_), len(x_)))
    sc2 = plt.contourf(x_, y_, z_, levels=128, cmap=CMAP,
                       vmin=min_val, vmax=max_val)
    # plt.colorbar(sc2, label="Vs [m/s]")
    plt.title("smoothed model")
    plt.gca().axes.yaxis.set_ticklabels([])

    # Colorbar
    f.subplots_adjust(right=0.8)
    cbar_ax = f.add_axes([0.85, 0.15, 0.02, 0.7])
    f.colorbar(sc2, cax=cbar_ax, label="Vs [m/s]")

    # Finish off the plot
    plt.subplots_adjust(wspace=0.25)
    plt.suptitle(f"{CHOICE} ({choiceval}~={val*1e-3:.2f}km)")

    plt.savefig(f"compare_{CHOICE}_{choiceval}_{int(val)}.png")
    # plt.show()

    plt.close()


if __name__ == "__main__":
    section = sys.argv[1]  # "crust"
    choicelist = ["x", "y", "z", "vp", "vs", "rho", "qmu", "qkappa"]

    #if section == "shallow":
    if False:
        # Read in data for both raw and smoothed models
        raw_data = np.loadtxt(f"raw_{section}_tomo_file.xyz")
        smt_data = np.loadtxt(f"smooth_{section}_tomo_file.xyz")

        # Subset data array for a given parameter
        raw_arr = raw_data[:, [0, 1, 2, choicelist.index(CHOICE)]]
        smt_arr = smt_data[:, [0, 1, 2, choicelist.index(CHOICE)]]

        # Get rid of border areas which have zero value and tell User new bounds
        raw_arr = raw_arr[np.where(raw_arr[:,-1] != 0)]
        smt_arr = smt_arr[np.where(smt_arr[:,-1] != 0)]

        for xval in [raw_arr[:,0].min(), raw_arr[:,0].mean(), raw_arr[:,0].max()]:
            plot2d(raw_arr, smt_arr, choice=0, val=xval)

        for yval in [raw_arr[:,1].min(), raw_arr[:,1].mean(), raw_arr[:,1].max()]:
            plot2d(raw_arr, smt_arr, choice=1, val=yval)

        if section == "mantle":
            zvals = np.arange(-44E3, -400E3, -4E3)
        elif section == "crust":
            zvals = np.arange(-7E3, -50E3, -1E3)
        elif section == "shallow":
            zvals = np.arange(2.25E3, -8E3, -250)

        for zval in zvals:
            plot2d(raw_arr, smt_arr, choice=2, val=zval)
    else:
        # Read in data for both raw and smoothed models
        data = np.loadtxt(f"tomography_model_{section}.xyz")

        # Subset data array for a given parameter
        arr = data[:, [0, 1, 2, choicelist.index(CHOICE)]]

        # Get rid of border areas which have zero value and tell User new bounds
        arr = arr[np.where(arr[:,-1] != 0)]

        # for xval in [arr[:,0].min(), arr[:,0].mean(), arr[:,0].max()]:
        for xval in [arr[:,0].mean()]:
            plot_single(arr, choice=0, val=xval, section=section)

        # for yval in [arr[:,1].min(), arr[:,1].mean(), arr[:,1].max()]:
        for yval in [arr[:,1].mean()]:
            plot_single(arr, choice=1, val=yval, section=section)

        if section == "shallow":
            zvals = [-2.5E3, -5E3, -7.5E3] 
        elif section == "crust":
            zvals = [-15E3, -25E3, -35E3]  
        elif section == "mantle":
            zvals = [-45E3, -90E3, -150E3]  
        for zval in zvals:
            plot_single(arr, choice=2, val=zval, section=section)
                   

