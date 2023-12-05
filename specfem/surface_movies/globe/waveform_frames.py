"""
Relatively simple plotting script to plot waveforms frame by frame to make a 
.gif of the waveform being recorded progressively. Looks nice in combination 
with a simulation movie that shows the source and receiver locations.
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from pyasdf import ASDFDataSet
from pyatoa import Manager
from subprocess import run


# Parameters
delta = .0325
frames = 616
t_offset = 20  # to ensure that the origin time is 0 and not a negative number
component = "Z"
input_file = "/Users/chow/Work/work/akatom/nakversion/datasets/AK_E24K.h5"
# rcv_station = "TA.C18K"
rcv_station = "TA.A21K"

with ASDFDataSet(input_file) as ds:
    mgmt = Manager(ds=ds)
    mgmt.load(rcv_station, "i01/s00")

mgmt.standardize().preprocess(normalize_to="syn")
obs = mgmt.st_obs.select(component=component)[0].data
obs /= obs.max()
syn = mgmt.st_syn.select(component=component)[0].data
syn /= syn.max()
times = mgmt.st_obs[0].times() - t_offset

# Loop through times and plot consecutively
for idx in np.linspace(0, len(times), frames + 1):
    idx = int(idx)  
    f, ax = plt.subplots(1, dpi=200)

    # Plot the progressive waveform and a tracking marker
    plt.plot(times[:idx], syn[:idx], c="r", linewidth=2, zorder=9, 
             label="synthetic")
    plt.scatter(times[idx], syn[idx], marker="o", linewidth=1, c="r",
                zorder=10, edgecolor="k")

    # Plot the entire observed seismogram in the background
    plt.plot(times, obs, c="k", linewidth=1, zorder=8, label="observed")

    # Make sure all the plots look the same
    plt.xlim([times.min(), times.max()])
    # Put a bit of a buffer on the Y axes
    plt.ylim([-1.15, 1.15])
    plt.xlabel("time [s]")
    plt.ylabel("normalized Green's function")
    plt.legend()
    plt.title(f"SOURCE: TA_E24K; RCV: TA_C18L; COMP: {component}")

    ax.tick_params(which="major", length=5, width=1.5, direction="in",
                   bottom=True, top=True, left=True, right=True)
    ax.tick_params(which="minor", length=2, width=1.5, direction="in",
                   bottom=True, top=True, left=True, right=True)

    for axis in ["top", "bottom", "left", "right"]:
        ax.spines[axis].set_linewidth(2)

    f.tight_layout()
    plt.savefig(f"./waveforms/wav_{int(idx):0>5}.png")
    plt.close()


# Make the gif
if False:
    os.chdir("./waveforms")
    call = f"convert -delay 0 -loop 1 *.png wavmov.gif"
    run(call.split(" "))

