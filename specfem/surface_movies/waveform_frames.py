"""
Relatively simple plotting script to plot waveforms frame by frame to make a 
.gif of the waveform being recorded progressively. Looks nice in combination 
with a simulation movie that shows the source and receiver locations.
"""
import os
import sys
from glob import glob
import numpy as np
import matplotlib.pyplot as plt


# Parameters
dt = 1  # seconds
t_offset = 20  # to ensure that the origin time is 0 and not a negative number
trial = False
component = "Z"
input_files = glob(f"./inputs/input_wav/??.???.HX{component}.semd")
assert(len(input_files) == 1)
input_file = input_files[0]
print(input_file)

# Read in the 2 column ascii and set offset
time, data = np.loadtxt(input_file).T
time += t_offset

# Loop through times and plot consecutive
idx = []
for t in np.arange(time.min(), time.max(), dt):
    print(t)
    # Assuming there is only one corresponding time step
    try:
        tstep = np.where(time==t)[0][0]
    except IndexError:
        # Exact value is not found (usually for 0) find nearest value
        import ipdb;ipbd.set_trace() 
        pass
    idx.append(tstep)

    # Plot the progressive waveform and a tracking marker
    f, ax = plt.subplots(1, figsize=[6, 3.])
    plt.plot(time[idx], data[idx], c="r", linewidth=2, zorder=1)
    plt.scatter(time[tstep], data[tstep], marker="o", linewidth=1, c="k",
                zorder=10)

    # Make sure all the plots look the same
    plt.xlim([time.min(), time.max()])
    data_max = max(abs(data.min()), abs(data.max()))
    plt.ylim([-1 * data_max, data_max])
    plt.xlabel("time [s]")
    plt.ylabel("displacement [m]")
    plt.title(os.path.basename(input_file))


    ax.tick_params(which="major", length=5, width=1.5, direction="in",
                   bottom=True, top=True, left=True, right=True)
    ax.tick_params(which="minor", length=2, width=1.5, direction="in",
                   bottom=True, top=True, left=True, right=True)

    for axis in ["top", "bottom", "left", "right"]:
        ax.spines[axis].set_linewidth(2)

    f.tight_layout()
    plt.savefig(f"./waveforms/wav_{int(t):0>3}.png")
    plt.close()

    if trial:
        sys.exit()


