"""
Plot the outputs of Seisflows output.stats
To do:
    Slope plot as a rose diagram
"""
import os
import glob
import numpy as np
import matplotlib.pyplot as plt

for fid in glob.glob("*"):
    tag = os.path.basename(fid)
    if ".py" in tag:
        continue
    vals = np.loadtxt(fid)
    plt.plot(vals, "ko-")
    plt.xlabel("Iteration")
    plt.ylabel(tag)
    plt.title(tag)
    f = plt.gcf()
    f.tight_layout()
    plt.show()
