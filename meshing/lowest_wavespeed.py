"""
Determine the lowest velocity for every depth slice of the initial model. This
allows getting at the minimum grid spacing required for the numerical mesh.
Need to be careful as the lowest velocity may be lowered through the inversion.

Data is: x, y, z, vp[m/s], vs[m/s], rho[kg/m**3], Qp, Qs
"""
import os
import json
import numpy as np
import matplotlib.pyplot as plt
from checkerboardiphy import xyz_reader


depths = {}
tags = ["mantle", "crust", "shallow"]
for t in tags:
    fid = f"tomography_model_{t}.xyz"
    header, data = xyz_reader(fid, save=True)
    unique_depths = np.unique(data[:, 2])
    min_velocities, max_velocities = [], []
    for d in unique_depths:
        print(f"{d} km")
        idx = np.where(data[:, 2] == d)[0]
        min_vp = data[idx, 3].min() 
        max_vp = data[idx, 3].max() 
        min_vs = data[idx, 4].min()
        max_vs = data[idx, 4].max()
        min_velocities.append(min(min_vs, min_vp))
        max_velocities.append(max(max_vs, max_vp))
    
    plt.plot(min_velocities, unique_depths, 'ko-')
    plt.xlabel("Min. Velocity [km/s]")
    plt.ylabel("Depths [km]")
    plt.savefig(f"{t}.png")
    plt.close()

        
