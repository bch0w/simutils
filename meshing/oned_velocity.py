"""
Derive a 1D velocity model from a 3D velocity model formatted in XYZ
"""
import numpy as np


with open("output_1d.txt", "w") as f:
    for level in ["shallow", "crust", "mantle"]:
        fid = f"tomography_model_{level}.xyz"
        header = np.load(fid + ".npz")
        data = np.load(fid + ".npy")
        zvals = sorted(np.unique(data[:, 2]))[::-1]
        for zval in zvals:
            zdata = np.where(data[:, 2] == zval)[0]
            vp_min = data[zdata, 3].min()
            vp_max = data[zdata, 3].max()
            vs_min = data[zdata, 4].min()
            vs_max = data[zdata, 4].max()

            print(f"{zval}\t{vp_min}\t{vp_max}\t{vs_min}\t{vs_max}\n")
            f.write(f"{zval}\t{vp_min}\t{vp_max}\t{vs_min}\t{vs_max}\n")
        
    
