"""
Plot BM base model file
"""
import sys
import numpy as np
import matplotlib.pyplot as plt


fid = sys.argv[1]
arr = np.loadtxt(fid, skiprows=5)
r, rho, vp, vs = arr.T


plt.plot(vs, r*1e-3, "ro-")
plt.show()

