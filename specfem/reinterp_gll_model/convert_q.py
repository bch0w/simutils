"""
SPECFEM interpolation output Q values in mu and kappa, need to convert to Qp 
and Qs using equation from Dahlen and Tropm Eq. 9.59
"""
import os
import numpy as np
from glob import glob


for fid in glob("mu_*.xyz"):
    data = np.loadtxt(fid)
    vp = data[:, 3]
    vs = data[:, 4]
    qmu = data[:, 6]
    qkp = data[:, 7]

    # Qs == qmu
    qs = qkp

    # Qp from relationship with other quants
    f = 4 / 3 * (vs / vp) ** 2 
    qp = 1 / ((1 - f) / qkp + f / qmu)

    print(f"qp: {qp.min()} -> {qp.max()}")
    print(f"qmu: {qmu.min()} -> {qmu.max()}")

    data[:, 6] = qp
  
    fidout = fid.replace("mu_kappa_", "")
    np.savetxt(fidout, data, 
               "%10.3f %11.3f %10.3f %8.3f %8.3f %8.3f %8.3f %8.3f")



