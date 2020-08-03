"""
Values from Ristau 2008 SRL Table 1 North Island, to be added to the Meshfem3D 
Mesh_Par_file
"""
import numpy as np 
import sys
sys.path.append("../tools/")
from qps2qkm import qps2qkm 

# Values copied directly from Table 1 of Ristau (2008)
thickness = np.array([1, 2, 2, 10, 18, 31, 31.4, None])
vp_km_per_s = np.array([3., 4.4, 5.3, 6., 7.4, 7.78, 7.94, 8.08])
vs_km_per_s = np.array([1.7, 2.54, 3., 3.5, 4.3, 4.39, 4.51, 4.52])
density = np.array([2.29, 2.57, 2.69, 2.72, 2.87, 2.91, 2.92, 3.04])

# "Qp and Qs generally do not have a significant effect on MT inversion modeling,
# kept constant", Meshfem3D requires Qk and Qp, need to convert
qp = np.array([400] * len(thickness))
qs = np.array([200] * len(thickness))
qk, qm = qps2qkm(qp, qs, vp_km_per_s, vs_km_per_s, method="qps2qkm")

# Convert to the correct units for Meshfem
for i, vp in enumerate(vp_km_per_s):
    vp *= 1E3
    vs = vs_km_per_s[i] * 1E3
    rho = density[i] * 1E3 
   
    print(f"{len(vp_km_per_s)-i}  {rho:.0f}  {vp:.0f}  {vs:.0f}  {qk[i]:.0f}  {qm[i]:.0f}  0  2")
