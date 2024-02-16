"""
SPECFEM3D_GLOBE projects the domain to a unit sphere. To make slices in 
ParaView we therefore need to determine the unit sphere radius that corresponds
with the depth of interest
"""
import sys
depth_km = float(sys.argv[1])
r_earth = 6371.  # km
r_sphere = (r_earth - depth_km) / r_earth
print(f"{r_sphere:.4f}")


