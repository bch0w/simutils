"""
SPECFEM3D_GLOBE projects the domain to a unit sphere. To make slices in 
ParaView we therefore need to determine the unit sphere radius that corresponds
with the depth of interest
"""
import sys


val = float(sys.argv[1])
try:
    # 'r' for reverse
    reverse = sys.argv[2]
except IndexError:
    reverse = False

r_earth = 6371.  # km
if reverse:
    r_sphere = val
    depth_km = r_earth + (r_sphere * (1 - r_earth))
    print(f"{depth_km:.2f}")
else:
    depth_km = val
    r_sphere = (r_earth - depth_km) / r_earth
    print(f"{r_sphere:.4f}")


