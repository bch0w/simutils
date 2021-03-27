"""
Slice model cuts a regular grid from the SPECFEM .bin files up to the highest
topography point. Consequently, values in the air are interpolated from surface 
points which is unrealistic. This script will set all points above topography
and bathymetry values to a specified value.
"""
import os
import sys
for p in ["/scale_wlg_persistent/filesets/project/nesi00263/PyPackages/"
          "simutils/meshing/topography",
          "/scale_wlg_persistent/filesets/project/nesi00263/PyPackages/"
          "simutils/meshing"]:
    sys.path.append(p)
import numpy as np
from checkerboardiphy import xyz_reader, parse_data_to_header, write_xyz
from generate_topo import (read_nc, convert_coordinates, cut_topography, 
                           interpolate_points)


work_dir = "/home/chowbr/tomo/subduction/set_values_above_topo"
xyz_fid = os.path.join(work_dir, "tomography_model_shallow.xyz")
topo_fid = "/home/chowbr/primer/auxiliary/topo/e140s10.nc"
fid_out = os.path.join(work_dir, "tomography_model_shallow_topo_zero.xyz")

h, d = xyz_reader(xyz_fid, save=True)

xmin = float(h["orig_x"])
xmax = float(h["end_x"])
ymin = float(h["orig_y"])
ymax = float(h["end_y"])
dx = float(h["spacing_x"])

# Open, cut and interpolate the topography to the correct dimensions
topo_utm = convert_coordinates(data=read_nc(topo_fid), utm_projection=-60)
topo_cut = cut_topography(data=topo_utm, x_or_lon_min=xmin - dx, 
                          x_or_lon_max=xmax + dx, y_or_lat_min=ymin - dx, 
                          y_or_lat_max=ymax + dx
                          )
topo_intp = interpolate_points(data=topo_cut, x_min=xmin, x_max=xmax, 
                               y_min=ymin, y_max=ymax, spacing_m=dx
                               )

import ipdb;ipdb.set_trace()
# Use the topography to determine values in the velocity model in air 
data_out = d.copy()
x, y, z = topo_intp.T
indices = None
with open("indices.txt", "a") as f:
    for i, (x_, y_, z_) in enumerate(zip(x, y, z)):
        idx = np.where((data_out[:, 0] == x_) & 
                       (data_out[:, 1] == y_) & 
                       (data_out[:, 2] > z_))[0]
        for idx_ in idx:
            f.write(f"{idx_}\n")

        if indices is None:
            indices = idx
        else:
            indices = np.concatenate((indices, idx))
        print(f"{i}/{len(x)}")

indices = np.unique(indices)
print(f"{len(indices)} / {len(x_)} values zerod")
data_out[indices, 3:] *= 0


print(f"{c} / {len(x)} values changed")
head_out = parse_data_to_header(data_out)
write_xyz(head_out, data_out, fid_out)






