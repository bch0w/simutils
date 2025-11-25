"""
Extract Venus Topo to equally spaced latlon grid
"""
import numpy as np
import pyshtools as pysh
import matplotlib.pyplot as plt

clm = pysh.datasets.Venus.VenusTopo719() / 1E3  # unit km
# grid = clm.expand(grid="DH")  # DH: equal lat lon, DH2: double sample in lon, GLQ: gauss-legendre quadrature
# coeffs = grid.to_array()

# set param
longitude = 180
dlat = 1.

# line of constant longitude
lat = np.arange(-90., 90., 1)
lon = np.ones(len(lat)) * longitude

# extract values
vals = clm.expand(lat=lat, lon=lon)

# plot
plt.plot(lat, vals, "ko-", markersize=1)
plt.xlabel("Latitude (north ->)")
plt.ylabel("Elevation [km]")
plt.title(f"VenusTopo719 (lon={longitude})")
plt.grid()

# plot the full dataset for comparison
clm.coeffs[0, 0, 0] = 0.
clm.coeffs[0, 2, 0] = 0.
grid = clm.expand()
fig, ax = grid.plot(colorbar="bottom", cb_label="Elevation [km]", 
                    cmap="viridis", cmap_limits=[-4,4, 0.1])

plt.show()


