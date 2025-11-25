# Plot Venus Topo
import pyshtools as pysh
import matplotlib.pyplot as plt

pysh.utils.figstyle(rel_width=0.75)
clm = pysh.datasets.Venus.VenusTopo719() / 1E3  # unit km
clm.coeffs[0, 0, 0] = 0.
clm.coeffs[0, 2, 0] = 0.
grid = clm.expand()
fig, ax = grid.plot(colorbar="bottom", cb_label="Elevation [km]", cmap="viridis", cmap_limits=[-4,4, 0.1])
