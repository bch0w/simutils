"""
Apply 3D stochastic perturbations following Frankel and Clayton 1986 [1]. 
Creates a 2D Gaussian distributed RNG velocity field, Fourier transforms into
the wavenumber domain, applies a choice of spectral filter, and then inverse
transforms. The final velocity field contains stochastic perturbations based on
the initial velocity field and chosen filter.

.. References::

    [1] Frankel, Arthur, and Robert W. Clayton. "Finite difference simulations 
        of seismic scattering: Implications for the propagation of short‐period 
        seismic waves in the crust and models of crustal heterogeneity."
        Journal of Geophysical Research: Solid Earth 91.B6 (1986): 6465-6489.
    [2] https://discourse.vtk.org/t0/\
            numpy-3d-array-into-vtk-data-types-for-volume-rendering/3455/3 
"""
import sys
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

from scipy.fft import fftn, ifftn
from numpy.fft import fftfreq, fftshift, ifftshift


# ==============================================================================
#                               PARAMETERS
# ==============================================================================
# PERTURBATION
a = 80  # km
norm_min = -1
norm_max = 1

# INITIAL CONDITIONS
seed = 123  
mean_vel = 1  # km/s
std_vel = 0.1  

# GRID SPACING [km]
dx = 1  
dy = 1 
dz = 1
xmin = 0
xmax = 210
ymin = 0
ymax = 200
zmin = 0
zmax = 200

# PLOTTING
plot = True
cmap = "seismic_r"
# ==============================================================================

# Allow User cmdline overwrite of choice
try:
    choice = sys.argv[1]
except IndexError:
    pass

# Define the space domain
x = np.arange(xmin, xmax, dx)  # x axis range
y = np.arange(ymin, ymax, dy)  # y axis range
z = np.arange(zmin, zmax, dz)  # y axis range
[X, Y, Z] = np.meshgrid(x, y, z, indexing="ij")

# Create the 2D random velocity field with Gaussian distribution
np.random.seed(seed)
S = np.random.normal(mean_vel, std_vel, X.shape)

# 3D FFT into wave number domain. Shift to get zero frequency at the center
FS = fftshift(fftn(S))

# Generate the frequency domain so that we can use it to index the perturbation
# Shift it so that zero frequency is at the center, so that when we apply the
# perturbation it acts on the correct axes
kx = fftshift(fftfreq(len(x), d=dx))
ky = fftshift(fftfreq(len(y), d=dy))
kz = fftshift(fftfreq(len(z), d=dz))
[KX, KY, KZ] = np.meshgrid(kx, ky, kz, indexing="ij")

# Generate perburations as a function of the wavenumber domain (Table 1 [1])
# Normalizations provide appropriate scaling (Appendix B [1])
KR = (KX**2 + KY**2 + KZ**2) ** 0.5  # radial wavenumber

# Von Karman spectral filter
perturbation = a**2 / (1 + KR**2 * a**2)

# Multiply the Gaussian with the RNG in the wavenumber domain 
FSP = FS * perturbation 

# Inverse Fourier transform to get back to space domain
FSPS = ifftshift(FSP) # Shift perturbed wavenumber spectra back
S_pert = ifftn(FSPS)  # Inverse FFT to space domain

# Normalize the final array from a to b so that it can be used for perturbations
nmin = -1
nmax = 1
arr = np.abs(S_pert)
S_pert = ((nmax - nmin) * (arr-arr.min()) / (arr.max()-arr.min())) + nmin

# Plot volumetric cube with PyVista
data = pv.wrap(S_pert)
data.plot(volume=True, cmap="seismic_r")
