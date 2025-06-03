"""
Apply 2D stochastic perturbations following Frankel and Clayton 1986 [1]. 
Creates a 2D Gaussian distributed RNG velocity field, Fourier transforms into
the wavenumber domain, applies a choice of spectral filter, and then inverse
transforms. The final velocity field contains stochastic perturbations based on
the initial velocity field and chosen filter.

.. Parameters::

    choice (str): spectral filter to use, 'gaussian', 'exponential', 'vonkarman'
    seed (int): seed for the random number generator
    mean_vel (float): mean velocity for velocity seed
    std_vel (float): standard deviation for velocity seed

    a (int): Correlation distance with relevant characteristics at ka > 1 
        where k is the wave-number
    sigma_c (float): standard deviation for spectral filter normalization
    norm_min (float): lower bound for final normalization
    norm_max (float): upper bound for final normalization

    dx (float): grid spacing x
    dy (float): grid spacing y
    xmin (float): min x-axis value
    xmax (float): max x-axis value
    ymin (float): min y-axis value
    ymax (float): max y-axis value

    plot_summary (bool): plot a 4-panel figure showing the creation process
    plot_final (bool): plot the final velocity model only
    cmap_space (str): matplotlib colormap for space domain
    cmap_wavenumber (str): matplotlib colormap for wavenumber domain

.. References::

    [1] Frankel, Arthur, and Robert W. Clayton. "Finite difference simulations 
        of seismic scattering: Implications for the propagation of short‐period 
        seismic waves in the crust and models of crustal heterogeneity."
        Journal of Geophysical Research: Solid Earth 91.B6 (1986): 6465-6489.
    [2] https://stackoverflow.com/questions/75903584/\
            how-to-construct-a-vtk-file-with-point-data-from-numpy-arrays
"""
import sys
import matplotlib.pyplot as plt
import numpy as np
import vtk

from scipy.fft import fftn, ifftn
from numpy.fft import fftfreq, fftshift, ifftshift

from vtk.numpy_interface import algorithms as algs
from vtk.numpy_interface import dataset_adapter as dsa


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
[X, Y, Z] = np.meshgrid(x, y, z)

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
[KX, KY, KZ] = np.meshgrid(kx, ky, kz)

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

# Export the data array as VTK so we can view in ParaView [2]
polydata = vtk.vtkPolyData()
pts = vtk.vtkPoints()
points = algs.make_vector(X.ravel(), Y.ravel(), Z.ravel())
pts.SetData(dsa.numpyTovtkDataArray(points, "Points"))
polydata.SetPoints(pts)
vectors = algs.make_vector(S_pert[:, :, 0].ravel(),
                           S_pert[:, :, 1].ravel(),
                           S_pert[:, :, 2].ravel()
                           )
breakpoint()
polydata.GetPointData().SetScalars(dsa.numpyTovtkDataArray(vectors, "Velocity"))

# Write the VTK file
writer = vtk.vtkPolyDataWriter()
writer.SetFileName("test.vtk")
writer.SetInputData(polydata)
writer.Update()
