"""
Apply 2D stochastic perturbations following Frankel and Clayton 1986 [1]. 
Creates a 2D Gaussian distributed RNG velocity field, Fourier transforms into
the wavenumber domain, applies a choice of spectral filter, and then inverse
transforms. The final velocity field contains stochastic perturbations based on
the initial velocity field and chosen filter.

.. Parameters::

    choice (str): spectral filter to use, 'gaussian', 'exponential', 'vonkarman'
    seed (int): seed for the random number generator
    vel_central (float): mean velocity for velocity seed
    vel_std (float): standard deviation for velocity seed

    a (int): Correlation distance with relevant characteristics at ka > 1 
        where k is the wave-number
    sigma_c (float): standard deviation for spectral filter normalization

    dx (float): grid spacing x
    dy (float): grid spacing y
    xmin (float): min x-axis value
    xmax (float): max x-axis value
    ymin (float): min y-axis value
    ymax (float): max y-axis value

    plot (bool): plot a 4-panel figure showing the creation process
    cmap_space (str): matplotlib colormap for space domain
    cmap_wavenumber (str): matplotlib colormap for wavenumber domain

.. References::

    [1] Frankel, Arthur, and Robert W. Clayton. "Finite difference simulations 
        of seismic scattering: Implications for the propagation of short‐period 
        seismic waves in the crust and models of crustal heterogeneity."
        Journal of Geophysical Research: Solid Earth 91.B6 (1986): 6465-6489.
"""
import sys
import matplotlib.pyplot as plt
import numpy as np

from scipy.fft import fftn, ifftn
from numpy.fft import fftfreq, fftshift, ifftshift


# ==============================================================================
#                               PARAMETERS
# ==============================================================================
# INITIAL CONDITIONS
choice = "vonkarman"  # gaussian, exponential, vonkarman
seed = 12345  # RNG seed to make sure it's consistent
vel_central = 6  # central velocity value [m/s]
vel_std = 0.2  # std. deviation from central velocity

# PERTURBATION
a = 80  # correlation distance [m]
sigma_c = 0.5  # std. deviation for perturbations

# GRID SPACING
dx = .1  # grid spacing x
dy = .1  # grid spacing y
xmin = 0
xmax = 210
ymin = 0
ymax = 200

# PLOTTING
plot = True
cmap_space = "seismic_r"  # colormap for plots
cmap_wavenumber = "viridis"  # colormap for plots
# ==============================================================================

# Allow User cmdline overwrite of choice
try:
    choice = sys.argv[1]
except IndexError:
    pass

# 2D random numbers in space
x = np.arange(xmin, xmax, dx)  # x axis range
y = np.arange(ymin, ymax, dy)  # y axis range

[X, Y] = np.meshgrid(x, y)
rng = np.random.default_rng(seed)
min_vel = vel_central - (vel_central * vel_std)
max_vel = vel_central + (vel_central * vel_std)
S = rng.uniform(min_vel, max_vel, X.shape)


# 2D FFT into wave number domain. Shift to get zero frequency at the center
FS = fftshift(fftn(S))

# Generate the frequency domain so that we can use it to index the perturbation
# Shift it so that zero frequency is at the center, so that when we apply the
# perturbation it acts on the correct axes
kx = fftshift(fftfreq(len(x), d=dx))
ky = fftshift(fftfreq(len(y), d=dy))
[KX, KY] = np.meshgrid(kx, ky)

# Generate perburations as a function of the wavenumber domain (Table 1 [1])
# Normalizations provide appropriate scaling (Appendix B [1])
KR = (KX**2 + KY**2) ** 0.5  # radial wavenumber; k_r = (k_x**2 + k_y**2)**(1/2)
if choice == "gaussian":
    perturbation = (a**2 / 2) * np.exp(-1 * KR**2 * a**2 / 4) 
    normalization = sigma_c ** 2
elif choice == "exponential": 
    perturbation = a**2 / (1 + KR**2  * a**2) ** (3/2)
    normalization = sigma_c ** 2
elif choice == "vonkarman":
    perturbation = a**2 / (1 + KR**2 * a**2)
    normalization = 0.29 * sigma_c  ** 2

# Multiply the Gaussian with the RNG in the wavefnumber domain 
FSP = FS * normalization * perturbation

# Inverse Fourier transform to get back to space domain
FSPS = ifftshift(FSP) # Shift perturbed wavenumber spectra back
S_pert = ifftn(FSPS)  # Inverse FFT to space domain


# 4-Panel plot to show all the steps of the process
if plot:
    f, axs = plt.subplots(2, 2, figsize=(10, 10))

    # Colorbar control parameters
    shrink = 0.6
    pad = 0.025

    # Original space domain defined by RNG velocities
    im1 = axs[0][0].imshow(S, cmap=cmap_space, 
                     extent=[x.min(), x.max(), y.min(), y.max()])
    axs[0][0].set_title(f"RNG Velocity (seed={seed})")
    axs[0][0].set_xlabel("X")
    axs[0][0].set_ylabel("Y")
    plt.colorbar(im1, ax=axs[0][0], shrink=shrink, pad=pad, label="Vp [km/s]")
    
    # Wavenumber domain without perturbation
    im2 = axs[0][1].imshow(np.log(np.abs(FS**2)), cmap=cmap_wavenumber, 
                     extent=[kx.min(), kx.max(), ky.min(), ky.max()])
    axs[0][1].set_title(f"Wavenumber (FFT velocity)")
    axs[0][1].set_xlabel("KX")
    axs[0][1].set_ylabel("KY")
    plt.colorbar(im2, ax=axs[0][1], shrink=shrink, pad=pad,
                 label="ln(abs(P(k)^2))")
    
    # Wavenumber domain with chosen spectra filter
    im3 = axs[1][0].imshow(np.log(np.abs(FSP**2)), cmap=cmap_wavenumber, 
                     extent=[kx.min(), kx.max(), ky.min(), ky.max()])
    axs[1][0].set_title(f"{choice.capitalize()} Wavenumber Filter")
    axs[1][0].set_xlabel("KX")
    axs[1][0].set_ylabel("KY")
    plt.colorbar(im3, ax=axs[1][0], shrink=shrink, pad=pad,
                 label="ln(abs(P(k)^2))")
    
    # Final space domain with stochastic perturbation
    im4 = axs[1][1].imshow(np.log(np.abs(S_pert**2)), cmap=cmap_space, 
                     extent=[x.min(), x.max(), y.min(), y.max()])
    axs[1][1].set_title(f"Final Velocity Model (a={a}m)")
    axs[1][1].set_xlabel("X")
    axs[1][1].set_ylabel("Y")
    plt.colorbar(im4, ax=axs[1][1], shrink=shrink, pad=pad, label="Vp [km/s]")

    f.tight_layout()
    plt.subplots_adjust(wspace=0.2, hspace=-0.2)
    plt.savefig(f"figures/{choice}_filter_a{a}_s{seed}.png")
    plt.show()
