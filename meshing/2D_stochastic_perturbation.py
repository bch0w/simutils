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
"""
import sys
import matplotlib.pyplot as plt
import numpy as np

from scipy.fft import fftn, ifftn
from numpy.fft import fftfreq, fftshift, ifftshift


# ==============================================================================
#                               PARAMETERS
# ==============================================================================
# PERTURBATION
choice = "vonkarman"  # gaussian, exponential, vonkarman
a = 160  # km
sigma_c = 0.1  
norm_min = -1
norm_max = 1

# INITIAL CONDITIONS
seed = 123  
mean_vel = 1  # km/s
std_vel = 0.1  
distribution = "normal"  # random, normal

# GRID SPACING [km]
dx = .1  
dy = .1 
xmin = 0
xmax = 210
ymin = 0
ymax = 200

# PLOTTING
plot_summary = False
plot_final = True
cmap_space = "seismic_r"
cmap_wavenumber = "viridis"
# ==============================================================================

# Allow User cmdline overwrite of choice
try:
    choice = sys.argv[1]
except IndexError:
    pass

# Define the space domain
x = np.arange(xmin, xmax, dx)  # x axis range
y = np.arange(ymin, ymax, dy)  # y axis range
[X, Y] = np.meshgrid(x, y)

# Create the 2D random velocity field with Gaussian distribution
if distribution == "random":
    rng = np.random.default_rng(seed=seed)
    min_vel = mean_vel - (mean_vel * std_vel)
    max_vel = mean_vel + (mean_vel * std_vel)
    S = rng.uniform(min_vel, max_vel, X.shape)
elif distribution == "normal":
    np.random.seed(seed)
    S = np.random.normal(mean_vel, std_vel, X.shape)

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
k_nyq = 1 / (2 * min(dx, dy))  # Nyquist wavenumber 

if choice == "gaussian":
    perturbation = (a**2 / 2) * np.exp(-1 * KR**2 * a**2 / 4) 
    normalization = sigma_c ** 2
elif choice == "exponential": 
    perturbation = a**2 / (1 + KR**2  * a**2) ** (3/2)
    normalization = sigma_c ** 2
elif choice == "vonkarman":
    perturbation = a**2 / (1 + KR**2 * a**2)
    normalization = 0.29 * sigma_c  ** 2

# Multiply the Gaussian with the RNG in the wavenumber domain 
FSP = normalization * FS * perturbation 

# Inverse Fourier transform to get back to space domain
FSPS = ifftshift(FSP) # Shift perturbed wavenumber spectra back
S_pert = ifftn(FSPS)  # Inverse FFT to space domain

# Normalize the final array from a to b so that it can be used for perturbations
nmin = -1
nmax = 1
arr = np.abs(S_pert)
S_pert = ((nmax - nmin) * (arr-arr.min()) / (arr.max()-arr.min())) + nmin

# 4-Panel plot to show all the steps of the process
if plot_summary:
    f, axs = plt.subplots(2, 2, figsize=(10, 10))

    # Colorbar control parameters
    shrink = 0.6
    pad = 0.025

    # Original space domain defined by RNG velocities
    im1 = axs[0][0].imshow(S, cmap=cmap_space, 
                     extent=[x.min(), x.max(), y.min(), y.max()])
    axs[0][0].set_title(f"RNG Velocity ('{distribution}' seed={seed})")
    axs[0][0].set_xlabel("X [km]")
    axs[0][0].set_ylabel("Y [km]")
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
    axs[1][0].set_title(f"'{choice.capitalize()}' Wavenumber Filter")
    axs[1][0].set_xlabel("KX")
    axs[1][0].set_ylabel("KY")
    plt.colorbar(im3, ax=axs[1][0], shrink=shrink, pad=pad,
                 label="ln(abs(P(k)^2))")
    
    # Final space domain with stochastic perturbation
    im4 = axs[1][1].imshow(S_pert, cmap=cmap_space, 
                           extent=[x.min(), x.max(), y.min(), y.max()])
    axs[1][1].set_title(f"Normalized Perturbed Model (a={a}km)")
    axs[1][1].set_xlabel("X [km]")
    axs[1][1].set_ylabel("Y [km]")
    plt.colorbar(im4, ax=axs[1][1], shrink=shrink, pad=pad, label="Vp [km/s]")

    f.tight_layout()
    f.suptitle(f"{choice.capitalize()} Stochastic Perturbation", fontsize=16)
    plt.subplots_adjust(wspace=0.2, hspace=-0.2)
    plt.savefig(f"figures/{choice}_filter_a{a}_s{seed}_4panel.png")
    plt.show()

# Plot perturbed velocity model
if plot_final:
    f, ax = plt.subplots(1, figsize=(8, 8))
    im = ax.imshow(S_pert, cmap=cmap_space, 
                   extent=[x.min(), x.max(), y.min(), y.max()])
    plt.title(f"{choice.capitalize()} Perturbation (a={a}km; seed={seed})")
    plt.xlabel(f"X [km] (dx={dx})")
    plt.ylabel(f"Y [km] (dy={dy})")
    plt.colorbar(im, shrink=0.8, pad=0.025, label="Vp [km/s]")
    plt.grid(alpha=0.25)
    
    plt.tight_layout()
    plt.savefig(f"figures/{choice}_filter_a{a}_s{seed}_final.png")
    plt.show()
