"""
Test out 2D stochastic perturbation following Frankel 1968

.. Parameters::

    a (int): Correlation distance with relevant characteristics at ka > 1
             where k is the wave-number


"""
import sys
import matplotlib.pyplot as plt
import numpy as np

from scipy.fft import fftn, ifftn
from numpy.fft import fftfreq, fftshift, ifftshift


# PARAMETERS
try:
    choice = sys.argv[1]
except IndexError:
    choice = "vonkarman"  # perturbation type [gaussian, exponential, vonkarman]
a = 80  # correlation distance [m]
seed = 12345  # RNG seed to make sure it's consistent
vel_central = 6  # central velocity value [m/s]
vel_std = 0.1  # std. deviation from central velocity
dx = .1  # grid spacing x
dy = .1  # grid spacing y
x = np.arange(0, 210, dx)  # x axis range
y = np.arange(0, 200, dy)  # y axis range
cmap_space = "seismic_r"  # colormap for plots
cmap_wavenumber = "viridis"  # colormap for plots
plot = [4]  # plots to make, 1: space, 2: wavenumber, 3: perturbation, 4: final


# 2D random numbers in space
[X, Y] = np.meshgrid(x, y)
rng = np.random.default_rng(seed)
min_vel = vel_central - (vel_central * vel_std)
max_vel = vel_central + (vel_central * vel_std)
S = rng.uniform(min_vel, max_vel, X.shape)

# Plot the space domain with RNG velocity
if 1 in plot:
    plt.imshow(S, cmap=cmap_space, extent=[x.min(), x.max(), y.min(), y.max()])
    plt.colorbar()
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title(f"Original RNG (dx={dx}, dy={dy})")
    plt.show()

# 2D FFT into wave number domain. Shift to get zero frequency at the center
FS = fftshift(fftn(S))

# Generate the frequency domain so that we can use it to index the perturbation
# Shift it so that zero frequency is at the center, so that when we apply the
# perturbation it acts on the correct axes
kx = fftshift(fftfreq(len(x), d=dx))
ky = fftshift(fftfreq(len(y), d=dy))
[KX, KY] = np.meshgrid(kx, ky)

# Plot the wavenumber domain with zero frequency at the center. Only plot the
# real part of the wavenumber domain
if 2 in plot:
    # Wrap the frequency domain so that 0 is in the middle
    plt.imshow(np.log(np.abs(FS**2)), cmap=cmap_wavenumber, 
               extent=[kx.min(), kx.max(), ky.min(), ky.max()])
    plt.colorbar()
    plt.xlabel("K_x")
    plt.ylabel("K_y")
    plt.title(f"2D FFT of RNG; wavenumber domain")
    plt.show()

# Generate perburations as a function of the wavenumber domain
# See reference Table 1 from Frankel and Clayton 1986
KR = (KX**2 + KY**2) ** 0.5  # radial wavenumber; k_r = (k_x**2 + k_y**2)**(1/2)
if choice == "gaussian":
    perturbation = (a**2 / 2) * np.exp(-1 * KR**2 * a**2 / 4) 
elif choice == "exponential": 
    perturbation = a**2 / (1 + KR**2  * a**2) ** (3/2)
elif choice == "vonkarman":
    perturbation = a**2 / (1 + KR**2 * a**2)

# Multiply the Gaussian with the RNG in the wavefnumber domain 
FS *= perturbation

if 3 in plot:
    plt.imshow(np.log(np.abs(FS**2)), cmap=cmap_wavenumber, 
               extent=[kx.min(), kx.max(), ky.min(), ky.max()])
    plt.colorbar()
    plt.xlabel("K_x")
    plt.ylabel("K_y")
    plt.title(f"{choice} perturbation, wavenumber domain")
    plt.show()

# Inverse Fourier transform to get back to space domain
FS = ifftshift(FS) # Shift perturbed wavenumber spectra back
S_pert = ifftn(FS)  # Inverse FFT to space domain

if 4 in plot:
    plt.imshow(np.log(np.abs(S_pert**2)), cmap=cmap_space, 
               extent=[x.min(), x.max(), y.min(), y.max()])
    plt.xlabel(f"X (dx={dx})")
    plt.ylabel(f"Y (dy={dy})")
    plt.title(f"{choice.capitalize()} Perturbation (a={a})")
    plt.colorbar()
    plt.show()



