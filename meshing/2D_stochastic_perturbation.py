"""
Test out 2D stochastic perturbation following Frankel 1968

.. Parameters::

    a (int): Correlation distance with relevant characteristics at ka > 1
             where k is the wave-number


"""
import matplotlib.pyplot as plt
import numpy as np
import scipy


# PARAMETERS
choice = "gaussian"
a = 25
n_rng = (2, 6)  # min max
dx = .1
dy = .1
x = np.arange(0, 101, dx)
y = np.arange(0, 101, dy)
[X, Y] = np.meshgrid(x, y)
plot = True
cmap = "seismic_r"

# 2D random numbers in space
rng = np.random.default_rng()
S = rng.uniform(n_rng[0], n_rng[1], X.shape)

# 2D FFT into wave number
FS = scipy.fft.fftn(S)

# Generate the frequency domain so that we can use it to index the perturbation
kx = np.fft.fftfreq(len(x), d=dx)
ky = np.fft.fftfreq(len(y), d=dy)
[KX, KY] = np.meshgrid(kx, ky)

# Wrap the frequency domain so that 0 is in the middle
if plot:
    FS_shift = np.log(np.abs(scipy.fft.fftshift(FS)))
    plt.imshow(FS_shift, cmap=cmap, 
               extent=[kx.min(), kx.max(), ky.min(), ky.max()])
    plt.colorbar()
    plt.xlabel("K_x")
    plt.ylabel("K_y")
    plt.title(f"{choice} perturbation, wavenumber domain")
    plt.show()

# Generate the perburation - Reference Table 1 from Frankel and Clayton 1986
if choice == "gaussian":
    # k_r is the radial wavenumber and equal to k_r = (k_x**2 + k_y**2)**(1/2)
    perturbation = (a**2 / 2) * np.exp(-1 * (KX**2 + KY**2) * a**2 / 4) 

# Multiply the Gaussian with the RNG in the wavefnumber domain 
FS *= perturbation

# Inverse transform back to space domain
S_pert = scipy.fft.ifftn(np.abs(FS))
S_pert_abs = np.log(np.abs(scipy.fft.fftshift(S_pert)))  # only take the real

if plot:
    plt.imshow(S_pert_abs, cmap=cmap, 
               extent=[x.min(), x.max(), y.min(), y.max()])
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title("Stochastic Perturbations - Space Domain")
    plt.colorbar()
    plt.show()



