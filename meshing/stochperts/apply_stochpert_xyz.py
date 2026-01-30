"""
Apply stochastic perturbations to an .xyz velocity  model

NEXT STEPS
write only the top layer and allow changing dz parameters and making new models
"""
import numpy as np
from scipy.fft import fftn, ifftn
from numpy.fft import fftfreq, fftshift, ifftshift


def cosine_taper_axis(data, axis=0, left_frac=0.05, right_frac=0.05):
    """
    Apply a 1-D raised-cosine taper along `axis` of `data` (ndarray).
    left_frac, right_frac: fraction (0..0.5) of axis length to taper at each end.
    Returns tapered array (same shape, dtype preserved where possible).

    Written by GPT-5 mini
    """
    data = np.asarray(data)
    n = data.shape[axis]
    taper = np.ones(n, dtype=float)

    # left taper
    l = int(np.floor(left_frac * n))
    if l > 0:
        theta = np.linspace(np.pi, 0.0, l, endpoint=True)
        taper[:l] = 0.5 * (1.0 + np.cos(theta))

    # right taper
    r = int(np.floor(right_frac * n))
    if r > 0:
        theta = np.linspace(0.0, np.pi, r, endpoint=True)
        taper[-r:] = 0.5 * (1.0 + np.cos(theta))

    # broadcast taper to data shape and apply
    shape = [1] * data.ndim
    shape[axis] = n
    taper_reshaped = taper.reshape(shape)

    return data * taper_reshaped


def make_perturbation(x, y, z, seed=123, a=1E3, hv_ratio=1, mean_vel=1.,
                      std_vel=0.1, nmin=-.1, nmax=.1,
                      left_frac=0.4, right_frac=0.0,
                      indexing="ij", plot=False):
    """
    Generate stochastic perturbation array based on unique indexing values

    .. note::

        For some reason the cosine taper applied to the bottom was shifting
        the min value by a small percentage (-1 to -.956). I just fudged it by 
        changing the original min value, that got the bounds back to [-1, 1]
        but if you change `left_frac` you may need to change `nmin`

    :param left_frac: fraction to start tapering off the bottom of the model
    :param right_frac: fraction to start tapering off the top of the model
    """
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    dz = z[1] - z[0]

    # MAKE PERTURBATION BRICK
    # Create the 3D random velocity field with Gaussian distribution
    [X, Y, Z] = np.meshgrid(x, y, z, indexing=indexing)
    np.random.seed(seed)  # gets applied to random call below
    S = np.random.normal(mean_vel, std_vel, Z.shape)

    # 3D FFT into wave number domain. Shift for zero freq. at the center
    FS = fftshift(fftn(S))

    # Generate the frequency domain so that we can use it to index the 
    # perturbation. Shift it so that zero frequency is at the center, so 
    # that when we apply the perturbation it acts on the correct axes
    kx = fftshift(fftfreq(len(x), d=dx))
    ky = fftshift(fftfreq(len(y), d=dy))
    kz = fftshift(fftfreq(len(z), d=dz))
    [KX, KY, KZ] = np.meshgrid(kx, ky, kz, indexing=indexing)

    # Scale the wavenumber vectors so that the perturbations come out
    # the same size, otherwise they will scale with the length of the 
    # domain which is not what we want
    scale_x = x.max() - x.min()
    scale_y = y.max() - y.min()
    scale_z = z.max() - z.min()
    scale_factor = max([scale_x, scale_y, scale_z])

    # Also allow horizontal-to-vertical scaling
    KX *= hv_ratio * (scale_x / scale_factor)
    KY *= hv_ratio * (scale_y / scale_factor)
    KZ *= (scale_z / scale_factor)

    # Generate perburations as a function of the wavenumber domain 
    # (Table 1 [1])
    # Normalizations provide appropriate scaling (Appendix B [1])
    KR = (KX**2 + KY**2 + KZ**2) ** 0.5  # radial wavenumber

    # Von Karman spectral filter in wavenumber domain
    perturbation = a**2 / (1 + KR**2 * a**2)

    # Multiply the Gaussian with the RNG in the wavenumber domain 
    FSP = FS * perturbation 

    # Inverse Fourier transform to get back to space domain
    FSPS = ifftshift(FSP) # Shift perturbed wavenumber spectra back
    S_pert = ifftn(FSPS)  # Inverse FFT to space domain

    # Normalize the final array from a to b 
    arr = np.abs(S_pert)
    S_pert = \
        ((nmax - nmin) * (arr-arr.min()) / (arr.max()-arr.min())) + nmin

    # Apply a 1-D cosine taper along the Z dimension (written by GPT-5 mini)
    if left_frac or right_frac:
        S_pert = cosine_taper_axis(S_pert, axis=2,  # z-axis
                                   left_frac=left_frac, right_frac=right_frac)

    if plot:
        import pyvista as pv
        data = pv.wrap(S_pert)
        data.plot(volume=True, cmap="seismic")

    return S_pert.ravel()


def main(path=None):
    """
    Main function to generate a 3D tomography model with stochastic perturbation
    and export it to a format compatible with SPECFEM3D_Cartesian.
    """
    with open(path, "r") as f:
        header = f.readlines()[:4]
    print(f"HEADER\n {''.join(header)}")

     # Load and distribute the data
    xyz = np.loadtxt(path, skiprows=4)
    if xyz.shape[1] == 6:  # only velocity
        x, y, z, vp, vs, rho = xyz.T
        qmu, qkappa = None, None
        q_bool = False
    else:  # velocity and attenuation
        x, y, z, vp, vs, rho, qmu, qkappa = xyz.T
        q_bool = True

    xvals = np.unique(x)
    yvals = np.unique(y)
    zvals = np.unique(z)

    perturbation = make_perturbation(x=xvals, y=yvals, z=zvals)

    # Perturb the model values
    vp = vp + (vp * perturbation)
    vs = vs + (vs * perturbation)
    rho = rho + (rho * perturbation)
    if q_bool:
        qmu = qmu + (qmu * perturbation)
        qkappa = qkappa + (qkappa * perturbation)

    with open("tomography_model_perturbed.xyz", "w") as f:
        # Header - min and max range values
        f.write(
            f"{xvals.min():.1f} {yvals.min():.1f} {zvals.min():.1f} "
            f"{xvals.max():.1f} {yvals.max():.1f} {zvals.max():.1f}\n"
            )

        # Header - spacing values
        f.write(f"{xvals[1]-xvals[0]:.1f} {yvals[1]-yvals[0]:.1f} "
                f"{zvals[1]-zvals[0]:.1f}\n")

        # Header - number of grid points
        f.write(
            f"{int(len(xvals)):d} {int(len(yvals)):d} {int(len(zvals)):d}\n")

        # Header - parameter min max values. Take either the perturbed 
        # version (preferred) or the 1D model value (if not perturbed)
        # !!! These cap the actual model values in SPECFEM so its important
        # !!! that we get them right
        for arr in [vp, vs, rho]:
            f.write(f"{arr.min():.1f} {arr.max():.1f} ")

        # Note that Q does not need to be included in the header
        f.write("\n")

        # Write actual values
        for i, (x_, y_, z_, vp_, vs_, rho_) in \
                            enumerate(zip(x, y, z, vp, vs, rho)):
            f.write(f"{x_:.1f} {y_:.1f} {z_:.1f} "
                    f"{vp_:.1f} {vs_:.1f} {rho_:.1f} ")
            # Optional Q values
            if q_bool:
                f.write(f"{qmu[i]:.1f} {qkappa[i]:.1f}")
            f.write("\n")


if __name__ == "__main__":
    main("tomography_model.xyz")
