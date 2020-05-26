"""
Specfem3D Cartesian uses an ellipsoidal/spherical Gaussian smoothing function
to remove short-wavelength features in models/kernels. Inputs for the smoothing
function are sigma_h and sigma_v which correspond to the standard deviation for 
the Gaussian function, which is defined in 1-D as:

    G(x) = (2*pi * o**2)**(-1/2) * (e ** (-x**2/(2*o**2)))

where o is used in place of sigma, the standard deviation or half-width of 
the Gaussian. Therefore the actual search area, or scalelength of the Gaussian 
will be larger than the given values of sigma_h and sigma_v.

From Tape et al. 2007, the following notes:
    -Select full-width (gamma) close or slightly smaller than the shortest
    wavelength of seismic waves used
    -
"""
import numpy as np


def convert_sigma_to_scalelength(sigma):
    """
    Scale length approximantely S ~ sigma * sqrt(8.0) for Gaussian smoothing,
    from Specfem3D smooth_sem.F90
    :type sigma: float
    :param sigma: standard deviation of the Gaussian operator
    :rtype: float
    :return: the approximate scale length of the Gaussian operator
    """
    return sigma * np.sqrt(8.)


def gaussian(x, sigma):
    """
    Return a Gaussian function
    """
    return (2 * np.pi * sigma**2)**(-1/2) * (np.e**(-x**2/2 * sigma**2))
    

