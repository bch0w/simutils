"""
Generate a 1D velocity model (e.g., PREM) as a 3D external tomography model
exportable to SPECFEM3D_Cartesian, then apply 3D stochastic perturbations
using a Von Karman spectral filter, creating a 1D velocity model with stochastic
3D perturbations that mimic small-scale heterogeneities and allow for 
high-frequency scattering in spectral element simulations


NEXT STEPS
write only the top layer and allow changing dz parameters and making new models
"""
import os
import matplotlib.pyplot as plt
import numpy as np
try:
    import pyvista as pv
except ImportError:
    print("PyVista is not installed, skipping 3D plotting")
    pass

from scipy.fft import fftn, ifftn
from numpy.fft import fftfreq, fftshift, ifftshift

# Bad habit but import all variables from the config file
from gptmcfg import *

# Define 1D PREM that can be interpolated and perturbed
MODELS = {
    "PREM": {
        "depth": 1E3 * np.array([
            -1.0, 3.0, 15., 24.4, 71., 80., 171., 220., 271., 371., 400.
        ]),  # km
        "vp": 1E3 * np.array([
            1.45, 5.80, 6.80, 8.11, 8.08, 8.08, 8.02, 7.99, 8.56, 8.66, 8.85
        ]),  # km/s
        "vs": 1E3 * np.array([
            1.0, 3.20, 3.90, 4.49, 4.47, 4.47, 4.44, 4.42, 4.62, 4.68, 4.75
        ]),  # km/s
        "rho": 1E3 * np.array([
            1.02, 2.6, 2.9, 3.38, 3.37, 3.37, 3.36, 3.36, 3.44, 3.47, 3.53
        ]),  # kg/m^3
        "qmu": np.array([
            0., 600., 600., 600., 600., 600., 80., 80., 143., 143., 143.
        ]),
        "qkappa": np.array([57323.] * 11)
        }
    }


def interp_1D_model(model, dz, zmin=None, zmax=None):
    """
    1D interpolation of the 1D model between major depth values to get gradients 
    in between rather than just step functions.
    """
    z = model["depth"]  # Convert to meters with positive up
    
    # User can set the bounds of the 1D model but if not just use entire model
    if zmin is None:
        zmin = z.min()
    if zmax is None:
        zmax = z.max()

    # Define the new depth values to intepolate against
    zs = np.arange(zmin, zmax + dz, dz)
    print(f"{len(zs)} total values along a 1D depth profile")

    # Create a new dictionary to store the interpolated values
    model_out = {}
    model_out["depth"] = zs

    # Interpolate each of the other values
    for key in model.keys():
        if key == "depth":
            continue
        y = model[key] 
        yinterp = np.interp(zs, z, y)
        model_out[key] = yinterp
    
    return model_out

# ==============================================================================
# Choose the Model
model = MODELS[model_choice]

# Paths (will be modified by 'tag')
fig_path = "./figures"
mod_path = "./created"

# Simple numbering scheme for tomography files
fids = [f"tomography_model_{_:0>2}.xyz" for _ in range(1, len(ZVALS)+1)]

# Misc parameters that we don't need to set
indexing = "ij"
order = "F"

# Plotting parameters
plot_all = None
show = False
plot_cube = False     # model
plot_brick = False   # perturbation
plot_profile = False  # 1D vel. profile
plot_2d = True       # 2d cross-sections
cmap = "viridis"
# ==============================================================================

# SETUP PROCEDURES
# Modify the tag and add final tag to paths
if perturbations:
    tag = tag+"_wpert"
else:
    tag = tag+"_nopert"

# Modify the paths to include the tag subdir
fig_path = os.path.join(fig_path, tag)
mod_path = os.path.join(mod_path, tag)

# Make directories
for path_ in [fig_path, mod_path]:
    if not os.path.exists(path_):
        os.makedirs(path_)

fids = [os.path.join(mod_path, _) for _ in fids]

print(f"Tag is: {tag}")
print(f"Figures will be saved to: {fig_path}")
print(f"Models will be saved to: {mod_path}\n")

# Overwrite Flags
if plot_all is not None:
    if plot_all:
        plot_cube, plot_brick, plot_profile, plot_2d = True, True, True, True
    else:
        plot_cube, plot_brick, plot_profile, plot_2d = False, False, False, False

# Drop attenuation if needed
if not include_q:
    print("droping 'Q' values from model")
    MODEL = {key: val for key, val in model.items() if not key.startswith("q")}

# BEGIN PROCESSING HERE
for l, zvals in enumerate(ZVALS):
    zmin, zmax = zvals
    print(f"{fids[l]} from {zmin} to {zmax}")

    # Assign the grid spacing and lower value
    dx, dy, dz = DX[l], DY[l], DZ[l]

    # Define the brick in which the perturbation exists
    x = np.arange(xmin, xmax + dx, dx)  # x axis range
    y = np.arange(ymin, ymax + dx, dy)  # y axis range
    z = np.arange(zmin, zmax + dz, dz)  # z axis range
    z_pert = np.arange(zmin_pert, zmax_pert, dz)  # z axis range, perturbation
    
    # Check if the perturbation exists in this layer
    pert_idx_start, pert_idx_end = None, None
    if perturbations:
        # Determine where to start the perturbation brick, either the brick 
        # lives inside the layer, or it starts at the top of the layer
        if zmin_pert <= zmax:
           try:
               # Brick lives inside the layer
               pert_idx_start = np.where(z == zmin_pert)[0][0]
           except IndexError:
               # Brick starts at the top of the layer
               pert_idx_start = np.where(z == zmin)[0][0]
        if zmax_pert >= zmin:
           try:
               # Brick lives inside the layer
               pert_idx_end = np.where(z == zmax_pert)[0][0]
           except IndexError:
               # Brick starts at the top of the layer
               pert_idx_end = np.where(z == zmax)[0][0]
    
    if pert_idx_start is None:
        print("No perturbations defined for this layer...")

    # Generate the perturbation if it exists in this layer
    if pert_idx_start is not None:
        print(f"generating perturbation brick from {z[pert_idx_start]} to "
              f"{z[pert_idx_end]}")
        # Figure out how to shift the indices of the perturbation brick so that
        # we can use the same indexing of the depth model to access perturbation

        # Print some information
        print(f"dX = {(xmax - xmin) * 1E-3:.2f} km\n"
              f"dY = {(ymax - ymin) * 1E-3:.2f} km\n"
              f"dZ = {(zmax - zmin) * 1E-3:.2f} km\n"
              f"dZ_pert = {(zmax_pert - zmin_pert) * 1E-3:.2f} km\n"
              "\n"
              f"nX = {len(x)}\n"
              f"nY = {len(y)}\n"
              f"nZ_pert = {len(z_pert)}\n"
              )

        # MAKE PERTURBATION BRICK
        # Create the 3D random velocity field with Gaussian distribution
        [X, Y, Z_PERT] = np.meshgrid(x, y, z_pert, indexing=indexing)
        np.random.seed(seed)
        S = np.random.normal(mean_vel, std_vel, Z_PERT.shape)

        # 3D FFT into wave number domain. Shift for zero frequency at the center
        FS = fftshift(fftn(S))

        # Generate the frequency domain so that we can use it to index the 
        # perturbation. Shift it so that zero frequency is at the center, so 
        # that when we apply the perturbation it acts on the correct axes
        kx = fftshift(fftfreq(len(x), d=dx))
        ky = fftshift(fftfreq(len(y), d=dy))
        kz = fftshift(fftfreq(len(z_pert), d=dz))
        [KX, KY, KZ] = np.meshgrid(kx, ky, kz, indexing=indexing)

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
        S_pert = ((nmax - nmin) * (arr-arr.min()) / (arr.max()-arr.min())) + nmin

        # Plot volumetric cube with PyVista
        if plot_brick:
            data = pv.wrap(S_pert)
            data.plot(volume=True, cmap="seismic_r")
            plt.close("all")

    # GENERATE VELOCITY MODEL
    # Now we generate the velocity model cube. First interpolate to our `dz`
    [X, Y, Z] = np.meshgrid(x, y, z, indexing=indexing)
    print(f"layer '{fids[l]}' has {np.prod(X.shape)} points")

    # Interpolate the 1D model to the required DZ
    model = interp_1D_model(model=model, dz=dz, zmin=zmin, zmax=zmax)

    if plot_profile:
        # Just plot one example parameter 
        f, ax = plt.subplots(figsize=(6,8))
        plt.plot(model["vs"], model["depth"], "ko-")
        plt.gca().invert_yaxis()
        plt.xlabel("Vs [m/s]")
        plt.ylabel("Depth [m]")
        plt.title(f"Interpolated 1D Model {model_choice}\n{fids[l]}")
        plt.grid()
        plt.savefig(os.path.join(fig_path, f"1d_profile_{zmin}-{zmax}.png"))
        if show:
            plt.show()
        plt.close("all")

    # Build model for each of the model parameters. Add perturbations if needed
    model_dict = {}
    for parameter in model.keys():
        # Ignore depth values, 
        if parameter == "depth":
            continue
        
        print(f"making model for par: {parameter}")
        arr = model[parameter]

        # Generate the 3D cube shape
        cube = np.ones((len(x), len(y), len(arr)))

        # Fill in the cube by multiplying the correct depth with their assigned
        # value. This is pretty brute force but it works
        for i, val in enumerate(arr):
            # Apply perturbation ontop of cube if we are in the correct layer
            if (parameter in perturb) and \
                 pert_idx_start is not None and \
                    (pert_idx_start <= i < pert_idx_end): 
                j = i - pert_idx_start  # set correct index for perturbation
                cube[:,:,i] *= (val + val * S_pert[:,:,j])
            # If not, then we just assign the correct value
            else:
                cube[:,:,i] *= val  
        
        # Store the unraveled 1D array for writing
        # Flip the values because we will flip all the spatial coordinates so 
        # that deepest points are listed first in the XYZ model
        model_dict[parameter] = np.flip(cube.ravel(order=order))

        # Plot test depth value and parameter
        if plot_2d:
            im = plt.imshow(cube[:,:,1], cmap=cmap, 
                            extent=np.array([xmin, xmax, ymin, ymax]) * 1E-3)
            plt.colorbar(im, shrink=0.8, pad=0.025, label=parameter)
            plt.grid()
            plt.xlabel("X [km]")
            plt.ylabel("Y [km]")
            ax = plt.gca()
            # ax.ticklabel_format(style="sci", axis="both", scilimits=(0,0))
            ax.ticklabel_format(style="plain", axis="both")
            ax.set_aspect("equal")
            plt.title(f"Vs [m/s] at {arr[1]}m")
            plt.tight_layout()
            plt.savefig(os.path.join(fig_path, 
                                     f"2d_plot_{parameter}_{zmin}-{zmax}.png")
                                     )
            if show:
                plt.show()
            plt.close("all")

    # Plot volumetric cube with PyVista
    if plot_cube:
        data = pv.wrap(cube[:,:,::-1])  # Flip the cube so Z positive is down
        data.plot(volume=True, cmap="rainbow_r")
        plt.close("all")

    # EXPORT MODEL
    # Ravel turns the gridded data into 1D arrays that can be written
    # We flip the arrays so that when written, the largest negative value starts
    # Z axis is multiplied by -1 so that positive Z is defined as up
    x_list = np.flip(X.ravel(order=order))
    y_list = np.flip(Y.ravel(order=order))
    z_list = np.flip(-1 * Z.ravel(order=order))

    with open(fids[l], "w") as f:
        # Header - min and max range values
        f.write(f"{x_list.min():.1f} {y_list.min():.1f} {z_list.min():.1f} "
                f"{x_list.max():.1f} {y_list.max():.1f} {z_list.max():.1f}\n")

        # Header - spacing values
        f.write(f"{dx:.1f} {dy:.1f} {dz:.1f}\n")

        # Header - number of grid points
        f.write(f"{int(len(x)):d} {int(len(y)):d} {int(len(z)):d}\n")

        # Header - parameter min max values
        f.write(f"{model['vp'].min():.1f} {model['vp'].max():.1f} "
                f"{model['vs'].min():.1f} {model['vs'].max():.1f} "
                f"{model['rho'].min():.1f} {model['rho'].max():.1f} ")
        # Note that Q does not need to be included in the header
        f.write("\n")

        for i in range(len(x_list)):
            f.write(f"{x_list[i]:.1f} {y_list[i]:.1f} {z_list[i]:.1f} "
                    f"{model_dict['vp'][i]:.1f} {model_dict['vs'][i]:.1f} "
                    f"{model_dict['rho'][i]:.1f} ")
            # Optional Q values
            if include_q:
                f.write(f"{model_dict['qmu'][i]:.1f} "
                        f"{model_dict['qkappa'][i]:.1f}")
            f.write("\n")
        
