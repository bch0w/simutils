"""
Generate a 1D velocity model (e.g., PREM) as a 3D external tomography model
exportable to SPECFEM3D_Cartesian, then apply 3D stochastic perturbations
using a Von Karman spectral filter, creating a 1D velocity model with stochastic
3D perturbations that mimic small-scale heterogeneities and allow for 
high-frequency scattering in spectral element simulations


NEXT STEPS
write only the top layer and allow changing dz parameters and making new models
"""
import argparse
import os
import json
import sys
import traceback
import matplotlib.pyplot as plt
import numpy as np
try:
    import pyvista as pv
except ImportError:
    print("PyVista is not installed, skipping 3D plotting")
    pass

from scipy.fft import fftn, ifftn
from numpy.fft import fftfreq, fftshift, ifftshift
from pyproj import Proj


MODELS = {
    # 1D PREM from IRIS EMC (warning might be too slow)
    "PREM": {
        "depth": 1E3 * np.array([
            -3.0, 0, 3.0, 15., 24.4, 71., 80., 171., 220., 271., 371., 400.
        ]),  # km
        "vp": 1E3 * np.array([
            1.45, 1.45, 5.80, 6.80, 8.11, 8.08, 8.08, 8.02, 7.99, 8.56, 8.66, 8.85
        ]),  # km/s
        "vs": 1E3 * np.array([
            1.0, 1.0, 3.20, 3.90, 4.49, 4.47, 4.47, 4.44, 4.42, 4.62, 4.68, 4.75
        ]),  # km/s
        "rho": 1E3 * np.array([
            1.02, 1.02, 2.6, 2.9, 3.38, 3.37, 3.37, 3.36, 3.36, 3.44, 3.47, 3.53
        ]),  # kg/m^3
        "qmu": np.array([
            600., 600., 600., 600., 600., 600., 600., 80., 80., 143., 143., 143.
        ]),
        # QK in PREM is inifinite, these values are inspired by Olsen et al 2018
        "qkappa": np.array([
            350., 350., 350., 350., 350., 350., 350., 200., 200., 200., 200., 200.
        ]),
        },
    # 1D PREM from SPECFEM3D
    "PREM_fast": {
        "depth": 1E3 * np.array([
            -3.0, 0, 3.0, 15., 24.4, 71., 80., 171., 220., 271., 371., 400.
        ]),  # km
        "vp": 1E3 * np.array([
            3.0, 3.0, 5.80, 6.80, 8.11, 8.08, 8.08, 8.02, 7.99, 8.56, 8.66, 8.85
        ]),  # km/s
        "vs": 1E3 * np.array([
            2.0, 2.0, 3.20, 3.90, 4.49, 4.47, 4.47, 4.44, 4.42, 4.62, 4.68, 4.75
        ]),  # km/s
        "rho": 1E3 * np.array([
            1.02, 1.02, 2.6, 2.9, 3.38, 3.37, 3.37, 3.36, 3.36, 3.44, 3.47, 3.53
        ]),  # kg/m^3
        "qmu": np.array([
            600., 600., 600., 600., 600., 600., 600., 80., 80., 143., 143., 143.
        ]),
        # QK in PREM is inifinite, these values are inspired by Olsen et al 2018
        "qkappa": np.array([
            350., 350., 350., 350., 350., 350., 350., 200., 200., 200., 200., 200.
        ]),
        },
    # HOMOGENEOUS HALFSPACE values from SPECFEM3D_Cartesian HH example
    "homogeneous_halfspace": {
        "depth":  np.array([0, 6371E3]),
        "vp":     np.array([2.8E3, 2.8E3]),
        "vs":     np.array([1.5E3, 1.5E3]),
        "rho":    np.array([2.3E3, 2.3E3]),
        "qmu":    np.array([300., 300.]),
        "qkappa": np.array([2444.4, 2444.4]),
        }
    }


def set_parameters(fid):
    """
    Load parameters from configuration file
    """
    try:
        with open(fid, "r") as f:
            parameters = json.load(f)
    except Exception as e:
        print("\n\nERROR READING JSON FILE") 
        traceback.print_exc()
        print("ERROR READING JSON FILE\n\n")
        sys.exit(-1)

    return parameters


def parse_args():
    """
    Allow plotting to be toggled from the command line
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("fid", type=str, help="file id")
    parser.add_argument("-A", "--all", action="store_true",
                        help="plot all figures")
    parser.add_argument("-m", "--map", action="store_true",
                        help="plot map view for a test parameter at a "
                             "random depth")
    parser.add_argument("-c", "--cross", action="store_true",
                        help="plot cross section through middle of model")
    parser.add_argument("-C", "--cube", action="store_true",
                        help="plot cube defining the whole 3D perturbed model, "
                             "require PyVista")
    parser.add_argument("-b", "--brick", action="store_true",
                        help="plot brick defining only the perturbation, "
                             "require PyVista")
    parser.add_argument("-p", "--profile", action="store_true",
                        help="plot 1D profile of the underlying velocity model")
    parser.add_argument("--cmap", default="RdYlBu", type=str, nargs="?",
                        help="colormap for matplotlib plots (map, cross)")
    parser.add_argument("--levels", default=31, type=int, nargs="?",
                        help="number of levels for colormap")
    parser.add_argument("--show", action="store_true",
                        help="show GUI for plots")
                                              
    return parser.parse_args()    


def lonlat_utm(lon, lat, utm_zone=None):                 
    """                                                                          
    convert latitude and longitude coordinates to UTM projection                 
    From Pyatoa                                                                  
                                                                                 
    :type lon: float
    :param lon: longitude value in WGS84 
    :type lat: float                                                  
    :param lat: latude value in WGS84 
    :type utm_zone: int                                                          
    :param utm_zone: UTM zone for conversion from WGS84                          
    :type inverse: bool                                                          
    :param inverse: if inverse == False, latlon => UTM, vice versa.              
    :rtype x: float                                                       
    :return x: x coordinate in UTM              
    :rtype y: float                                                       
    :return y: y coordinate in UTM         
    """                                                                          
    # Determine if the projection is north or south                              
    if utm_zone < 0:                                                             
        direction = "south"                                                      
    else:                                                                        
        direction = "north"                                                      
    # Proj doesn't accept negative zones                                         
    utm_zone = abs(utm_zone)                                                                       

    # Initiate a Proj object and convert the coordinates                         
    my_proj = Proj(proj="utm", zone=abs(utm_zone), south=True, ellps="WGS84",
                   datum="WGS84", units="m", no_defs=True)
    x, y = my_proj(lon, lat, inverse=False)            
                                                                                 
    return x, y         


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


def main(model_choice=None, tag=None, path=None,
         # Model definition parameters
         utm=None, xmin=None, xmax=None, ymin=None, ymax=None, 
         ZVALS=None, DX=None, DY=None, DZ=None, 
         # Perturbation flags
         perturbations=None, include_q=None, 
         # Perturbation control parameters
         seed=None, a=None, hv_ratio=None, mean_vel=None, std_vel=None, 
         nmin=None, nmax=None, zmin_pert=None, zmax_pert=None, perturb=None
         ):
    """
    Main function to generate a 3D tomography model with stochastic perturbation
    and export it to a format compatible with SPECFEM3D_Cartesian.

    :type fid: str
    :param fid: Path to the configuration file containing parameters.
    :type path: str
    :param path: Path to the directory where the model will be saved, 
        subdirectories will be created for figures and models.
    :type model_choice: str
    :param model_choice: The choice of the 1D model to use (e.g., "PREM").
    :type tag: str
    :param tag: A tag for output naming
    :type perturbations: bool
    :param perturbations: Whether to apply stochastic perturbations to the 
        model or just make a standard model
    :type include_q: bool
    :param include_q: Whether to include Q values in the model xyz files
    :type seed: int
    :param seed: Random seed for reproducibility of perturbations
    :type a: float
    :param a: Von Karman spectral filter parameter, controls the vertical scale 
        of perturbations. Use `hv_ratio` to select the horizontal scaling, 
        select `hv_ratio`=1 to get a_h == a_v
    :type hv_ratio: int
    :param hv_ratio: ratio of perturbation length in the horizontal direction
        to the vertical direction. 1 means a is the same for both
    :type mean_vel: float
    :param mean_vel: Mean velocity for the Gaussian perturbation distribution
    :type std_vel: float
    :param std_vel: Standard deviation of the velocity perturbation
    :type nmin: float
    :param nmin: Minimum value for the velocity perturbation normalization
    :type nmax: float
    :param nmax: Maximum value for the velocity perturbation normalization
    :type utm: int
    :param utm: Whether the coordinates are in UTM or not. If a value is given
        then the coordinates are assumed to be in lon/lat and will be
        converted to UTM in meters
    :type xmin: float
    :param xmin: Minimum x-coordinate for the model grid (in meters)
    :type xmax: float
    :param xmax: Maximum x-coordinate for the model grid (in meters)
    :type ymin: float
    :param ymin: Minimum y-coordinate for the model grid (in meters)
    :type ymax: float
    :param ymax: Maximum y-coordinate for the model grid (in meters)
    :type zmin_pert: float
    :param zmin_pert: Minimum z-coordinate for the perturbation brick [m]
    :type zmax_pert: float
    :param zmax_pert: Maximum z-coordinate for the perturbation brick [m]
    :type dx: float
    :param dx: Grid spacing in the x-direction (in meters)
    :type dy: float
    :param dy: Grid spacing in the y-direction (in meters)
    :type dz: float
    :param dz: Grid spacing in the z-direction (in meters)
    :type ZVALS: list of tuples
    :param ZVALS: List of tuples defining the z-coordinate ranges for each layer
        in the model. If None, defaults to zvals.
    :type DX: list
    :param DX: List of x-grid spacings for each layer. If None, defaults to dx.
    :type DY: list
    :param DY: List of y-grid spacings for each layer. If None, defaults to dy.
    :type DZ: list
    :param DZ: List of z-grid spacings for each layer. If None, defaults to dz.
    :type perturb: list
    :param perturb: List of parameters in `model` to apply perturbations to. 
        If None, defaults to all parameters in the model except 'depth'.
    """
    args = parse_args()
    parameters = set_parameters(args.fid)

    # Distribute the parameters to local variables because that's how the code
    # is written and I'm too lazy to change it to dictionary access
    model_choice = model_choice or parameters["model_choice"]
    tag = tag or os.path.basename(args.fid).split(".")[0]
    path = path or parameters["path"]

    # Model grid parameters    
    utm = utm or parameters["utm"]
    xmin = xmin or parameters["xmin"]
    xmax = xmax or parameters["xmax"]
    ymin = ymin or parameters["ymin"]
    ymax = ymax or parameters["ymax"]
    DX = DX or parameters["DX"]
    DY = DY or parameters["DY"]
    DZ = DZ or parameters["DZ"]
    ZVALS = ZVALS or parameters["ZVALS"]

    # Perturbation flags
    perturbations = perturbations or parameters["perturbations"]
    hv_ratio = hv_ratio or parameters["hv_ratio"]
    include_q = include_q or parameters["include_q"]

    # Perturbation parameters
    seed = seed or parameters["seed"]
    a = a or parameters["a"]
    mean_vel = mean_vel or parameters["mean_vel"]
    std_vel = std_vel or parameters["std_vel"]
    nmin = nmin or parameters["nmin"]
    nmax = nmax or parameters["nmax"]
    perturb = perturb or parameters["perturb"]
    zmin_pert = zmin_pert or parameters["zmin_pert"]
    zmax_pert = zmax_pert or parameters["zmax_pert"]

    # ==========================================================================
    # Set some additional parameters that are not in the config file
    if perturbations:
        tag = tag+"_wpert"
    else:
        tag = tag+"_nopert"

    # Note that these paths will be modified by `tag`
    mod_path = os.path.join(path, "created", tag)
    fig_path = os.path.join(mod_path, "figures")

    # Simple numbering scheme for tomography files
    fids = [f"tomography_model_{_:0>2}.xyz" for _ in range(1, len(ZVALS)+1)]
    fids = [os.path.join(mod_path, _) for _ in fids]

    # Misc parameters that we don't need to set
    indexing = "ij"
    order = "F"

    # Choose the Model
    model = MODELS[model_choice]

    # Make directories
    if not os.path.exists(mod_path):
        os.makedirs(mod_path)
    if any([args.cube, args.brick, args.profile, args.map, args.cross, args.all]):
        if not os.path.exists(fig_path):
            os.makedirs(fig_path)

    print(f"Tag is: {tag}")
    print(f"Figures will be saved to: {fig_path}")
    print(f"Models will be saved to: {mod_path}\n")

    # Drop attenuation if needed
    if not include_q:
        print("dropping 'Q' values from model")
        model = {key: val for key, val in model.items() 
                 if not key.startswith("q")}
        
    # Convert coordinates to UTM from lat/lon if needed
    if utm:
        xmin, ymin = lonlat_utm(xmin, ymin, utm_zone=utm)
        xmax, ymax = lonlat_utm(xmax, ymax, utm_zone=utm)

    print(f"vertical perturbation correlation length   a_v={a}m")
    print(f"horizontal perturbation correlation length a_h={a*hv_ratio}m\n")

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
        z_pert = np.arange(zmin_pert, zmax_pert, dz)  # z axis range, pert
        
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
            print("making meshgrid")
            [X, Y, Z_PERT] = np.meshgrid(x, y, z_pert, indexing=indexing)
            np.random.seed(seed)
            S = np.random.normal(mean_vel, std_vel, Z_PERT.shape)

            # 3D FFT into wave number domain. Shift for zero freq. at the center
            FS = fftshift(fftn(S))

            # Generate the frequency domain so that we can use it to index the 
            # perturbation. Shift it so that zero frequency is at the center, so 
            # that when we apply the perturbation it acts on the correct axes
            kx = fftshift(fftfreq(len(x), d=dx))
            ky = fftshift(fftfreq(len(y), d=dy))
            kz = fftshift(fftfreq(len(z_pert), d=dz))
            [KX, KY, KZ] = np.meshgrid(kx, ky, kz, indexing=indexing)

            # Scale the wavenumber vectors so that the perturbations come out
            # the same size, otherwise they will scale with the length of the 
            # domain which is not what we want
            scale_x = xmax - xmin
            scale_y = ymax - ymin
            scale_z = zmax_pert - zmin_pert
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

            # Plot volumetric cube with PyVistaw

            if args.brick or args.all:
                data = pv.wrap(S_pert)
                data.plot(volume=True, cmap=args.cmap)
                if args.show:
                    plt.show()
                plt.close("all")

        # GENERATE VELOCITY MODEL
        # Now we generate the velocity model cube. First interpolate to our `dz`
        [X, Y, Z] = np.meshgrid(x, y, z, indexing=indexing)
        print(f"layer '{fids[l]}' has {np.prod(X.shape)} points")

        # Interpolate the 1D model to the required DZ
        model = interp_1D_model(model=model, dz=dz, zmin=zmin, zmax=zmax)

        if args.profile or args.all:
            # Just plot one example parameter 
            f, ax = plt.subplots(figsize=(6,8))
            plt.plot(model["vs"], model["depth"], "ko-")
            plt.gca().invert_yaxis()
            plt.xlabel("Vs [m/s]")
            plt.ylabel("Depth [m]")
            plt.title(f"Interpolated 1D Model {model_choice}\n{fids[l]}")
            plt.grid()
            plt.savefig(os.path.join(fig_path, f"1d_profile_{zmin}-{zmax}.png"))
            if args.show:
                plt.show()
            plt.close("all")

        # Build model for each of the model parameters. Add perturbations 
        model_dict = {}
        map_plotted, cross_plotted = False, False
        for parameter in model.keys():
            # Ignore depth values, 
            if parameter == "depth":
                continue
            
            print(f"making model for par: {parameter}")
            arr = model[parameter]

            # Generate the 3D cube shape
            cube = np.ones((len(x), len(y), len(arr)))

            # Fill in the cube by multiplying the correct depth with their 
            # assigned value. This is pretty brute force but it works
            for i, val in enumerate(arr):
                # Apply perturbation ontop of cube if we are in the correct layer
                if (parameter in perturb) and \
                    pert_idx_start is not None and \
                        pert_idx_end is not None and \
                        (pert_idx_start <= i < pert_idx_end): 
                    j = i - pert_idx_start  # set correct index for perturbation
                    cube[:,:,i] *= (val + val * S_pert[:,:,j])
                # If not, then we just assign the correct value
                else:
                    cube[:,:,i] *= val  
            
            # Store the unraveled 1D array for writing
            # Flip the values because we will flip all the spatial coordinates 
            # so that deepest points are listed first in the XYZ model
            model_dict[parameter] = np.flip(cube.ravel(order=order))

            # Plot test depth value and parameter
            if (args.map or args.all) and not map_plotted:
                # Plot map view
                cmap = plt.get_cmap(args.cmap, args.levels)
                im = plt.imshow(
                        cube[:,:,1], cmap=cmap, 
                        extent=np.array([0, xmax-xmin, 0, ymax-ymin]) * 1E-3
                        )
                plt.colorbar(im, shrink=0.8, pad=0.025, label=parameter)
                plt.grid()
                plt.xlabel("X [km]")
                plt.ylabel("Y [km]")
                ax = plt.gca()
                # ax.ticklabel_format(style="sci", axis="both", scilimits=(0,0))
                ax.ticklabel_format(style="plain", axis="both")
                ax.set_aspect("equal")
                plt.title(f"{parameter} at {arr[1]:.2f}m")
                plt.tight_layout()
                plt.savefig(os.path.join(fig_path, 
                                        f"2d_plot_{parameter}_{zmin}-{zmax}.png")
                                        )
                if args.show:
                    plt.show()
                plt.close("all")
                map_plotted = True
                sys.exit()
            
            if (args.cross or args.all) and not cross_plotted:
                # Plot cross section
                breakpoint()
                cmap = plt.get_cmap(args.cmap, args.levels)
                im = plt.imshow(
                        cube[:,:,1], cmap=cmap, 
                        extent=np.array([0, xmax-xmin, 0, ymax-ymin]) * 1E-3
                        )
                plt.colorbar(im, shrink=0.8, pad=0.025, label=parameter)
                plt.grid()
                plt.xlabel("X [km]")
                plt.ylabel("Y [km]")
                ax = plt.gca()
                # ax.ticklabel_format(style="sci", axis="both", scilimits=(0,0))
                ax.ticklabel_format(style="plain", axis="both")
                ax.set_aspect("equal")
                plt.title(f"{parameter} at {arr[1]:.2f}m")
                plt.tight_layout()
                plt.savefig(os.path.join(fig_path, 
                                        f"2d_plot_{parameter}_{zmin}-{zmax}.png")
                                        )
                if args.show:
                    plt.show()
                plt.close("all")
                cross_plotted = True

        # Plot volumetric cube with PyVista
        if args.cube or args.all:
            data = pv.wrap(cube[:,:,::-1])  # Flip the cube so Z positive down
            data.plot(volume=True, cmap=args.cmap)
            if args.show:
                plt.show()
            plt.close("all")

        # EXPORT MODEL
        # Ravel turns the gridded data into 1D arrays that can be written
        # We flip the arrays so that when written, the largest negative value 
        # starts Z axis is multiplied by -1 so that positive Z is defined as up
        x_list = np.flip(X.ravel(order=order))
        y_list = np.flip(Y.ravel(order=order))
        z_list = np.flip(-1 * Z.ravel(order=order))

        with open(fids[l], "w") as f:
            # Header - min and max range values
            f.write(
                f"{x_list.min():.1f} {y_list.min():.1f} {z_list.min():.1f} "
                f"{x_list.max():.1f} {y_list.max():.1f} {z_list.max():.1f}\n"
                )

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
            

if __name__ == "__main__":
    main()
