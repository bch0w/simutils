"""
Generate a 1D velocity model (e.g., PREM) as a 3D external tomography model
exportable to SPECFEM3D_Cartesian. 
"""
import numpy as np
import matplotlib.pyplot as plt


def PREM():
    """
    Approximately defines the 1D PREM model. Values were listed in the 
    Bormann (2002) NMSOP Table 1.1. The values are not exact, but are close 
    enough
    """
    prem = {
        "depth": np.array([
            0.0, 3.0, 15., 24.4, 71., 80., 171., 220., 271., 371., 400.
        ]),
        "vp": np.array([
            1.45, 5.80, 6.80, 8.11, 8.08, 8.08, 8.02, 7.99, 8.56, 8.66, 8.85
        ]),
        "vs": np.array([
            1.0, 3.20, 3.90, 4.49, 4.47, 4.47, 4.44, 4.42, 4.62, 4.68, 4.75
        ]),
        "rho": np.array([
            1.02, 2.6, 2.9, 3.38, 3.37, 3.37, 3.36, 3.36, 3.44, 3.47, 3.53
        ]),
        "qmu": np.array([
            0., 600., 600., 600., 600., 600., 80., 80., 143., 143., 143.
        ]),
        "qkappa": np.array([57323.] * 11)
        }
    
    prem["depth"] *= 1E3  # Convert to meters with positive up

    return prem

def AK135F(depth_cutoff_km=80, plot=False):
    """
    Read and parse AK135F model as gathered directly from IRIS EMC
    """
    ak135f = np.loadtxt("/Users/prof/Data/models/AK135F/AK135F_AVG.csv",
                        delimiter=",", dtype=float)
    depth_km, density, vp_kms, vs_kms, qk, qm = ak135f.T

    # Determine where the cutoff depth index is
    idxs = np.where(depth_km < depth_cutoff_km)

    depth = depth_km[idxs] * 1E3  # m, -Z up
    density = density[idxs] * 1E3
    vp = vp_kms[idxs] * 1E3  # m
    vs = vs_kms[idxs] * 1E3  # m
    qp = qk[idxs]  # this is the wrong name but .xyz file takes qk and qm
    qs = qm[idxs]

    # Turn the water and mud layer into solid layer. This matches what is done
    # for AK135-f Syngine model, where 0-10km is a single-valued layer
    # Solid layer starts at index 4, above that is water and mud
    density[:4] = density[4]
    vp[:4] = vp[4]
    vs[:4] = vs[4]
    qp[:4] = qp[4]
    qs[:4] = qs[4]

    if plot:
        plt.plot(vp, depth, "bo-", label="Vp")
        plt.plot(vs, depth, "ro-", label="Vs")
        plt.show()
        a=1/0
    
    # Remove values below depths we don't care about
    ak135f = {"depth": depth, "vp": vp, "vs": vs, 
              "rho": density, "qp": qp, "qs": qs}

    return ak135f

def extract_1d_profile(xyz_fid, xselect, yselect):
    """
    Extract a 1D velocity profile from a 3D velocity model. Allows generating
    arbitrary 1D profiles depending on the model. Must be loaded from a
    SPECFEM .xyz external tomography model

    NK6 in UTM = 502485.224644004, 4575658.10875230
    """
    x, y, z, vp, vs, rho, *q = np.loadtxt(xyz_fid, skiprows=4).T
    if q:
        qp, qs = q
    else:
        qp, qs = None, None

    # Find nearest model x and y values to the selected points
    xclosest = x[(np.abs(x - xselect)).argmin()]
    yclosest = y[(np.abs(y - yselect)).argmin()]
    
    # Select indices that only match the chosen x/y position
    idxs = np.where((x==xclosest) & (y==yclosest))[0]

    # If no attenuation in model then we input a homogeneous Q model
    if q:
        qp = q[0][idxs]
        qs = q[1][idxs]
    else:
        qp_0 = 350
        qs_0 = 350
        qp = np.array([qp_0] * len(idxs))
        qs = np.array([qs_0] * len(idxs))
    
    # Interpolate onto given depth values
    model = {"depth": z[idxs], "vp": vp[idxs], "vs": vs[idxs],
             "rho": rho[idxs], "qp": qp, "qs": qs
             }
    dx = np.unique(x)[1] - np.unique(x)[0]
    dy = np.unique(y)[1] - np.unique(y)[0]
    dz = np.unique(z)[1] - np.unique(z)[0]

    return model, dx, dy, dz
    

def interp_1D_model(model, dz, plot=True):
    """
    1D interpolation of the 1D model to get regular step sizes
    """
    # Define the new depth values to intepolate against
    z = model["depth"]  # Assume -Z up
    zs = np.arange(z.min(), z.max() + dz, dz)
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

    # Re-define Z-axis, +Z up
    model_out["depth"] *= -1  # m

    # Flip every array so that it conforms with the structure required by SPECFEM
    model_out = {key: val[::-1] for key, val in model_out.items()}

    if plot:
        plt.plot(model_out["vp"], model_out["depth"], "bo-", label="Vp")
        plt.plot(model_out["vs"], model_out["depth"], "ro-", label="Vs")
        plt.show()

    return model_out


def make_model(model, X, Y, fid="tomography_model.xyz"):
    """
    Generates a 1D model for SPECFEM3D_Cartesian. The model is defined by 
    the depth, vp, vs, rho, qmu, and qkappa values. The model is then 
    exported to a file in the SPECFEM3D_Cartesian format.
    """
    # Flip the Z axis because positive is up which means arange flipped it 
    # previously
    Z = model["depth"]

    with open(fid, "w") as f:
        # Header - min and max range values
        f.write(f"{X.min():.1f} {Y.min():.1f} {Z.min():.1f} "
                f"{X.max():.1f} {Y.max():.1f} {Z.max():.1f}\n")
        # Header - spacing values
        f.write(f"{X[1]-X[0]:.1f} {Y[1]-Y[0]:.1f} {Z[1]-Z[0]:.1f}\n")
        # Header - number of grid points
        f.write(f"{int(len(X)):d} {int(len(Y)):d} {int(len(Z)):d}\n")
        # Header - parameter min max values
        f.write(f"{model['vp'].min():.1f} {model['vp'].max():.1f} "
                f"{model['vs'].min():.1f} {model['vs'].max():.1f} "
                f"{model['rho'].min():.1f} {model['rho'].max():.1f}\n")
        
        for i, z in enumerate(Z):
            for y in Y:
                for x in X:
                    f.write(f"{x:.1f} {y:.1f} {z:.1f} "
                            f"{model["vp"][i]:.1f} {model["vs"][i]:.1f} "
                            f"{model["rho"][i]:.1f} {model["qp"][i]:.1f} "
                            f"{model["qs"][i]:.1f}\n"
                            )
     

if __name__ == "__main__":
    # User input parameter - all units in meters
    # ABC MODEL
    xmin = 425_000.
    xmax = 656_000.
    ymin = 4_490_000.
    ymax = 4_680_000.
    
    # CORRIDOR MODEL
    # xmin = 492_485.2
    # xmax = 557_079.1
    # ymin = 4_565_658.1
    # ymax = 4_950_643.3

    extract = False  # if True, pull 1D model from a 3D tomography model

    print(f"X={(xmax-xmin)*1E-3:.2f}km")
    print(f"Y={(ymax-ymin)*1E-3:.2f}km")

    # Make model here
    if not extract:
        dx = 0.25E3  # m
        dy = 0.25E3 
        dz = .1E3  # small enough so we can get the sharp jumps across layers
        model = AK135F()
    else:
        # NK6 location for SimBlast paper
        xselect = 502485
        yselect = 4575658
        model, dx, dy, dz = extract_1d_profile(xyz_fid="tomography_model.xyz",
                                               xselect=xselect, yselect=yselect)
        
    # Creates the horizontal grid of points
    X = np.arange(xmin - dx, xmax + dx, dx)
    Y = np.arange(ymin - dx, ymax + dy, dy)

    interp_model = interp_1D_model(model, dz=dz)
    make_model(interp_model, X=X, Y=Y)
