"""
Generate a 1D velocity model (e.g., PREM) as a 3D external tomography model
exportable to SPECFEM3D_Cartesian. 
"""
import numpy as np


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
    
    prem["depth"] *= 1E3  # Convert to meters

    return prem

def interp_1D_model(model, dz):
    """
    1D interpolation of the 1D model between major depth values to get gradients 
    in between rather than just step functions.
    """
    # Define the new depth values to intepolate against
    z = model["depth"]
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
    
    return model_out

def make_model(model, X, Y, fid="tomography_model.xyz"):
    """
    Generates a 1D model for SPECFEM3D_Cartesian. The model is defined by 
    the depth, vp, vs, rho, qmu, and qkappa values. The model is then 
    exported to a file in the SPECFEM3D_Cartesian format.
    """
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
                            f"{model["rho"][i]:.1f} {model["qmu"][i]:.1f} "
                            f"{model["qkappa"][i]:.1f}\n"
                            )
     

if __name__ == "__main__":
    # User input parameter - all units in meters
    dx = 50E3 
    dy = 50E3 
    dz = 1E3  

    xmin = 245.750E3
    xmax = 890.650E3
    ymin = 4487.550E3 
    ymax = 5050.670E3

    # Creates the horizontal grid of points
    X = np.arange(xmin - dx, xmax + dx, dx)
    Y = np.arange(ymin - dx, ymax + dy, dy)

    # Make model here
    model = PREM()
    interp_model = interp_1D_model(model, dz=dz)
    make_model(interp_model, X=X, Y=Y)
