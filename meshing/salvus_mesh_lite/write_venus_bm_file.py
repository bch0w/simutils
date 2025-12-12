"""
Create input .bm files used for Salvus mesher for AxiSEM3D. Make custom tweaks
to the underlying velocity models in the code so that they are preserved
and not have to be done manually each time

Notes
- at r=6051.9km is the START of the atmosphere (vs=0 here)
"""
import os
import numpy as np
from glob import glob


# Venus specific
vradius_idx = 501  # index for first atmo layer at 6051.9 km
sponge_amt = 0.01  # for atmo, add extra atmo as absorbing boundary

# Paths
fpath = "/Users/chow/Work/research/venusseis/simulations/dumoulin2017_1d_models"
opath = "bm_files"
model_paths = glob(os.path.join(fpath, "*.out"))


# Loop on models
for path_ in model_paths:
    model = os.path.basename(path_).split(".")[0]  # e.g., V1-Tc.out
    print(model)

    # Write out values from Dumoulin model
    radius, temp, _, rho, vp, vs = np.loadtxt(path_, skiprows=1).T

    # Different configurations
    for key in ["w_atmo", "no_atmo"]:
        # Keep atmosphere, add sponge layer above
        if key == "w_atmo":  
            rmax = radius.max()
            buffer = rmax * sponge_amt

            print(f"atmo sponge adding {buffer}km buffer")
            # Extrapolate last atmo value out to the buffer value
            radius = np.append(radius, [rmax + buffer])
            rho = np.append(rho, [rho[-1]])
            vp = np.append(vp, [vp[-1]])
            vs = np.append(vs, [vs[-1]])

        # Strip the atmosphere from the mesh/model
        elif key == "no_atmo":
            print("stripping atmo layer from model")
            radius = radius[:vradius_idx]
            rho = rho[:vradius_idx]
            vp = vp[:vradius_idx]
            vs = vs[:vradius_idx]

        fid_out = f"{model}_{key}.bm"
        path_out = os.path.join(opath, fid_out)

        # Write file header
        with open(path_out, "w") as f:
            f.write(f"NAME {fid_out.split('.')[0]}\n")
            f.write(f"ANELASTIC T\n")
            f.write(f"ANISOTROPIC F\n")
            f.write(f"UNITS m\n")
            f.write(f"COLUMNS radius rho vp vs qka qmu\n")

            # Write file data points, starting from the core and moving up
            for i, r in enumerate(radius):
                f.write("    ")
                if rho[i] > 1:  # density in atmosphere can go below 1
                    f.write(f"{r*1E3:7.0f}.  {rho[i]:8.2f}  {vp[i]*1E3:8.2f}  "
                            f"{vs[i]*1E3:8.2f}")
                else:
                    f.write(f"{r*1E3:7.0f}.  {rho[i]:8.2E}  {vp[i]*1E3:8.2f}  "
                            f"{vs[i]*1E3:8.2f}")
                # Append attenuation based on radius
                if r > 6051.9 and vs[i] == 0:  # atmosphere 
                    f.write(f"  {0:8.1f}  {0:8.1f}")
                elif vs[i] == 0:
                    f.write(f"  {10000:8.1f}  {0:8.1f}")  # core
                else:
                    f.write(f"  {10000:8.1f}  {300:8.1f}")  # rest
                f.write("\n")
                
