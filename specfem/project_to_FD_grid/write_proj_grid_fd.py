"""
Utility function to write on 'proj_grid_fd.txt' for xproject SPECFEM script
given a User-requested finite-difference (FD) grid projection. 

.. parameters::

    origin_x: starting X location of FD grid in unit meters 
    origin_y: starting Y location of FD grid in unit meters 
    origin_z: starting Z location of FD grid in unit meters. This should be the 
        BOTTOM of your mesh or desired FD grid
    end_x/end_y/end_z: ending location of the FD grid in unit meters. `end_z`
        shoudl be the TOP fo your mesh or desired FD grid
    nx/ny/nz: number of samples in a given direction (x, y, z)
    hx/hy/hz: sampling rate in a given direction (x, y, z) in unit meters
    choice: 'n' - define the FD grid based on nx/ny/nz, ignore hx/hy/hz
            'h' - define the FD grid based on hx/hy/hz, ignore nx/ny/nz

An example output file will look like this:

ox oy oz
hx hy hz
nx ny nz

Where ox/y/z are the origin of the FD grid, hx/y/z defines the sampling rate
and nx/y/z defines the number of points in each direction.
"""
import sys

# >>> Template, copy this and fill in your own values based on Mesh_Par_file
fd_grid = dict(
    choice = "n",
    # Origin (0, 0, 0) of your Mesh
    origin_x =   
    origin_y = 
    origin_z = 
    # Top-right corner of your mesh
    end_x = 
    end_y = 
    end_z = 
    # Choose either to define your grid by nx/ny/nz
    nx = 
    ny = 
    nz = 
    # Or to define your grid by hx/hy/hz
    hx = 
    hy = 
    hz = 
    )

# The following is an example from the Chow et al. (2022a) model published on 
# IRIS EMC
nzatom_crust = dict(
    origin_x = 171000.,
    origin_y = 5286000.,
    origin_z = -50000.,
    end_x = 631000.,
    end_y = 5902000.,
    end_z = -7000.,
    nx = 100.,
    ny = 100.,
    nz = 500.,
    hx = 1000.,
    hy = 1000.,
    hz = 500.,
    choice = "h",
    )

# File writing starts here
d = fd_grid

# Set by number of points here
if d["choice"] == "n":
    nx = d["nx"]
    ny = d["ny"]
    nz = d["nz"]
    # Calculates h values based on desired start and end
    hx = (d["end_x"] - d["origin_x"]) // nx
    hy = (d["end_y"] - d["origin_y"]) // ny
    hz = (d["end_z"] - d["origin_z"]) // nz
elif d["choice"] == "h":
    hx = d["hx"]
    hy = d["hy"]
    hz = d["hz"]
    # Indexing starts at 1 (?)
    nx = 1 + (d["end_x"] - d["origin_x"]) // hx
    ny = 1 + (d["end_y"] - d["origin_y"]) // hy
    nz = 1 + (d["end_z"] - d["origin_z"]) // hz

print(f"hx={hx}, hy={hy}, hz={hz}")
print(f"nx={nx}, ny={ny}, nz={nz}")

with open("fd_proj_grid.txt", "w") as f:
    f.write(f"{d['origin_x']:.2f} {d['origin_y']:.2f} {d['origin_z']:.2f}\n")
    f.write(f"{int(hx)} {int(hy)} {int(hz)}\n")
    f.write(f"{int(nx)} {int(ny)} {int(nz)}\n")




