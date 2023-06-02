"""
Utility function to write on proj_grid_fd.txt for xproject SPECFEM script
given some mesh dimensions and desired FD grid size

# Original Values
origin_x = 322697.19
origin_y = 3663978.75
origin_z = -30000.

end_x = 416983.19
end_y = 3773609.43
end_z = 300.

origin: starting location in unit meters 
end: ending location in unit meters
n?: number of samples in a given direction
h?: sampling rate in a given direction
choice: 'n' for hard-coding n?, 'h' for hard-coding h?
"""
test_case = dict(
    choice = "n",
    origin_x = 322500.00,
    origin_y = 3663500.00,
    origin_z = -30000.,
    end_x = 416500.00,
    end_y = 3773500.00,
    end_z = 500.,
    nx = 100.,
    ny = 100.,
    nz = 500.,
    hx = 1000.,
    hy = 1000.,
    hz = 50.,
    )

nzatom = dict(
    origin_x = 170000.,
    origin_y = 5286000.,
    origin_z = -400000.,
    end_x = 635000.,
    end_y = 5905000.,
    end_z = -44000.,
    nx = 100.,
    ny = 100.,
    nz = 500.,
    hx = 4000.,
    hy = 4000.,
    hz = 4000.,
    choice = "h",
    )


# Assign a shorter variable name for dict
d = nzatom

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




