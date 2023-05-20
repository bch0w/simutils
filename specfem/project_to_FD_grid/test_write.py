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
"""
# Define mesh start and end here (output_mesher.txt)
origin_x = 322500.00
origin_y = 3663500.00
origin_z = -30000.

end_x = 416500.00
end_y = 3773500.00
end_z = 500.

# Set by number of points here
if True:
    # Define FD grid here (n?=numpoints in FD grid for ? dimension)
    nx = 10.
    ny = 10. 
    nz = 5.

    # nx += 1 
    # ny += 1
    # nz += 1

    # Calculates h values based on desired start and end
    hx = (end_x - origin_x) // nx
    hy = (end_y - origin_y) // ny
    hz = (end_z - origin_z) // nz
else:
    # Define grid spacing here (units of meters)
    hx = 1000
    hy = 1000
    hz = 50

    # Indexing starts at 1 (?)
    nx = 1 + (end_x - origin_x) // hx
    ny = 1 + (end_y - origin_y) // hy
    nz = 1 + (end_z - origin_z) // hz


print(f"hx={hx}, hy={hy}, hz={hz}")
print(f"nx={nx}, ny={ny}, nz={nz}")

with open("fd_proj_grid.txt", "w") as f:
    f.write(f"{origin_x:.2f} {origin_y:.2f} {origin_z:.2f}\n")
    f.write(f"{int(hx)} {int(hy)} {int(hz)}\n")
    f.write(f"{int(nx)} {int(ny)} {int(nz)}\n")




