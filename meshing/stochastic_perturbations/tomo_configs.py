# ==========================================================================
#                               PARAMETERS
# ==========================================================================
choice = PREM

# GRID SPACING LISTS [m]
# d? lenght should be 1-len(zvals)
DX = [0.5E3, 25E3]
DY = [0.5E3, 25E3]
DZ = [0.5E3, 10E3]
ZVALS = [0, 10E3, 200E3]  # positive down, we will flip this later

# DEFINE FULL DOMAIN [m]
xmin = 245.750E3
xmax = 890.650E3
ymin = 4487.550E3
ymax = 5050.670E3

# PERTURBATION
a = 5E3  # m
nmin = -0.1
nmax = 0.1
zmin_pert = 0.  # depth extent of the perturbation
zmax_pert = 10.E3
seed = 123
mean_vel = 1  # km/s
std_vel = 0.1

# PLOTTING
plot_cube = False    # model
plot_brick = False  # perturbation
plot_2d = False
cmap = "viridis"

# EXPORT
# fid length should match 1-len(ZVALS)
fids = [f"tomography_model_{_:0>2}.xyz" for _ in range(1, 5)]

# MISC
indexing = "ij"
order = "F"
