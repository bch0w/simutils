# NK Test model 6/13/25
model_choice = "PREM"
perturbations = False

# GRID SPACING LISTS [m] (d? lenght should be 1-len(zvals))
# We don't need the double layer but it helps to match the existing mesh
DX = [25E3, 25E3]
DY = [25E3, 25E3]
DZ = [2E3, 10E3]
ZVALS = [0, 10E3, 200E3]  # positive down, we will flip this later

# DEFINE FULL DOMAIN [m]
xmin = 245.750E3
xmax = 890.650E3
ymin = 4487.550E3
ymax = 5050.670E3

# PERTURBATION PARAMETERS
a = 5E3  # m
nmin = -0.1
nmax = 0.1
zmin_pert = 0.  # depth extent of the perturbation
zmax_pert = 10.E3
seed = 123
mean_vel = 1  # km/s
std_vel = 0.1
