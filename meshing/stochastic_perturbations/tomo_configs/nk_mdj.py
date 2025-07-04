# NK Test model 6/13/25
model_choice = "PREM"
tag = "nk_mdj"
perturbations = True
include_q = True

# GRID SPACING LISTS [m] (d? lenght should be 1-len(zvals))
DX = [0.125E3, 0.5E3, 1.0E3] 
DY = [0.125E3, 0.5E3, 1.0E3] 
DZ = [0.125E3, 0.5E3, 1.0E3] 
ZVALS = [[0, 5E3], [5E3, 10E3], [10E3, 100E3]]  # positive down, we will flip this later

# DEFINE FULL DOMAIN [m]
xmin = 491650.
xmax = 563385.
ymin = 4594265.
ymax = 4949940.

# PERTURBATION PARAMETERS
a = 1E3  # m
nmin = -0.1
nmax = 0.1
perturb = ["vp", "vs", "rho"]  # parameters to apply perturbations to
zmin_pert = -1E3  # depth extent of the perturbation
zmax_pert = 10.E3
seed = 123
mean_vel = 1  # km/s
std_vel = 0.1
