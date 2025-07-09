# NK Test model 6/13/25
model_choice = "PREM"
tag = "nk_mdj_lowres"
perturbations = False
include_q = True

# GRID SPACING LISTS [m] (d? lenght should be 1-len(zvals))
DX = [1E3, 1E3, 5.0E3] 
DY = [1E3, 1E3, 5.0E3] 
DZ = [0.5E3, 1E3, 5.0E3] 
ZVALS = [[-2.5E3, 5E3], [5E3, 11E3], [11E3, 100E3]]  # positive down, we will flip this later

# DEFINE FULL DOMAIN [m]
xmin = 483230.
xmax = 563272.
ymin = 4560978.
ymax = 4961045.

# PERTURBATION PARAMETERS
a = 1E3  # m
nmin = -0.1
nmax = 0.1
perturb = ["vp", "vs", "rho"]  # parameters to apply perturbations to
zmin_pert = -2.5E3  # depth extent of the perturbation
zmax_pert = 10.E3
seed = 123
mean_vel = 1  # km/s
std_vel = 0.1
