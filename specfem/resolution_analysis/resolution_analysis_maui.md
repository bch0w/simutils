# Resolution Analysis on Maui

Bryant Chow 
5/15/23

The following document describes how to run the resolution analysis performed in Chow et al. (2022a,b) 
using SPECFEM3D_Cartesian on Maui. These analyses include a Zeroth Moment test and Point Spread
Function test, outlined in various publications from Andreas Fichtner et al. (see background).

## Background

The general idea here is that we are approximating the Hessian (derivative of the gradient) to get
some information about resolution of given features.

Mathematically, we are calculating the "action Hessian", which is the finite difference between two
gradients. That is, the Hessian `H(m)` can be calculated as `H(m) = g(m + dm) - g(m)`, where `g(m)`
is the gradient of your final model, and `g(m + dm)` is the gradient calculated using a perturbed
version of your final model. The choice of perturbation `dm` determines what type of test we 
are running.

Practically, this involves N additional simulations per test, where N is the total number of events. 
You should already have `g(m)` from your final simulation run. That means we just need to 
calculate `g(m+dm)` by perturbing our final model `m` by some quantity `dm`, and calculating the gradient 
`g(m+dm)` using an adjoint simulation. We then take the difference of the two gradients to retrieve `H(m)`.

### Zeroth Moment Test (Volumentric Point Spread Function)

"The zeroth moment, M(0)(x), is equal to the integral over the PSF for a point perturbation at position x. 
It follows that M(0)(x) is small when the misfit χ is nearly unaffected by a point perturbation at x. In 
contrast, M(0)(x) is large when the PSF has a high amplitude, a large spatial extent, or both. In this sense, 
the spatial variability of M(0)(x) reflects the relative weight of neighbouring PSFs. Information on the 
resolution length is not contained in the zeroth moment."

In the Zeroth Moment test our perturbation `dm` is a homogeneous offset to all values in our final model `m`.
Based on Fichtner et al. (2015), perturbation within a few percent should provide a meaningful choice. In
Chow et al. (2022a), we used 5\% perturbation of the given average velocity.

### Point Spread Function Test

"The Hessian represents our blurred perception of a point-localized perturbation at position 
y in a linearized tomographic inversion. The effect of the off-diagonal elements Hij|i≠j is 
to introduce unwanted updates of model parameters mj|i≠j that have initially not been perturbed." 
-- Fichter et al. (2011)

In a point spread function test, we perturb our final model with a point/spike perturbation and attempt
to recover it. The size and shape of the given perturbation depend on the numerical resolution of
your inversion setup, and also the feature of interest that one is probing. The theoretical lower limit 
of the size of the spike is around half a wavelength. In Chow et al. (2022a,b) we took a more conservative 
approach and generated perturbation of similar size as the features we were interested in.


#### Size of PSF 
"Most of the requirements imposed by numerical methods can be met by choosing the width of the random 
perturbations to be **around half a wavelength**... Wider perturbations, for instance, in response to numerical 
restrictions, will generally increase the estimated resolution lengths, thus making them more conservative." 

#### Amplitude of Perturbation
"...travel times measured by cross correlation are linearly related to seismic velocity 
perturbations of up to 10% [Mercerat and Nolet, 2013], meaning that **random velocity perturbations in the percent range 
would be a meaningful choice."** 

**References**:
1. [Fichter et al. (2011)](https://academic.oup.com/gji/article/187/3/1604/616815)
2. [Fichter et al. (2015)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015JB012106)
3. [Chow et al. (2022a)](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2021JB022865)
4. [Chow et al. (2022b)](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2021JB022866)

## Perturbing Velocity Model

To perform these tests, you will need to perturb your velocity model (`m + dm`) and run simulations through this new 
model. This is straighforward for the Zeroth moment, more involved for PSF.

### Zeroth Moment

From the working directory of your final velocity model, you can run `seisflows debug` and copy the following
Python code block, which will load in your model, perturb it, and write it back out in SPECFEM format.

1. Go to the SeisFlows working directory where your final model is
2. Run `seisflows debug` to open up interactive mode
3. Copy-paste the following code block into the terminal
4. The resulting perturbed model should be saved in *PATH.OUTPUT*

```python
"""
To be run during SeisFlows debug mode to generate a perturbed velocity model
`m + dm` for Zeroth Moment resolution analysis. Reads in the existing model
from the optimization library, perturbs it and saves the resulting 
perturbed model in SPECFEM-ready format.
"""
import os

# Choose perturbation amount here
perturbation = 0.05  # percentage

# Update model parameters by perturbing
m = optimize.load("m_new")
n = int(len(m) / 2)  # assuming two model parameters (Vp, Vs)

# Apply perturbation to actual model values
dm = m * perturbation  # perturbed model
mdm = m + dm

# Quick comparison of model values
print("original model")
print(f"min_0={m[:n].min()}, max_0={m[:n].max()}")
print(f"min_1={m[n:].min()}, max_1={m[n:].max()}")

print("perturbed model")
print(f"min_0={mdm[:n].min()}, max_0={mdm[:n].max()}")
print(f"min_1={mdm[n:].min()}, max_1={mdm[n:].max()}")

# Save perturbations to requisite directories. Do one at a 
# time Vs then Vp
for i, dir_ in enumerate([f"output/dvs_{int(perturbation*1E2):d}pct", 
                          f"output/dvp_{int(perturbation*1E2):d}pct"]):
    if not os.path.exists(dir_):
        os.mkdir(dir_)
    os.chdir(dir_)
    mdm_ = mdm.copy()
    if i == 0:
        mdm_[:n] *= 0  # zero out all Vp values
    elif i == 1:
        mdm_[n:] *= 0  # zero out all Vs values
    solver.save(save_dict=solver.split(mdm_), path="./")
    os.chdir("..")
```

### Point Spread Function

To generate a point spread function, I wrote the [following Python script](https://github.com/bch0w/simutils/blob/master/meshing/point_local.py). 
The idea behind the script is that it takes an external tomography model (.xyz files), and zeros out all values within the domain, except for 
a user-defined 3D Gaussian ellipsoid at a given location.

Once these .xyz files are created, you must manually generate databases using SPECFEM to get them into SPECFEM-ready format. Once these
PSF database files are ready, you can read them in and perturb your final velocity model by running the following code block in 
SeisFlows debug mode.


```python
"""
To be copy-pasted into the SeisFlows debug mode to generate the correct 
perturbation files for use in the point spread tests
"""
import os
assert(os.path.basename(os.getcwd()) == "original"), "wrong directory"
vp_offset = 3000  # set by point_local.py because otherwise SPECFEM 
vs_offset = 1500  # complains when most velocity values are zero (unphysical)
perturbation = 0.05

# Load perturbation and remove offsets
dm = solver.merge(solver.load('./'))
n = int(len(dm)/2)
dm[:n] -= vp_offset
dm[n:] -= vs_offset
print(f"min_0={dm[:n].min()}, max_0={dm[:n].max()}")
print(f"min_1={dm[n:].min()}, max_1={dm[n:].max()}")

# Save as template for future use
os.chdir("..")
if not os.path.exists("template"):
    os.mkdir("template")
os.chdir("template")
solver.save(save_dict=solver.split(dm), path="./")
os.chdir("..")

# Apply perturbation to actual model values
m = solver.merge(solver.load(PATH.MODEL))
dm *= perturbation
dm *= m
print(f"min_0={dm[:n].min()}, max_0={dm[:n].max()}")
print(f"min_1={dm[n:].min()}, max_1={dm[n:].max()}")

# Save perturbations to requisite directories
for i, dir_ in enumerate([f"dvs_{int(perturbation*1E2):d}pct", 
                          f"dvp_{int(perturbation*1E2):d}pct"]):
    if not os.path.exists(dir_):
        os.mkdir(dir_)
    os.chdir(dir_)
    dm_ = dm.copy()
    if i == 0:
        dm_[:n] *= 0  # dvs
    elif i == 1:
        dm_[n:] *= 0  # dvp
    print(f"{dir_}: min={dm_.min()}, max={dm_.max()}")
    solver.save(save_dict=solver.split(dm_), path="./")
    os.chdir("..")
```

## Workflow Steps

We will use SeisFlows to run the test. Start a **new** working directory (e.g., *resolution*) and:

1) Copy the parameter file `utils/parameters.yaml` to your new working directory
2) Set *PATH.GLOBAL* to be a directory in your working directory: *resolution/global*
3) Copy gradient `g(m)` files from final model to: *resolution/global/eval1/gradient* 
4) Create empty directory: *resolution/global/eval2/gradient* 
5) Copy your final model `m` to: *PATH.MODEL*
6) Perturb your final velocity model `m+dm`  (see section above)
7) Copy your perturbed model `m+dm` to *PATH.PERTURB*
8) Run `seisflows submit` which will run a forward simulation, misfit calculation, and adjoint simulation to calcualte `g(m+dm)`
9) The resulting Hessian (`H(m)`) will be saved in *PATH.OUTPUT*

