"""
08/31/22 SPECFEM2D 7e019b62

Tape2007 Checkerboard model (EXAMPLES/Tape2007/DATA/model_velocity.dat_checker)
does not build correctly in SPECFEM2D. The underlying velocity model (which 
shows 4 large lobes) gets jumbled during the read as coordinates are mismatched.

This function corrects this by finding out SPECFEM2D's internally assumed
coordinate structure, and matching values of the input model one by one.

.. note::
    You will need to run the example with Par_file parameter `SAVE_MODEL`=='gll'
    in order to get the files 'proc000000_[xz].bin'. The GLL model has the best 
    precision. The ASCII output model (.dat) rounds off decimals in comparison.
"""
import numpy as np
from math import isclose
from seisflows.tools.specfem import read_fortran_binary

# Truncate at 0 because thankfully the model coordinates are very unique so
# we only need integer values to match
TRUNCATE = 0

def trunc(values, decs=0):
    """https://stackoverflow.com/questions/42021972/
        truncating-decimal-digits-numpy-array-of-floats
    """
    return np.trunc(values*10**decs)/(10**decs)

input_model = np.genfromtxt("model_velocity.dat_checker")
_, x_input, z_input, rho_input, vp_input, vs_input = input_model.T

# Truncate floats for easier matching
x_input = trunc(x_input, decs=TRUNCATE)
z_input = trunc(z_input, decs=TRUNCATE)

# These define SPECFEM's assumed internal coordinate structure
x = read_fortran_binary("proc000000_x.bin")
z = read_fortran_binary("proc000000_z.bin")

# Building a new model from the coordinates of the 'internal' model, matched
# against the values of the 'input' model
rho, vp, vs = [], [], []  # the newly ordered model values

# Looping through internal model assumed coordinate structure
for i, (x_, z_) in enumerate(zip(x, z)):
    x_ = trunc(x_, TRUNCATE)
    z_ = trunc(z_, TRUNCATE)

    # Find line number corresponding to x, z coordinates
    # !!! There are doubles in the input model, only take first, assuming
    # !!! that they are the SAME
    idx = np.where((x_input == x_) & (z_input == z_))[0][0]
    rho.append(rho_input[idx])
    vp.append(vp_input[idx])
    vs.append(vs_input[idx])

# Write out the new model
_, x_out, z_out, *_ = input_model.T
with open("model_velocity.dat_checker_UPDATED", "w") as f:
    for i in range(len(x_out)):
        f.write(f"\t{i}\t{x_out[i]:10.4f}\t{z_out[i]:10.4f}\t"
                f"{rho[i]:8.4f}\t{vp[i]:8.4f}\t{vs[i]:8.4f}\n")

