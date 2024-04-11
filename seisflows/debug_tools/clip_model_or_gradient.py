"""
To be used in SeisFlows Debug model to clip model or gradient values. To use, 
change desired threshold values in `thresholds` and then run:

$ seisflows debug

and copy-paste this script into the debug module. Modified files will be stored
in `path_out`
"""
import os
from seisflows.tools.model import Model

# !!! SET MIN AND MAX VALUES HERE !!!
thresholds = {
        "vp_kernel": {"min_val": -1e-15,
                      "max_val": 1e-15,
                      },
        "vs_kernel": {"min_val": -1e-15,
                      "max_val": 1e-15,
                      }
        }

# Paths are relative to the SeisFlows working directory
path_in = os.path.join("scratch", "eval_grad", "misfit_kernel")
path_out = os.path.join("output", "GRADIENT_CLIPPED")

# Read in and clip the model
m = Model(path=path_in)
for name, threshold in thresholds.items():
    for i, proc in enumerate(m.model[name]):
        m.model[name][i] = proc.clip(threshold["min_val"], threshold["max_val"])
  
# Write the newly clipped model out
if not os.path.exists(path_out):
    os.makedirs(path_out)

m.write(path=path_out)
