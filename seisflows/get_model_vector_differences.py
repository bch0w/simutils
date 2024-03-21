"""
For synthetic inversions, find the model difference between each updated model
in the inversion and the target model to determine model "improvement". Works
with a SeisFlows output directory where updated model parameters ONLY are 
stored. We do not care about model parameters that are not updated
"""
import numpy as np
import matplotlib.pyplot as plt
from seisflows.tools.model import Model

diffs = []
pars = ["vsh"]  #, "vsv", "vph", "vpv"]

m_true = Model(path="MODEL_TRUE", parameters=pars).vector


for i in range(0, 8):
    if i == 0:
        path = "MODEL_INIT"
    else:
        path = f"MODEL_0{i}"
    print(path)
    m = Model(path=path, parameters=pars).vector
    diff = np.sqrt((m_true - m)**2).sum()
    diffs.append(diff)

print(diffs)
plt.plot(diffs, "ko-")
plt.savefig("model_diffs.png")
