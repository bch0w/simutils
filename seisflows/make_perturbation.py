"""
To be copy-pasted into the SeisFlows debug mode to generate the correct 
perturbation files for use in the point spread tests
"""
import os
assert(os.path.basename(os.getcwd()) == "original"), "wrong directory"
vp_offset = 3000
vs_offset = 1500
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
    

