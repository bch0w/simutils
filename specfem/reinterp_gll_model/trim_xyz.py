"""
Tomo files have some 0 values after running xproject_... in SPECFEM3D Cartesian, 
trim those away to retain only spatial values that correspond to non-zero model 
values
"""
import numpy as np

transfer = ["x", "y", "z", "vp", "vs", "rho", "qmu", "qkappa"]


def print_minmax(arr):
    """print min, max values for all columns"""
    print(arr.shape)
    for i in range(arr.shape[1]):
        print(transfer[i], arr[:,i].min(), arr[:,i].max())


for fid in ["og_tomography_model_mantle.xyz", 
            "og_tomography_model_crust.xyz",
            "og_tomography_model_shallow.xyz",]:
    data = np.loadtxt(fid)
    print(fid)
    print_minmax(data)

    # Use Vs as the check for 0 values, should be the same everywhere
    data = data[np.where(data[:,4] != 0)]  # drop 0 values
    data = data[np.where((data[:,0] >= 171E3) & (data[:,0] <= 631E3))]  # xlim
    data = data[np.where((data[:,1] >= 5286E3) & (data[:,1] <= 5902E3))]  # ylim
    print_minmax(data)

    np.savetxt(fid[3:], data, 
               "%10.3f %11.3f %10.3f %8.3f %8.3f %8.3f %8.3f %8.3f")


            
