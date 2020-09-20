"""
Generate a simple regular grid for mesh interpolation
"""
import numpy as np


def make_grid(x_min, x_max, dx, y_min, y_max, dy, z_min, z_max, dz, 
              fid="grid.txt"):
    """
    Take the input parameters and generate the grid 
    """
    x_reg = np.arange(x_min, x_max, dx)
    y_reg = np.arange(y_min, y_max, dy)
    x_grid, y_grid = np.meshgrid(x_reg, y_reg)
    x_out = x_grid.flatten()
    y_out = y_grid.flatten()

    with open(fid, "w") as f:
        for z in np.arange(z_min, z_max, dz):
            for x, y in zip(x_out, y_out):
                f.write(f"{x:16.3f}\t{y:16.3f}\t{z:16.3f}\n")
        

def define_regions(region):
    """
    Define fine-mesh grid for shallow structure
    NZ-North Mesh actual dimensions:
        x_min = 171311.859
        x_max = 633468.312
        y_min = 5286952.00
        y_max = 5904085.00

    Ebht19 interpolated by Carl Tape has three levels:
        shallow: 
            zmim=-8000., zmax=2642.156, dx=1, dy=1, dz=0.25
        crust: 
            zmim=-50000., zmax=-8000., dx=2, dy=2, dz=1
        mantle: 
            zmim=-400000., zmax=-50000., dx=8, dy=8, dz=4
    """
    kwargs = {"x_min": 170000.00,
              "x_max": 635000.00,
              "y_min": 5286000.00,
              "y_max": 5905000.00
              }
    
    if region == "shallow":
        kwargs.update({"z_min": -8000.00,
                       "z_max": 2642.156004,
                       "dx": 5 * 1E3, 
                       "dy": 5 * 1E3, 
                       "dz": 5 * 1E3  
                       })
    elif region == "crust":
        kwargs.update({"z_min": -50000.00,
                       "z_max": -8000.00,
                       "dx": 50 * 1E3,  # 2
                       "dy": 50 * 1E3,  # 2
                       "dz": 1 * 1E3  # 1
                       })
    elif region == "mantle":
        kwargs.update({"z_min": -400000.00,
                       "z_max": -50000.00,
                       "dx": 50 * 1E3,  # 8
                       "dy": 50 * 1E3,  # 8
                       "dz": 1 * 1E3  # 4
                       })
    else:
        raise NotImplementedError

    return kwargs


def make_all():
    """
    Main conveneince function to make each of the structured grids
    """
    for region in ["shallow", "crust", "mantle"]:
        make_grid(fid=f"{region}.xyz", **define_regions(region))


if __name__ == "__main__":
    # make_all()
    make_grid(fid="shallow.xyz", **define_regions("shallow"))



