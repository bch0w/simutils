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
    """
    kwargs = {"x_min": 171311.859,
              "x_max": 633468.312,
              "y_min": 5286952.00,
              "y_max": 5904085.00
              }
    kwargs = {"x_min": 172000.00,
              "x_max": 620000.00,
              "y_min": 5287000.00,
              "y_max": 5904000.00
              }
    
    if region == "shallow":
        kwargs.update({"z_min": -8000.00,
                       # "z_max": 2642.156004,
                       "z_max": -7000.00,
                       "dx": 10 * 1E3,  # 1
                       "dy": 10 * 1E3,  # 1
                       "dz": 1 * 1E3  # .25
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
    make_all()



