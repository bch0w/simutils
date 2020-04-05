"""
Create some synthetic topography based on a 2d sine function
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm



def create_topo(xmin, xmax, dx, kx, ymin, ymax, dy, ky, max_topo_m=1, plot=True, 
                save=None):
    """
    Create a two-dimensional sine function that defines the topography

    :type xmin: float
    :param xmin: minimum x-dimension of the topography
    :type xmax: float
    :param xmax: maximum x-dimension of the topography
    :type dx: float
    :param dx: grid spacing in the x direction
    :type ymin: float
    :param ymin: minimum y-dimension of the topography
    :type ymax: float
    :param ymax: maximum y-dimension of the topography
    :type dy: float
    :param dy: grid spacing in the y direction
    :type lambda_x: int
    :param lambda_x: number of wavelengths of topography across the x dimension
    :type lambda_y: int
    :param lambda_y: number of wavelengths of topography across the y dimension
    """
    def sin2d(x_, y_, kx_, ky_):
        """The 2D sine function"""
        return max_topo_m * np.sin(kx_ * x_) * np.sin(ky_ * y_)

    # Create a regular grid based on mesh dimensions
    x_reg = np.arange(xmin, xmax, dx)
    y_reg = np.arange(ymin, ymax, dy)
    x_grid, y_grid = np.meshgrid(x_reg, y_reg)

    # Create the sine functions with given number of oscillations
    z_sin = sin2d(x_grid, y_grid, 2*np.pi*kx/(xmax-xmin), 
                  2*np.pi*ky/(ymax-ymin))

    # Quick plot the topography for easier interpretation
    if plot:
        f = plt.figure()
        ax = f.gca(projection="3d")
        ax.plot_surface(x_grid, y_grid, z_sin, cmap=cm.coolwarm)
        plt.show()

    # Save the file into a format that Meshfem can recognize
    if save:
        z_out = z_sin.flatten()
        np.savetxt(f"{save}.dat", z_out, fmt="%.2f")


if __name__ == "__main__":
    """
    Parameters here:
    Make sure that kx and ky are not too large w.r.t number of points, 
    otherwise you may run into aliasing problems. I didn't do the proper math
    but you should be able to tell by the plot if the function is aliased 
    or not smooth enough
    """
    xmin = ymin = 0
    xmax = ymax =  1000
    dx = dy = 1
    kx = ky = 5
    max_topo_m = 200.
    plot = True
    save = "test"
    create_topo(xmin, xmax, dx, kx, ymin, ymax, dy, ky, max_topo_m, plot, save)


