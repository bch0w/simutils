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
    y_reg = np.arange(ymin, ymax+dy/2, dy)
    x_grid, y_grid = np.meshgrid(x_reg, y_reg)
    print(f"X: {len(x_reg)} / Y: {len(y_reg)}\n"
          f"X0: {xmin} / Y0: {ymin}\n"
          f"DX: {dx} / DY: {dy}")

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
        x_out = x_grid.flatten()
        y_out = y_grid.flatten()
        z_out = z_sin.flatten()

        # remove negative values from the topography
        z_out = z_out.clip(min=0)

        # Save a general purpose xyz file
        xyz = np.hstack((np.array([x_out]).T, np.array([y_out]).T,
                         np.array([z_out]).T))
        np.savetxt(f"{save}.xyz", xyz, fmt="%.7f\t%.7f\t%.4f")

        # Save the meshfem3D dat file
        np.savetxt(f"{save}.dat", z_out, fmt="%.4f")



if __name__ == "__main__":
    """
    Parameters here:
    Make sure that kx and ky are not too large w.r.t number of points, 
    otherwise you may run into aliasing problems. I didn't do the proper math
    but you should be able to tell by the plot if the function is aliased 
    or not smooth enough
    """
    # xmin = 642975.5657
    # xmax = 785291.8532
    # ymin = 5233238.5952
    # ymax = 5459851.1609
    # dx = 1000.
    # dy = 1000.
    # These are for the synthetic topo
    # xmin = 172.7
    # xmax = 174.5
    # ymin = -43.0
    # ymax = -41.0
    # dx = 0.0075319
    # dy = 0.0075319
    # This is for the real topo
    xmin = 172.7002728
    xmax = 174.5983116
    ymin = -43.0003455
    ymax = -40.8989454
    # ymax = -40.8914135
    dx = dy = 0.0075319
    kx = ky = 5
    max_topo_m = 5000.
    plot = False
    save = "kg_synth_topo"
    create_topo(xmin, xmax, dx, kx, ymin, ymax, dy, ky, max_topo_m, plot, save)


