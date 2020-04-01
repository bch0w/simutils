"""
Create some synthetic topography based on sine functions
"""
import numpy as np


def sin2d(xmin, xmax, dx, ymin, ymax, dy, lambda_x, lambda_y):
    """
    Create a two-dimensional sine function that defines the topography

    :type 
    """
    # Create a regular grid that we will interpolate onto
    x_reg = np.arange(xmin, xmax, dx)
    y_reg = np.arange(ymin, ymax, dy)

    x_grid, y_grid = np.meshgrid(x_reg, y_reg)

    # Create the sine functions with a tunable frequency
    x_sin = np.sin(np.linspace(-np.pi * lambda_x, np.pi * lambda_x, len(x_reg)))
    
