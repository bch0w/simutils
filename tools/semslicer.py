"""
A Python wrapper for the SPECFEM3D source file sem_model_slice.f90
"""


class Semslicer:
    """
    A wrapper class for the model slice capabilities of Specfem, used to
    interpolate the unstructured Specfem grid onto a regular XYZ coordinate
    system which can then be used as a new external tomography file, or for
    plotting capabilities.
    """
    def __init__(self, xmin, xmax, dx, ymin, ymax, dy, zmin, zmax, dz,):
        """

        """

    def set_f90_template(self, path_to_template="./"):
        """
        Set parameters in a template fortan script to determine
        :return:
        """
