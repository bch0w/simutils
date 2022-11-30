def read_specfem_vtk(path_to_vtk):
    """
    Removed from Pyatoa.utils.read, might be useful somewhere else later    

    Read the unstructured grid VTlK files that are output by Specfem for model,
    gradient, and kernel visualizations. Returns a header as a dictionary, and
    the lines of the data file. Useful for manipulating VTK files in place.

    :type path_to_vtk: str
    :param path_to_vtk: full path to the .vtk file to read
    :rtype: tuple (list, dict)
    :return: first index of tuple is a list read line by line from readlines(),
        second index of tuple is a dictionary with all relevant header 
        information from an unstructured_grid vtk file
    """
    with open(path_to_vtk, "r") as f:
        lines = f.readlines()

    # determine important line numbers, headers are specific to specfem3d output
    for i, line in enumerate(lines):
        if "POINTS" in line:
            points_n = int(line.split()[1])
            points_line = i + 1
        elif "CELL" in line and "CELL_TYPE" not in line:
            cells_n = int(line.strip().split()[1])
            cells_size = int(line.strip().split()[2])
            cells_line = i
        elif "CELL_TYPES" in line:
            cell_types_n = int(line.strip().split()[1])
            cell_types_line = i
        elif "POINT_DATA" in line:
            point_data_n = int(line.split()[1])
            point_data_line = i
        elif "SCALARS" in line:
            scalars = line.strip().split()[1]
            data_line = i + 2

    # easier returns in a dictionary
    header_dict = {"points_n": points_n, "points_line": points_line,
                   "cells_n": cells_n, "cells_size": cells_size,
                   "cells_line": cells_line, "cell_types_n": cell_types_n,
                   "cell_types_line": cell_types_line,
                   "points_data_n": point_data_n,
                   "points_data_line": point_data_line, "scalars": scalars,
                   "data_line": data_line
                   }

    return lines, header_dict

