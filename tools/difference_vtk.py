"""
Simple script to difference VTK files for the same model.
There is no real checking so user needs to be sure that the grid data is same
Assumes that the VTK files were made from the same mesh
"""
import os
import sys
import glob


def read_file(pathname):
    """
    read the file and return header info and data

    :type pathname: str
    :param pathname: full path to the .vtk file to read
    :rtype lines: list
    :return lines: data line by line from readlines()
    :rtype header_dict: dic
    :return header_dict: dictionary with all relevant header information
        from an unstructured_grid vtk file
    """
    with open(pathname, "r") as f:
        lines = f.readlines()
    
    # determine important line numbers, headers are specific to specfem3d output
    for i, line in enumerate(lines):
        if "POINTS" in line:
            points_n = int(line.split()[1])
            points_line = i + 1
        elif "CELL" in line and not "CELL_TYPE" in line:
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
            data_line = i+2

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


def difference_vtk(model_a_fid, model_b_fid, reverse=1, method="subtract", 
                   write=None):
    """
    read each model and scan line by line, difference all necessary values

    the difference is defined as: c = a - b
    to reverse, set reverse == -1

    :type model_?: str
    :param model_?: name of the model file
    :type path: str
    :param path: path holding the models to be read in
    :type reverse: int
    :param reverse: 1 for a-b, -1 for b-a
    :type write: str
    :param write: name of the output file to be written
    :rtype differences: list
    :return differences: list of the differences in values between a and b
    """
    # read files
    model_a, header_dict_a = read_file(model_a_fid)
    model_b, header_dict_b = read_file(model_b_fid)

    # check that the files have the same characteristics before parsing
    for key in header_dict_a.keys():
        if key == "scalars":
            continue
        elif header_dict_a[key] != header_dict_b[key]:
            sys.exit("{} not equal".format(key))

    # parse through models together and separate by len of line, skip header
    differences = []
    start = header_dict_a["data_line"]
    for a, b in zip(model_a[start:-1], model_b[start:-1]):
        try:
            difference = reverse * (float(a.strip()) - float(b.strip()))
            # this will give a percent difference rather than absolute diff
            if method == "divide":
                difference /= float(a.strip())
            differences.append(difference)
        except (ValueError):
            print("value error")

    # Write out the differences to a new file
    if write:
        with open(write, "w") as f:
            f.writelines(model_a[:header_dict_a["data_line"]])
            for diff in differences:
                if diff == 0:
                    f.write("{:13.5f}    \n".format(float(diff)))
                elif abs(diff) > 1:
                    f.write("{:13.5f}    \n".format(float(diff)))
                else:
                    f.write("{:13.10f}    \n".format(float(diff)))
            # Paraview will sometimes get stuck in a loop if there is no newline
            # at the end of the file
            f.write("\n")

    return differences


def dynamic_file_pick(basepath, method="all", globchoice="*"):
    """
    Dynamically pick the VTK files to difference
    """
    dynamic_pick = glob.glob(os.path.join(basepath, f'{globchoice}.vtk'))
    dynamic_pick.sort()

    # If the folder only contains two files
    if len(dynamic_pick) == 2:
        # If "model_init" make it first, if "model_true", make it last
        for i, pick in enumerate(dynamic_pick):
            if "model_init" in pick:
                model_a = pick
                model_b = dynamic_pick[1-i]
            elif "model_true":
                model_b = pick
                model_a = dynamic_pick[1-i]
            else:
                model_a = dynamic_pick[0]
                model_b = dynamic_pick[1]

        # Specfiy the output modle name based on the input model names
        fid_out = "diff_{}_and_{}.vtk".format(
                                    os.path.basename(model_a).split(".")[0],
                                    os.path.basename(model_b).split(".")[0]
                                          )
        fid_out = os.path.join(basepath, fid_out)

        return [model_a], [model_b], [fid_out]
    # If the folder contains multiple .vtk files
    else:
        # List out all the picks for the user
        for i, pick in enumerate(dynamic_pick):
            print(f"{i}: {pick}")
        
        # Select one model and difference all other models from it
        if method == "select_one":
            select = input("Select index to difference from: ")
            select = dynamic_pick[int(select)]
            model_a = [select] * (len(dynamic_pick) - 1)
            model_b, fid_out = [], []
            for pick in dynamic_pick:
                if pick == select:
                    continue
                else:
                    model_b.append(pick)
                    fid_out.append(os.path.join(
                        basepath, "diff_{}_and_{}.vtk".format(
                                        os.path.basename(select).split(".")[0],
                                        os.path.basename(pick).split(".")[0])
                                                    ))
            return model_a, model_b, fid_out
        # Select both models to difference
        elif method == "select":
            model_a = input("Select model_a index: ")
            model_a = dynamic_pick[int(model_a)] 

            model_b = input("Select model_b index: ")
            model_b = dynamic_pick[int(model_b)]

            fid_out = "diff_{}_and_{}.vtk".format(
                                    os.path.basename(model_a).split(".")[0],
                                    os.path.basename(model_b).split(".")[0]
                                                  )
            fid_out = os.path.join(basepath, fidout)

            return [model_a], [model_b], [fid_out]
        # Difference all models from one another
        elif method == "all":
            raise NotImplementedError
        import ipdb;ipdb.set_trace()
            

if __name__ == "__main__":
    basepath = './'
    base = 'vs_eberhart'
    model_a = f"{base}15.vtk"
    model_b = f"{base}19.vtk"
    dynamic_method = "select_one"
    globchoice = "*"

    # If the choice of models doesn't exist, pick based on the files available
    if not os.path.exists(os.path.join(basepath, model_a)):
        model_a, model_b, fid_out = dynamic_file_pick(basepath, dynamic_method,
                                                      globchoice
                                                      )
    else:
        fid_out = ["diff_{}_and_{}.vtk".format(
                                os.path.basename(model_a).split(".")[0],
                                os.path.basename(model_b).split(".")[0]
                                              )]
        model_a = [model_a]
        model_b = [model_b]
    # Difference Vtk files
    for a, b, f in zip(model_a, model_b, fid_out):
        print(f)
        differences = difference_vtk(a, b, reverse=1, write=f)


