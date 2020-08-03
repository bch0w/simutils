#!/usr/bin/env python3
"""
Simple script to difference VTK files for the same model.
There is no real checking so user needs to be sure that the grid data is same
Assumes that the VTK files were made from the same mesh
"""
import os
import sys
import glob
import numpy as np


def read_file(pathname):
    """
    Read the .vtk file and return header info and data

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


def pick_files(basepath, method="all", globchoice="*", diff_method="subtract"):
    """
    Dynamically pick the VTK files to difference, useful for differencing a 
    single .vtk file from a list of other files. For example when taking
    model differencing from the initial model. Returns lists containing
    model_a, model_b and the output file identifiers.

    :type basepath: str
    :param basepath: path to where all .vtk files live
    :type method: str
    :param method: method for dynamic picking
        1) select: choose model_a and model_b 
        2) select_one: choose model_a or model_b, difference all other files 
                       from the chosen model
    :type globchoice: str
    :param globchoice: custom wildcard for glob search
    :type diff_method: str
    :param diff_method: method for differencing the vtk files, see 
        'difference_vtk()'
    """
    pick_list = glob.glob(os.path.join(basepath, f'{globchoice}.vtk'))
    pick_list.sort()

    # List out all the picks for the user
    for i, pick in enumerate(pick_list):
        print(f"{i}: {pick}")
    print("\n")

    # Select one model and difference all other models from it
    if method == "select_one":
        choice = input("Are you choosing model [a] or [b]?: ")
        assert(choice in ["a", "b"]), "choice must be 'a' or 'b'"

        select = input(f"Select index for model_{choice}: ")
        select = pick_list[int(select)]
        select_tag = os.path.basename(select).split(".")[0]  # e.g. model_0001
        
        # Model_1 will be a list of itself, Model_2 needs to be filled
        model_1 = [select] * (len(pick_list) - 1)
        model_2, fid_out = [], []

        # Don't difference the selection from itself
        pick_list.remove(select)

        # Make sure output file id naming matches the order of selection
        out_dict = {"a": "{d}_{s}_{p}.vtk",
                    "b": "{d}_{p}_{s}.vtk"}
        for pick in pick_list:
            model_2.append(pick)
            pick_tag = os.path.basename(pick).split(".")[0]
            out_tag = os.path.join(basepath, 
                                   out_dict[choice].format(d=diff_method,
                                                           s=select_tag, 
                                                           p=pick_tag)
                                   )
            # Shorten file tags to make them easier to read
            out_tag = out_tag.replace("gradient_", "g")
            out_tag = out_tag.replace("model_", "m")

            fid_out.append(out_tag)
        if choice == "a":
            model_a, model_b = model_1, model_2
        elif choice == "b":
            model_b, model_a = model_1, model_2

        return model_a, model_b, fid_out

    # Select both models to difference
    elif method == "select":
        model_a = input("Select model_a index: ")
        model_a = pick_list[int(model_a)] 

        model_b = input("Select model_b index: ")
        model_b = pick_list[int(model_b)]

        fid_out = "{}_{}_and_{}.vtk".format(
                                diff_method,
                                os.path.basename(model_a).split(".")[0],
                                os.path.basename(model_b).split(".")[0]
                                              )
        fid_out = os.path.join(basepath, fid_out)

        return [model_a], [model_b], [fid_out]
    else:
        raise NotImplementedError


def difference_vtk(model_a_fid, model_b_fid, method="subtract", write=None):
    """
    read each model and scan line by line, difference all necessary values

    the difference is defined as: c = a - b

    :type model_?: str
    :param model_?: name of the model file
    :type path: str
    :param path: path holding the models to be read in
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
            a = float(a.strip())
            b = float(b.strip())
            if method in ["subtract", "divide"]:
                difference = a - b
                # this will give a percent difference rather than absolute diff
                if method == "divide":
                    difference /= a
            # Take the natural log of the the quotient of a and b, this gives
            # to first order approximation, the percent difference. Yoshi said
            # Albert Tarantola said, "always view models in log space"
            elif method == "log":
                difference = np.log(a / b)
            elif method == "poissons_ratio":
                difference = 0.5 * (a**2 - 2 * b**2) / (a**2 - b**2)

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


if __name__ == "__main__":
    """
    PARAMETER SET
    basepath (str): path to the .vtk files
    model_? (str): for manual file picking set in parameters
    dynamic_method (str): method for automatically selecting files
        available - select, select_one, None
    diff_method (str): method for differencing vtk values, 
        available - subtract, divide, log
    globchoice (str): wildcard for dynamic file picking
    """
    basepath = './'
    dynamic = False

    if dynamic:
        globchoice = input("Specific wildcard for selection? [e.g 'model_*']: ")
        if not globchoice:
            globchoice = "*"
        pick_method = input("Selection method? [select, select_one]: ")
        diff_method = input("Method? [subtract, divide, log, poissons_ratio]: ")
    else:
        globchoice = "model_*"
        pick_method = "select_one"
        diff_method = "log"


    # Small reminder to aid in selection process
    if diff_method == "poissons_ratio":
        print(f"{'='*80}\nmodel_A = Vp; model_B = Vs{'='*80}")
    elif diff_method == "log":
        print(f"{'='*80}\n"
              "For net model update / log difference:\n"
              "log(a/b) ~= (a-b) / b; 'a'\n"
              "model_a = FINAL;  model_b = INITIAL\n"
              f"{'='*80}")

    # Dynamic file picking
    model_a, model_b, fid_out = pick_files(basepath, pick_method, globchoice, 
                                           diff_method)

    # Difference VTK files
    for a, b, f in zip(model_a, model_b, fid_out):
        print(f"diff {a} and {b} with '{diff_method}'... {f}")
        differences = difference_vtk(a, b, write=f, method=diff_method)


