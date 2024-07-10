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
            scalars_line = i
            data_line = i+2

    # easier returns in a dictionary
    header_dict = {"points_n": points_n, "points_line": points_line,
                   "cells_n": cells_n, "cells_size": cells_size,
                   "cells_line": cells_line, "cell_types_n": cell_types_n,
                   "cell_types_line": cell_types_line,
                   "points_data_n": point_data_n,
                   "points_data_line": point_data_line, "scalars": scalars,
                   "scalars_line": scalars_line, "data_line": data_line
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
        3) auto: auto search for models based on diff_method
            for each of the given diff methods:
                *log:      select model init as model b, and any other
                           model tags as model a
                *divide:   automatically calculates Vp/Vs ratio for all
                           available Vp and Vs models
                *poissons: automatically calculates Poissons ratio
                           for all available Vp and Vs models

    :type globchoice: str
    :param globchoice: custom wildcard for glob search
    :type diff_method: str
    :param diff_method: method for differencing the vtk files, see 
        'difference_vtk()'
    """
    pick_list = glob.glob(os.path.join(basepath, f'{globchoice}.vtk'))
    pick_list.sort()

    # List out all the picks for the user
    if method != "auto":
        for i, pick in enumerate(pick_list):
            print(f"{i}: {os.path.basename(pick)}")
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

        model_a = [model_a]
        model_b = [model_b]
        fid_out = [fid_out]

    # Auto select models - Assuming we are in the working directory
    elif method == "auto":
        print("Auto selecting models")
        model_a, model_b, fid_out = [], [], []
        if diff_method == "log":
            for tag in ["vpv", "vph", "vsv", "vsh"]:
                model_init = f"model_init_{tag}.vtk"
                assert(os.path.exists(model_init)), f"{model_init} no existo"
            
                # Determine which other models are available besides init
                other_models = glob.glob(model_init.replace("init", "????"))
                other_models.remove(model_init)
                for other_model in sorted(other_models):
                    model_number = other_model.split("_")[1]
                    fid_out_ = f"update_{model_number}_{tag}.vtk"
                    if os.path.exists(fid_out_):
                        print(f"{fid_out_} exists, skipping...")
                        continue
                    model_a.append(other_model)  # FINAL MODEL
                    model_b.append(model_init)  # INITIAL MODEL
                    fid_out.append(fid_out_)
        elif diff_method in ["poissons", "divide"]:
            # We will gather Vp files first, then choose Vs based on Vp names
            vp_files = glob.glob("model_????_vp.vtk")
            assert vp_files, "No Vp files found, cannot calculate ratios"
            for vp_fid in sorted(vp_files):
                model_number = vp_fid.split("_")[1]
                if diff_method == "poissons":
                    fid_out_ = f"ratio_{model_number}_poissons.vtk"
                elif diff_method == "divide":
                    fid_out_ = f"ratio_{model_number}_vpvs.vtk"
                if os.path.exists(fid_out_):
                    print(f"{fid_out_} exists, skipping...")
                    continue
                vs_fid = vp_fid.replace("vp", "vs")
                assert(os.path.exists(vs_fid)), f"No Vs file found for {vp_fid}"
                model_a.append(vp_fid)
                model_b.append(vs_fid)
                fid_out.append(fid_out_)
        elif diff_method == "mu":
            vs_files = glob.glob("model_????_vs.vtk")
            assert vs_files, "No Vs files found, cannot calcualte mu"
            rho_file = glob.glob("*rho*.vtk")[0]
            for vs_fid in sorted(vs_files):
                model_number = vs_fid.split("_")[1]
                fid_out_ = f"modulus_{model_number}_shear.vtk"
                model_a.append(vs_fid)
                model_b.append(rho_file)
                fid_out.append(fid_out_)
    else:
        raise NotImplementedError

    return model_a, model_b, fid_out


def difference_vtk(model_a_fid, model_b_fid, method="subtract", write=None):
    """
    read each model and scan line by line, difference all necessary values

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
            if method in ["subtract", "pct"]:
                difference = a - b
                # this will give a percent difference rather than absolute diff
                if method == "pct":
                    difference /= a
            elif method == "divide":
                difference = a / b
            elif method == "add":
                difference = a + b
            # Take the natural log of the the quotient of a and b, this gives
            # to first order approximation, the percent difference. Yoshi said
            # Albert Tarantola said, "always view models in log space"
            elif method == "log":
                difference = np.log(a / b)
            elif method == "poissons":
                difference = 0.5 * (a**2 - 2 * b**2) / (a**2 - b**2)
            elif method == "mu":
                # Units of GPa iff a~[m/2] and b~[kg/m**3]
                difference = 1E-9 * (a**2 * b)

            differences.append(difference)
        except ValueError:
            print("value error")

    # Write out the differences to a new file
    if write:
        # Replace the name of the data so paraview identifies it differently
        old_value = header_dict_a["scalars"]
        scalars_line = header_dict_a["scalars_line"]
        # Kinda ugly but e.g. model_0001_vp -> model_0001_log
        new_value = "_".join(old_value.split("_")[:-1] + [method])
        model_a[scalars_line] = model_a[scalars_line].replace(old_value, 
                                                              new_value)

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


def manipulate_vtk(model_a_fid, c, method="subtract", write=None):
    """
    read each model and scan line by line, add, subtract, multiply, divide

    :type model_?: str
    :param model_?: name of the model file
    :type c: float
    :param c: value to manipulate vtk file with
    :type path: str
    :param path: path holding the models to be read in
    :type write: str
    :param write: name of the output file to be written
    :rtype differences: list
    :return differences: list of the differences in values between a and b
    """
    # read files
    model_a, header_dict_a = read_file(model_a_fid)

    # parse through models together and separate by len of line, skip header
    start = header_dict_a["data_line"]
    # Convert from strings to floats
    differences = np.array(list(map(float, model_a[start:-1])))
    c = float(c)

    if method == "subtract":    
        differences -= c
    elif method == "divide":
        differences /= c
    elif method == "multiply": 
        differences *= c
    elif method == "add":
        differences += c
    elif method == "norm":
        import ipdb;ipdb.set_trace()
        c = max([abs(max(differences)), abs(min(differences))])
        differences /= c

    # Write out the differences to a new file
    if write:
        # Replace the name of the data so paraview identifies it differently
        old_value = header_dict_a["scalars"]
        scalars_line = header_dict_a["scalars_line"]
        # Kinda ugly but e.g. model_0001_vp -> model_0001_log
        new_value = "_".join(old_value.split("_")[:-1] + [method])
        model_a[scalars_line] = model_a[scalars_line].replace(old_value, 
                                                              new_value)

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

def print_header(diff_method):
    """
    Display some information that should be relevant to the user in terms of how
    to use the various options.
    """
    print("=" * 80)
    if diff_method == "poissons":
        print(f"For Poissons ratio: model_a = Vp; model_b = Vs")
    elif diff_method == "log":
        print("For net model update / log difference:\n"
              "\nlog(a / b) ~= (a - b) / b\n"
              "\nmodel_a = FINAL;  model_b = INITIAL\n")
    elif diff_method == "divide":
        print("For division: diff = model_a / model_b")
    elif diff_method == "subtract":
        print("For subtraction: diff = model_a - model_b")
    elif diff_method == "pct":
        print("For percentage difference: diff = "
              "(model_a - model_b) / model_a")
    elif diff_method == "mu":
        print("For shear modulus mu: diff = a**2 * b\n"
              "model_a = Vs; model_b = rho (density)")
    print("=" * 80)


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
    basepath = os.getcwd()
    dynamic = True

    if dynamic:
        globchoice = "*"
        # globchoice = input("Specific wildcard for selection? [e.g 'model_*']: ")
        # if not globchoice:
        #     globchoice = "*"

        # Choose which method for picking files to diff. Allow both string and
        # index choosing of method
        available_pick = ["auto", "select", "select_one", "constant"]
        pick_method = input(f"Selection method? {available_pick}: ")
        try:
            pick_method = available_pick[int(pick_method)]
        except ValueError:
            assert(pick_method in available_pick), f"unknown: {pick_method}"
            pass

        # Choose which method for diff'ing files, allow string and index choice
        if pick_method != "constant":
            available_diff = ["log", "poissons", "divide", "pct", "subtract", 
                              "mu", "add"]
        else:
            available_diff = ["add", "subtract", "divide", "multiply", "norm"]
        diff_method = input(f"Method? {available_diff}: ")
        try:
            diff_method = available_diff[int(diff_method)]
        except ValueError:
            assert(diff_method in available_diff), f"unknown: {diff_method}"
            pass
    else:
        globchoice = "model_*"
        pick_method = "select_one"
        diff_method = "log"


    if pick_method != "constant":
        if pick_method != "auto":
            print_header(diff_method)

        # Dynamic file picking
        model_a, model_b, fid_out = pick_files(basepath, pick_method, globchoice, 
                                               diff_method)

        # Difference VTK files
        for a, b, f in zip(model_a, model_b, fid_out):
            print(f"diff {a} and {b} with '{diff_method}'... {f}")
            differences = difference_vtk(a, b, write=f, method=diff_method)
    else:
        model_a = input("model_a fid: ")
        constant = input("constant: ")
        fid_out = input("fid_out: ")
        manipulate_vtk(model_a, constant, diff_method, fid_out)


