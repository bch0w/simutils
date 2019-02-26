#!/usr/bin/env python3
"""
Simulation Utilities
A collection of utilities for working with Specfem3D Cartesian on an HPC cluster
These scripts and functions glue together the inputs and outputs of various
Fortran codes executed by Specfem3D. Originally a set of bash scripts to move
files around and create run-scripts, I thought it would be best to do everything
in Python as Bash is finicky.
"""
from __future__ import print_function

import os
import sys
import glob


def directories():
    """
    Return a collection of paths that most of the functions will want.
    These are hardcoded and should be standard among runs.
    :return: dict
    """
    # directories with user-information
    home = os.path.expanduser("~")
    primer = os.path.join(home, "primer")

    # current specfem3d master folder
    runfolder = os.getcwd()
    data = os.path.join(runfolder, "DATA")
    output_files = os.path.join(runfolder, "OUTPUT_FILES")

    directories_dictionary = {"home": home, "primer": primer, "data": data,
                              "runfolder": runfolder,
                              "output_files": output_files}

    return directories_dictionary


def dynamic_filenames(choice):
    """
    Names of files and folders are set in the Par_file or in constats.h or in
    constants_tomography.h. Rather than hardcoding them in here, dynamically
    fetch them incase the user decides to change them
    :param choice: str
    :return:
    """
    dir = directories()

    # grab information from the Par_file
    if choice == "forward":
        par_file = os.path.join(dir["data"], "Par_file")
        with open(par_file) as f:
            lines = f.readlines()
        for line in lines:
            if "TOMOGRAPHY_PATH" in line:  # e.g. DATA/tomo_files/
                tomo_path = line.strip().split()[-1]
                dir["tomo_path"] = tomo_path
            elif "LOCAL_PATH" in line:  # e.g. OUTPUT_FILES/DATABASES_MPI/
                local_path = line.strip().split()[-1]
                dir["local_path"] = local_path
        return dir

    # grab information from
    if choice == "adjoint":
        return


def edit_par_file(fid, choice):
    """
    Edits the Par_file for forward or adjoint simulations
    Encompasses the tasks of change_simulation_type.pl, with a few extras
    Assumes the Par_file will always be the same format (e.g. sim_type   = 1),
    but I don't think Specfem devels will change it during my PhD so okay.
    :param fid: str
    :param choice: str
    """
    with open(fid, "r+") as f:
        lines = f.readlines()
        if choice == "forward":
            simulation_type = "1"
            save_forward = ".true."
            attenuation = ".true."
        elif choice == "adjoint":
            simulation_type = "3"
            save_forward = ".false."
            attenuation = ".false."

        # Change parameters
        for i, line in enumerate(lines):
            if "SIMULATION_TYPE" in line:
                old_value = line.strip().split()[-1]
                lines[i].replace(old_value, simulation_type)
                print(lines[i].strip())
            elif "SAVE_FORWARD" in line:
                old_value = line.strip().split()[-1]
                lines[i].replace(old_value, save_forward)
                print(lines[i].strip())
            elif "ATTENUATION " in line and "_ATTENUATION" not in line:
                old_value = line.strip().split()[-1]
                lines[i].replace(old_value, attenuation)
                print(lines[i].strip())

        fid.writelines(lines)
        print("overwrote Par_file")


def build_forward(event_id):
    """
    Prepare the run folder for a forward simulation for a given event id
    TO DO:
    + add a check for proc**_Database files to see if mesher was run
    + add a check for proc**_external_mesh.bin files to see if generate_db run
    + tomo_files isn't really required past model update 1?
    :param event_id: str
    :return:
    """
    dir = dynamic_filenames(choice="forward")

    # check if all the prerequisite files are included
    cmtsolution = os.path.join(
        dir["primer"], "cmtsolution_files", "{}CMTSOLUTION".format(event_id)
    )
    cmtsolution_check = os.path.exists(cmtsolution)
    tomofiles_check = os.path.exists(dir["tomo_path"])
    databases_check = os.path.exists(dir["local_path"])
    stations_check = os.path.exists(os.path.join(dir["data"], "STATIONS"))
    mesh_check = os.path.exists(os.path.join(dir["runfolder"], "MESH"))

    # Check-exits to make sure the run folder is set correctly
    if not databases_check:
        print("{} does not exist".format(local_path))
        sys.exit()
    if not cmtsolution_check:
        print("{}CMTSOLUTION does not exist in primer".format(event_id))
        sys.exit()
    if not tomofiles_check:
        print("{tomofiles} does not exist in {data}".format(
            tomofiles=dir["tomo_path"], data=dir["data"])
        )
        sys.exit()
    if not stations_check:
        print("STATIONS does not exist in DATA/")
        sys.exit()
    if not mesh_check:
        print("MESH/ does not exist in RUNFOLDER")
        sys.exit()

    # set up the run folder
    print("symlinking CMTSOLUTION")
    os.symlink(cmtsolution, os.path.join(dir["data"],"CMTSOLUTION"))

    print("editing Par_file", end=".")
    edit_par_file(fid=os.path.join(dir["DATA"], "Par_file"), choice="forward")

    print("generating RUNFORWARD file")
    fid = os.path.join(
        dir["primer"], "simutils", "run_templates", "forward_simulation.sh")
    with open(fid, "r") as f:
        lines = f.readlines()


