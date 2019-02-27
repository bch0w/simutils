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
    primer = os.path.join(home, "tomo", "primer")

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
                lines[i] = line.replace(old_value, simulation_type)
                print("\t", lines[i].strip())
            elif "SAVE_FORWARD" in line:
                old_value = line.strip().split()[-1]
                lines[i] = line.replace(old_value, save_forward)
                print("\t", lines[i].strip())
            elif ("ATTENUATION " in line) and ("_ATTENUATION" not in line):
                old_value = line.strip().split()[-1]
                lines[i] = line.replace(old_value, attenuation)
                print("\t", lines[i].strip())
                break

        # delete original file, set cursor to beginning
        f.truncate(0)
        f.seek(0)
        f.writelines(lines)
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
    cmtsolution_destination = os.path.join(dir["data"], "CMTSOLUTION")
    if os.path.exists(cmtsolution_destination):
        os.remove(cmtsolution_destination)
    os.symlink(cmtsolution, cmtsolution_destination)

    print("editing Par_file")
    edit_par_file(fid=os.path.join(dir["data"], "Par_file"), choice="forward")

    print("generating forwardrun file")
    fid_in = os.path.join(
        dir["primer"], "simutils", "run_templates", "forward_simulation.sh")
    fid_out = os.path.join(dir["runfolder"], "forwardrun.sh".format(event_id))
    if os.path.exists(fid_out):
        os.remove(fid_out)

    def generate_runforward(fid_in, fid_out, event_id):
        """
        generate a forward run script to be called by sbatch
        """
        with open(fid_in, "r") as f_in:
            lines = f_in.readlines()
        for i, line in enumerate(lines):
            if "${EVENT_ID}" in line:
                lines[i] = line.replace("${EVENT_ID}", event_id)
        with open(fid_out, "w") as f_out:
            f_out.writelines(lines)

    generate_runforward(fid_in, fid_out, event_id)

    print("forward build complete")


def event_id_from_cmt(fid):
    """
    occasionally need to get the event_id from the cmtsolution file
    :param fid:
    :return:
    """
    with open(fid, "r") as f:
        lines = f.readlines()
    for line in lines:
        if "event name" in line:
            return line.strip().split()[-1]


def post_forward():
    """
    UNTESTED
    After a forward simulation, files need to be moved around so that local path
    and output_files can be free for a new forward or adjoint simulation.
    Make sure
    :return:
    """
    import shutil

    dir = dynamic_filenames(choice="forward")
    output = dir["output_files"]
    output_solver = os.path.join(output, "output_solver.txt")

    # check the output_solver.txt file to make sure simulation has finished
    if os.path.exists(output_solver):
        with open(output_solver) as f:
            text = f.read()
        if "End of the simulation" not in text:
            print("Simulation has not finished, disregarding")
            sys.exit()
    else:
        print("No output_solver.txt file, disregarding")
        sys.exit()

    # get event id from cmtsolution
    cmtsolution = os.path.join(out, "CMTSOLUTION")
    if not os.path.exists(cmtsolution):
        print("CMTSOLUTION doesn't exist in OUTPUT_FILES")
        sys.exit()
    event_id = event_id_from_cat(cmtsolution)

    # create storage folder for all the event specific outputs
    storage = os.path.join(output, "STORAGE", event_id)
    if not os.path.exists(storage):
        os.mkdir(storage)

    # start moving files from OUTPUT_FILES/ to STORAGE/${EVENT_ID}
    # .sem? and timestamp* files
    for quantity in ["*.sem?", "timestamp*"]:
        for fid in glob.glob(os.path.join(output, quantity)):
            shutil.move(fid, storage)

    # the random one-off output files
    for quantity in ["starttimeloop.txt", "sr.vtk", "output_list_sources.txt",
                     "output_list_stations.txt", "Par_file", "CMTSOLUTION",
                     "STATIONS", "output_solver.txt"]:
        shutil.move(os.path.join(output, quantity), storage)

    # proc*_save_forward_arrays.bin files
    for quantity in ["proc*_save_forward_arrays.bin", "proc*_absorb_field.bin"]:
        for fid in glob.glob(os.path.join(dir["local_path"], quantity)):
            shutil.move(fid, storage)

    print("post forward complete")


def build_adjoint(event_id):
    """
    UNTESTED
    Prepare the run folder for an adjoint simulation for a given event id
    TO DO
    SEM folder check?
    adjoint choice for dynamic filenames?
    :param event_id: str
    :return:
    """
    dir = dynamic_filenames(choice="adjoint")
    storage = os.path.join(dir["output_files"], "STORAGE", event_id)

    # make sure cmtsolution file is correct, e.g. if another forward run was
    # made between the prerequisite forward for this adjoint...
    # not sure if this is necessary
    print("symlinking CMTSOLUTION")
    cmtsolution = os.path.join(
        dir["primer"], "cmtsolution_files", "{}CMTSOLUTION".format(event_id)
    )
    cmtsolution_destination = os.path.join(dir["data"], "CMTSOLUTION")
    if os.path.exists(cmtsolution_destination):
        os.remove(cmtsolution_destination)
    os.symlink(cmtsolution, cmtsolution_destination)

    print("editing Par_file")
    edit_par_file(fid=os.path.join(dir["data"], "Par_file"), choice="adjoint")

    print("generating adjointrun file")
    fid_in = os.path.join(
        dir["primer"], "simutils", "run_templates", "adjoint_simulation.sh")
    fid_out = os.path.join(dir["runfolder"], "adjointrun.sh".format(event_id))
    if os.path.exists(fid_out):
        os.remove(fid_out)

    def generate_adjointrun(fid_in, fid_out, event_id):
        """
        generate a forward run script to be called by sbatch
        """
        with open(fid_in, "r") as f_in:
            lines = f_in.readlines()
        for i, line in enumerate(lines):
            if "${EVENT_ID}" in line:
                lines[i] = line.replace("${EVENT_ID}", event_id)
        with open(fid_out, "w") as f_out:
            f_out.writelines(lines)

    generate_adjointrun(fid_in, fid_out, event_id)

    print("checking OUTPUT_FILES")
    # check if proc*_save_forward_arrays.bin files are symlinks/ can be removed
    save_forward_check = glob.glob(
        os.path.join(dir["local_path"], "proc*_save_forward_arrays.bin")
    )
    if save_forward_check:
        if os.path.islink(save_forward_check[0]):
            print("removing symlink proc*_save_forward_arrays.bin")
            for fid in save_forward_check:
                os.remove(fid)
            print("symlinking proc**_save_forward_arrays.bin")
            files = glob.glob(
                os.path.join(storage, "proc*_save_forward_arrays.bin"))
            for fid in files:
                destination = os.path.join(dir["local_path"],
                                           os.path.basename(fid))
                os.symlink(fid, destination)
        else:
            print("proc**_save_forward_arrays.bin files are real, please move")

    # check if proc*_absorb_field.bin files can be removed
    absorb_check = glob.glob(
        os.path.join(dir["local_path"], "proc*_save_forward_arrays.bin")
    )
    if absorb_check:
        if os.path.islink(absorb_check[0]):
            print("removing symlink proc*_absorb_field.bin")
            for fid in absorb_check:
                os.remove(fid)
            print("symlinking proc**_absorb_field.bin")
            files = glob.glob(
                os.path.join(storage, "proc*_absorb_field.bin"))
            for fid in files:
                destination = os.path.join(dir["local_path"],
                                           os.path.basename(fid))
                os.symlink(fid, destination)
        else:
            print("proc**_absorb_field.bin files are real, please move")

    print("symlinking SEM/ files")
    sem_folder = os.path.join(dir["runfolder"], "SEM")
    if os.path.exists(sem_folder):
        os.remove(sem_folder)

    input_sem = os.path.join(dir["runfolder"], "INPUT_SEM", event_id)
    os.symlink(input_sem, sem_folder)

    print("symlinking STATIONS_ADJOINT")
    stations_adjoint = os.path.join(dir["data"], "STATIONS_ADJOINT")
    if stations_adjoint:
        os.remove(stations_adjoint)
    input_stations_adjoint = os.path.join(input_sem, "STATIONS_ADJOINT")
    os.symlink(input_stations_adjoint, stations_adjoint)

    print("adjoint build complete")


def post_adjoint():
    """
    UNTESTED
    After an adjoint simulation, files need to be moved around to get ready
    for more simulations, or for model update proceedings.
    TO DO
    should we delete the one off outputs? I think they're the same as forward
    :return:
    """
    import shutil

    dir = dynamic_filenames(choice="forward")
    output = dir["output_files"]
    output_solver = os.path.join(output, "output_solver.txt")

    # check the output_solver.txt file to make sure simulation has finished
    if os.path.exists(output_solver):
        with open(output_solver) as f:
            text = f.read()
        if "End of the simulation" not in text:
            print("Simulation has not finished, disregarding")
            sys.exit()
    else:
        print("No output_solver.txt file, disregarding")
        sys.exit()

    # get event id from cmtsolution and set storage
    cmtsolution = os.path.join(out, "CMTSOLUTION")
    if not os.path.exists(cmtsolution):
        print("CMTSOLUTION doesn't exist in OUTPUT_FILES")
        sys.exit()
    event_id = event_id_from_cat(cmtsolution)
    storage = os.path.join(output, "STORAGE", event_id)

    # start moving files from OUTPUT_FILES/ to STORAGE/${EVENT_ID}
    for fid in glob.glob(os.path.join(output, "timestamp*")):
        shutil.move(fid, storage)

    # the random one-off output files. add an adjoint tag because naming
    for quantity in ["starttimeloop.txt", "sr.vtk", "output_list_sources.txt",
                     "output_list_stations.txt", "Par_file", "CMTSOLUTION",
                     "STATIONS", "output_solver.txt"]:
        quantity = quantity.split('.')[0] + "_adjoint" + quantity.split('.')[1]
        shutil.move(os.path.join(output, quantity), storage)

    # proc*_save_forward_arrays.bin and proc*_absorb_field.bin files
    for quantity in ["proc*_save_forward_arrays.bin", "proc*_absorb_field.bin"]:
        for fid in glob.glob(os.path.join(dir["local_path"], quantity)):
            shutil.move(fid, storage)

    print("post adjoint complete")


def pre_precoondition_sum():
    """
    Summing kernels requires a few folders and symlinks to be set up beforehand
    TO DO
    kernels_list.txt is listed in constants_tomography.h, dynamically get?
    :return:
    """
    dir = directories()
    input_kernels = os.path.join(dir["runfolder"], "INPUT_KERNELS")
    storage = os.path.join(dir["output_files"], "STORAGE")
    kernels_list = os.path.join(dir["runfolder"], "kernels_list.txt")
    output_sum = os.path.join(dir["runfolder"], "OUTPUT_SUM")

    import pdb;pdb.set_trace()

    # check-make directories
    if not os.path.exists(input_kernels):
        os.mkdir(input_kernels)

    if not ospath.exists(output_sum):
        os.mkdir(output_sum)

    # writing event numbers into kernels list while symlinking files to
    # INPUT_KERNELS/ directory
    with open(kernels_list, "r+") as f:
        f.trunacte(0)
        f.seek(0)
        for event in glob.glob(os.path.join(storage, "*")):
            event_id = os.path.basename(event)
            event_dir = os.path.join(input_kernels, event_id)
            os.mkdir(event_dir)
            os.symlink(event, event_dir)
            f.write("{}\n".format(event_id))


if __name__ == "__main__":
    pre_precoondition_sum()
