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


def dynamic_filenames():
    """
    Names of files and folders are set in the Par_file or in constats.h or in
    constants_tomography.h. Rather than hardcoding them in here, dynamically
    fetch them incase the user decides to change them
    :param choice: str
    :return:
    """
    drc = directories()

    # grab information from the Par_file
    par_file = os.path.join(drc["data"], "Par_file")
    with open(par_file) as f:
        lines = f.readlines()
    for line in lines:
        if "TOMOGRAPHY_PATH" in line:  # e.g. DATA/tomo_files/
            tomo_path = line.strip().split()[-1]
            drc["tomo_path"] = tomo_path
        elif "LOCAL_PATH" in line:  # e.g. OUTPUT_FILES/DATABASES_MPI/
            local_path = line.strip().split()[-1]
            drc["local_path"] = local_path
    return drc


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
            model = "default"
        elif choice == "adjoint":
            simulation_type = "3"
            save_forward = ".false."
            attenuation = ".false."
            model = "default"
        elif choice == "update":
            simulation_type = "1"
            save_forward = ".false."
            attenuation = ".false."
            model = "gll"
  
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
            elif "MODEL " in line and not "_MODEL" in line:
                old_value = line.strip().split()[-1]
                lines[i] = line.replace(old_value, model)
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
    drc = dynamic_filenames()

    # check if all the prerequisite files are included
    cmtsolution = os.path.join(
        drc["primer"], "cmtsolution_files", "{}CMTSOLUTION".format(event_id)
    )
    cmtsolution_check = os.path.exists(cmtsolution)
    tomofiles_check = os.path.exists(drc["tomo_path"])
    databases_check = os.path.exists(drc["local_path"])
    stations_check = os.path.exists(os.path.join(drc["data"], "STATIONS"))
    mesh_check = os.path.exists(os.path.join(drc["runfolder"], "MESH"))

    # Check-exits to make sure the run folder is set correctly
    if not databases_check:
        print("{} does not exist".format(local_path))
        sys.exit()
    if not cmtsolution_check:
        print("{}CMTSOLUTION does not exist in primer".format(event_id))
        sys.exit()
    if not tomofiles_check:
        print("{tomofiles} does not exist in {data}".format(
            tomofiles=drc["tomo_path"], data=drc["data"])
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
    cmtsolution_destination = os.path.join(drc["data"], "CMTSOLUTION")
    if os.path.exists(cmtsolution_destination):
        os.remove(cmtsolution_destination)
    os.symlink(cmtsolution, cmtsolution_destination)

    print("editing Par_file")
    edit_par_file(fid=os.path.join(drc["data"], "Par_file"), choice="forward")

    print("generating forwardrun file")
    fid_in = os.path.join(
        drc["primer"], "simutils", "run_templates", "forward_simulation.sh")
    fid_out = os.path.join(drc["runfolder"], "forwardrun.sh".format(event_id))
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
    working - may need more checks if already run
    After a forward simulation, files need to be moved around so that local path
    and output_files can be free for a new forward or adjoint simulation.
    Make sure
    :return:
    """
    import shutil

    drc = dynamic_filenames()
    output = drc["output_files"]
    output_solver = os.path.join(output, "output_solver.txt")
    input_sem = os.path.join(drc["runfolder"], "INPUT_SEM")

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
    cmtsolution = os.path.join(output, "CMTSOLUTION")
    if not os.path.exists(cmtsolution):
        print("CMTSOLUTION doesn't exist in OUTPUT_FILES")
        sys.exit()
    event_id = event_id_from_cmt(cmtsolution)

    # create storage folder for all the event specific outputs
    storage = os.path.join(output, "STORAGE")
    event_storage = os.path.join(storage, event_id)
    if not os.path.exists(event_storage):
        if not os.path.exists(storage): 
            os.mkdir(storage)
        os.mkdir(event_storage)

    # start moving files from OUTPUT_FILES/ to STORAGE/${EVENT_ID}
    # .sem? and timestamp* files
    for quantity in ["*.sem?", "timestamp*"]:
        for fid in glob.glob(os.path.join(output, quantity)):
            shutil.move(fid, event_storage)

    # the random one-off output files
    for quantity in ["starttimeloop.txt", "sr.vtk", "output_list_sources.txt",
                     "output_list_stations.txt", "Par_file", "CMTSOLUTION",
                     "output_solver.txt"]:
        try:
            shutil.move(os.path.join(output, quantity), event_storage)
        except FileNotFoundError:
            print("{} not found".format(quantity))
            continue
    
    # proc*_save_forward_arrays.bin files
    for quantity in ["proc*_save_forward_arrays.bin", "proc*_absorb_field.bin"]:
        for fid in glob.glob(os.path.join(drc["local_path"], quantity)):
            shutil.move(fid, event_storage)
    
    # create folder in INPUT_SEM/
    event_sem = os.path.join(input_sem, event_id)
    if not os.path.exists(event_sem):
        if not os.path.exists(input_sem):
            os.mkdir(input_sem)
        os.mkdir(event_sem) 

    print("post forward complete")


def build_adjoint(event_id):
    """
    working
    Prepare the run folder for an adjoint simulation for a given event id
    TO DO
    SEM folder check?
    adjoint choice for dynamic filenames?
    symlinking save_forward_arrays was throwing errors randomly, maybe just move
        the files rather than symlinking?
    :param event_id: str
    :return:
    """
    drc = dynamic_filenames()
    storage = os.path.join(drc["output_files"], "STORAGE", event_id)
    
    # make sure cmtsolution file is correct, e.g. if another forward run was
    # made between the prerequisite forward for this adjoint...
    # not sure if this is necessary
    print("symlinking CMTSOLUTION")
    cmtsolution = os.path.join(
        drc["primer"], "cmtsolution_files", "{}CMTSOLUTION".format(event_id)
    )
    cmtsolution_destination = os.path.join(drc["data"], "CMTSOLUTION")
    if os.path.exists(cmtsolution_destination):
        os.remove(cmtsolution_destination)
    os.symlink(cmtsolution, cmtsolution_destination)

    print("editing Par_file")
    edit_par_file(fid=os.path.join(drc["data"], "Par_file"), choice="adjoint")

    print("generating adjointrun file")
    fid_in = os.path.join(
        drc["primer"], "simutils", "run_templates", "adjoint_simulation.sh")
    fid_out = os.path.join(drc["runfolder"], "adjointrun.sh".format(event_id))
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
                lines[i] = line.replace("${EVENT_ID}", event_id+"_adj")
        with open(fid_out, "w") as f_out:
            f_out.writelines(lines)

    generate_adjointrun(fid_in, fid_out, event_id)

    print("checking OUTPUT_FILES")
    # check if proc*_save_forward_arrays.bin files are symlinks/ can be removed
    save_forward_check = glob.glob(
        os.path.join(drc["local_path"], "proc*_save_forward_arrays.bin")
    )
    if save_forward_check:
        if os.path.islink(save_forward_check[0]):
            print("removing symlink proc*_save_forward_arrays.bin")
            for fid in save_forward_check:
                os.remove(fid)
        else:
            print("proc**_save_forward_arrays.bin files are real, please move")
            sys.exit()
    
    print("symlinking proc**_save_forward_arrays.bin")
    files = glob.glob(
        os.path.join(storage, "proc*_save_forward_arrays.bin"))
    for fid in files:
        destination = os.path.join(drc["local_path"],
                                   os.path.basename(fid))
        os.symlink(fid, destination)

    # check if proc*_absorb_field.bin files can be removed
    absorb_check = glob.glob(
        os.path.join(drc["local_path"], "proc*_absorb_field.bin")
    )
    if absorb_check:
        if os.path.islink(absorb_check[0]):
            print("removing symlink proc*_absorb_field.bin")
            for fid in absorb_check:
                os.remove(fid)
        else:
            print("proc**_absorb_field.bin files are real, please move")
            sys.exit()

    print("symlinking proc**_absorb_field.bin")
    files = glob.glob(
        os.path.join(storage, "proc*_absorb_field.bin"))
    for fid in files:
        destination = os.path.join(drc["local_path"],
                                   os.path.basename(fid))
        os.symlink(fid, destination)
    
    print("symlinking SEM/ files")
    sem_folder = os.path.join(drc["runfolder"], "SEM")
    if os.path.exists(sem_folder):
        os.remove(sem_folder)

    input_sem = os.path.join(drc["runfolder"], "INPUT_SEM", event_id)
    os.symlink(input_sem, sem_folder)

    print("symlinking STATIONS_ADJOINT")
    stations_adjoint = os.path.join(drc["data"], "STATIONS_ADJOINT")
    if os.path.exists(stations_adjoint) or os.path.islink(stations_adjoint):
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

    drc = dynamic_filenames()
    output = drc["output_files"]
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
    cmtsolution = os.path.join(output, "CMTSOLUTION")
    if not os.path.exists(cmtsolution):
        print("CMTSOLUTION doesn't exist in OUTPUT_FILES")
        sys.exit()
    event_id = event_id_from_cmt(cmtsolution)
    storage = os.path.join(output, "STORAGE", event_id)

    # start moving files from OUTPUT_FILES/ to STORAGE/${EVENT_ID}
    for fid in glob.glob(os.path.join(output, "timestamp*")):
        dst = os.path.join(storage, os.path.basename(fid))
        shutil.move(fid, dst)

    print("moving output_solver.txt to storage")
    quantity = "output_solver.txt"
    quantity_new = quantity.split('.')[0] + "_adjoint" + quantity.split('.')[1]
    if os.path.exists(quantity):
        shutil.move(os.path.join(output, quantity), 
                    os.path.join(storage, quantity_new)
                   )

    # proc*_save_forward_arrays.bin and proc*_absorb_field.bin files
    # TO DO: move kernels and remove symlinks to forward saves
    print("removing forward array symlinks")
    for quantity in ["proc*_save_forward_arrays.bin", "proc*_absorb_field.bin"]:
        for fid in glob.glob(os.path.join(drc["local_path"], quantity)):
            if os.path.islink(fid):
                os.remove(fid)

    print("moving kernels to storage")
    for fid in glob.glob(os.path.join(drc["local_path"], "*kernel.bin")):
        new_fid = os.path.join(storage, os.path.basename(fid))
        shutil.move(fid, new_fid)

    print("post adjoint complete")


def precondition_sum():
    """
    working
    Summing kernels requires a few folders and symlinks to be set up beforehand
    TO DO
    kernels_list.txt is listed in constants_tomography.h, dynamically get?
    :return:
    """
    drc = directories()
    input_kernels = os.path.join(drc["runfolder"], "INPUT_KERNELS")
    storage = os.path.join(drc["output_files"], "STORAGE")
    kernels_list = os.path.join(drc["runfolder"], "kernels_list.txt")
    output_sum = os.path.join(drc["runfolder"], "OUTPUT_SUM")

    # check-make directories
    if not os.path.exists(input_kernels):
        os.mkdir(input_kernels)

    if not os.path.exists(output_sum):
        os.mkdir(output_sum)

    # writing event numbers into kernels list while symlinking files to
    # INPUT_KERNELS/ directory
    with open(kernels_list, "r+") as f:
        f.truncate(0)
        f.seek(0)
        for event in glob.glob(os.path.join(storage, "*")):
            event_id = os.path.basename(event)
            event_dir = os.path.join(input_kernels, event_id)
            if not os.path.exists(event_dir):
                os.mkdir(event_dir)
            f.write("{}\n".format(event_id))
            for fid in glob.glob(os.path.join(
                                 storage, event_id, "*kernel.bin")):
                dst = os.path.join(event_dir, os.path.basename(fid))
                try:
                    os.symlink(fid, dst)
                except OSError as e:
                    continue


def model_update():
    """
    Specfem's xmodel_update requires a few directories to be made in the 
    OUTPUT_FILES/ directory, and to have summed and smoothed kernels located
    in some input_kernels directory. Attenuation must also be turned off
    in the par_file
    TO DO
    move mode_update pathnames into dynamic filenames and call from here and post
    :return:
    """
    drc = directories()
    output_sum = os.path.join(drc["runfolder"], "OUTPUT_SUM")
    model_update = os.path.join(
        drc["runfolder"], "src", "tomography", "model_update.f90")
    edit_par_file(fid=os.path.join(drc["data"], "Par_file"), choice="update")

    def strip_markers(string):
        """hacky way to remove the markers from 'text/' 
        """
        return string[1:-2]

    # get model_update pathnames
    with open(model_update, "r") as f:
        lines = f.readlines()
    for line in lines:
        if ":: INPUT_KERNELS_DIR_NAME" in line:
            drctry = strip_markers(line.strip().split()[-1])
            input_kernels_dir = os.path.join(drc["output_files"], drctry)
        elif ":: LOCAL_PATH_NEW_NAME" in line:
            drctry = strip_markers(line.strip().split()[-1])
            local_path_new = os.path.join(drc["output_files"], drctry)
        elif ":: OUTPUT_STATISTICS_DIR_NAME" in line:
            drctry = strip_markers(line.strip().split()[-1])
            output_statistics_dir = os.path.join(drc["output_files"], drctry)
    
    # make directories if they don't exit
    for drctry in [input_kernels_dir, local_path_new, output_statistics_dir]: 
        if not os.path.exists(drctry):
            os.mkdir(drctry)
   
    # flush symlinks in INPUT_KERNELS_DIR_NAME/
    old_kerns = glob.glob(os.path.join(input_kernels_dir, "*kernel_smooth.bin"))
    for ok in old_kerns:
        if os.path.islink(ok):
            os.remove(ok)   
 
    # symlink kernels into input_kernels_dir
    kernels = glob.glob(os.path.join(output_sum, "*kernel_smooth.bin"))
    if not kernels:
        sys.exit("Kernels must be smoothed before model update")
    for kernel in kernels:
        kernel_new = os.path.join(input_kernels_dir, os.path.basename(kernel))
        os.symlink(kernel, kernel_new)

    print("ready for model_update")
   

def post_model_update():
    """
    After a model update, we need to clean the run folder for the next forward
    simulations, and make sure the old outputs aren't overwritten by the next
    swath of simulations
    TO DO:
    check-stop if M00_OUTPUT_FILES exists, otherwise you start making edits to 
    the new output_files directory
    """ 
    import shutil
    drc = dynamic_filenames()
    slurm = os.path.join(drc["output_files"], "SLURM")    
    edit_par_file(fid=os.path.join(drc["data"], "Par_file"), choice="update")
    
    # move slurm files into a separate folder
    print("moving slurm* files")
    if not os.path.exists(slurm):
        os.mkdir(slurm)
    slurmfiles = glob.glob(os.path.join(drc["runfolder"], "slurm-*.out"))
    for sfile in slurmfiles:
        sfile_new = os.path.join(slurm, os.path.basename(sfile))
        shutil.move(sfile, sfile_new)
   
    # change the name of OUTPUT_FILES/ if M00_OUTPUT_FILES doesn't exist
    # if it exists, assume that this has already been run and continue
    old_output_files = os.path.join(drc["runfolder"], "M00_OUTPUT_FILES")
    if not os.path.exists(old_output_files):
        print("moving OUTPUT_FILES/ to M00_OUTPUT_FILES/")
        shutil.move(drc["output_files"], old_output_files)
        
        # set up new OUTPUT_FILES/
        os.mkdir(drc["output_files"])
        surface_h = os.path.join(old_output_files, "surface_from_mesher.h")
        values_h = os.path.join(old_output_files, "values_from_mesher.h")
        for fid in [surface_h, values_h]:    
            shutil.copyfile(fid, os.path.join(
                            drc["output_files"], os.path.basename(fid))
                           )
    else:
        query = input("M00_OUTPUT_FILES exists, "
                      "continue to populate new OUTPUT_FILES? [y/(n)]")
        if query != "y":
            sys.exit()
    
    # create new local path
    print("creating new local_path")
    old_local_path = os.path.join(old_output_files, 
                                  os.path.basename(drc["local_path"])
                                 )
    local_path = drc["local_path"]
    try:
        os.mkdir(local_path)
    except OSError:
        pass
    
    # mv attenuation.bin, Database files
    print("copying *attenuation and *Database files to new local_path")
    for tag in ["*attenuation.bin", "*Database"]:
        fids = glob.glob(os.path.join(old_local_path, tag))
        for fid in fids:
            new_fid = os.path.join(local_path, os.path.basename(fid))
            shutil.copyfile(fid, new_fid)
    
    import pdb;pdb.set_trace()
     
    # mv and rename vp_new.bin, vs_new.bin and rho_new.bin files
    # !!! TO DO remove the hardcoded mesh_files_m01 here
    print("moving and renaming mesh files to new local_path")
    mesh_files = os.path.join(old_output_files, "mesh_files_m01")
    for tag in ["*vp_new.bin", "*vs_new.bin", "*rho_new.bin"]:
        fids = glob.glob(os.path.join(mesh_files, tag))
        for fid in fids:
            # get rid of the _new tag when renaming
            new_tag = os.path.basename(fid).split('_')
            new_tag = "{}_{}.bin".format(new_tag[0], new_tag[1])
            new_fid = os.path.join(local_path, new_tag)
            shutil.copyfile(fid, new_fid)
   
    print("ready for xgenerate_databases") 
      

def check_status():
    """
    For all directories in the current path, check what stage of the run cycle
    you're in by looking at the current state of outputs
    """ 
    def read_par_file(fid):
        """similar to edit par_file except only read in the choices
        """
        with open(fid, "r") as f:
            lines = f.readlines()
      
            # Change parameters
            for i, line in enumerate(lines):
                if "SIMULATION_TYPE" in line:
                    simulation_type = line.strip().split()[-1]
                elif "SAVE_FORWARD" in line:
                    save_forward = line.strip().split()[-1]
            if simulation_type == "1":
                if save_forward == ".true.":
                    status = "forward"
                elif save_forward == ".false."
                    status = "update"
            elif simulation_type == "3":
                status = "adjoint"
 
        return status
    
  
    folders = glob.glob("*")
    # check statuses for each run folder in the parent directory
    for folder in folders:
        if os.path.isdir(folder):
            os.chdir(folder)
            drc = directories()
            status = read_par_file(os.path.join(drc["data"], "Par_file")
            if status == "forward":
                output_solver = os.path.join(
                                    drc["output_files"], "output_solver.txt") 
                if os.path.exists(output_solver):
                else:
                    
                
                           

if __name__ == "__main__":
    # simple argument distribution
    available_funcs = ["build_forward", "post_forward", "build_adjoint",
                       "post_adjoint", "precondition_sum", "model_update",
                       "post_model_update"]
    try:
        func = sys.argv[1]
    except IndexError:
        sys.exit("argument 1 must be in\n {}".format(available_funcs))

    if func not in available_funcs:
        sys.exit("argument 1 must be in\n {}".format(available_funcs))

    try:
        event_id = sys.argv[2]
    except IndexError:
        if "build" in func:
            sys.exit("argument 2 must be event_id")

    if func == "build_forward":
        build_forward(event_id)
    elif func == "build_adjoint":
        build_adjoint(event_id)
    elif func == "post_forward":
        post_forward()
    elif func == "post_adjoint":
        post_adjoint()
    elif func == "precondition_sum":
        pre_precondition_sum()
    elif func == "model_update":
        model_update()
    elif func == "post_model_update":
        post_model_update()

