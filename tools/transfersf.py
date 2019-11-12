"""
This function transfers all the important files from a seisflows run directory
somewhere else for storage (e.g. off a nobackup dir onto a project dir).
Retains output log files and information about the inversion, but no machinery 
to restart the inversion
"""
import os
import glob
import shutil

def check(flist, path, fid=None):
    """
    I was tired of writing this over and over
    """
    if fid:
        # allows wildcards
        if "*" in fid:
            for f in glob.glob(os.path.join(path, fid)):
                if os.path.exists(f):
                    flist.append(f)
        else:
            if os.path.exists(os.path.join(path, fid)):
                flist.append(os.path.join(path, fid))
    else:
        if os.path.exists(path):
            flist.append(path)
        
    return flist


def gather_all(sf_dir, slim=True):
    """
    Gather file ids to transfer, this assumes that all the seisflows
    path names are standard and have not been edited

    :type sf_dir: str
    :param sf_dir: Seisflows run directory path
    :type slim: bool
    :param slim: if True, only copy light files (i.e. .txt, .json etc.), if 
        False, copy all files regardless of size (i.e. .h5, .png, .pdf)
    """
    flist = []
    # written if you want to rename log files, but you can also just transfer
    # everything and sort it out later

    # for fid in ["output.log", "error.log"]:
    #     log = os.path.join(sf_dir, fid)
    #     # if there is are multiple log files combine before transferring
    #     log_prior = log + '_prior'
    #     if os.path.exists(log_prior):
    #         log_temp = log + '_temp'
    #         os.rename(log, log_temp)
    #         with open(log, 'w') as f_out:
    #             for fid_in in [log_prior, log_temp]:
    #                 with open(fid_in) as f_in:
    #                     f_out.write(f_in.read())
    #                 os.remove(fid_in) 
    #     # append paths to a list
    #     flist.append(log) 
   
    # other text files in the main directory that don't have repeats 
    for fid in ["output.log*", "error.log*", "output.optim"]:
        flist = check(flist, sf_dir, fid)
        
    # output directory contains saved states of seisflows, and paths, params
    flist = check(flist, sf_dir, "output.stats") 
    flist = check(flist, os.path.join(sf_dir, "output"), "*.p")
    flist = check(flist, os.path.join(sf_dir, "output"), "*.json")

    # get pyatoa config and outputs
    pyatoa_io = os.path.join(sf_dir, "pyatoa.io")
    flist = check(flist, pyatoa_io, "*.json")
    if not slim:
        flist = check(flist, os.path.join(pyatoa_io, "data"), "*.h5")
        flist = check(flist, os.path.join(pyatoa_io, "figures", "composites"))
    flist = check(flist, os.path.join(pyatoa_io, "figures"), "output_optim.png")

    # get solver data and output files
    solvers = glob.glob(os.path.join(sf_dir, "scratch", "solver", "*"))
    solvers.sort()
    if solvers:
        mainsolver = solvers[0]
    flist = check(flist, mainsolver, "solver.log") 
    for fid in ["Par_file", "STATIONS"]:
        flist = check(flist, os.path.join(mainsolver, "DATA"), fid)
    for fid in ["output_generate_databases.txt", "output_solver.txt", "*.h"]:
        flist = check(flist, os.path.join(mainsolver, "OUTPUT_FILES"), fid)

    return flist


def transfer(transfer_from, transfer_to, slim=True):
    """
    transfer files from a seisflows directory to some other directory,
    retain the path structure
    """
    flist = gather_all(transfer_from, slim)
    for src in flist:
        # get the relative path
        f_rel = os.path.relpath(src, transfer_from)

        # break it into directory and filename
        f_rel_base = os.path.basename(f_rel)
        f_rel_dir = os.path.dirname(f_rel)

        # create the relative path in the destination
        f_dst_dir = os.path.join(transfer_to, f_rel_dir)
        os.makedirs(f_dst_dir, exist_ok=True)
       
        # os.rename only works in the same filesystem, shutil allows for
        # cross-device linking 
        dst = os.path.join(f_dst_dir, f_rel_base)
        print(f"{src} -> {dst}\n") 
        shutil.move(src, dst)


if __name__ == "__main__":
    slim = True
    transfer_from = os.getcwd()   
    transfer_to = os.path.basename(transfer_from)

    # allow for customized filenames
    transfer_to_temp = input(f"filename {transfer_to} okay? ([y]/n): ")
    if transfer_to_temp:
        transfer_to = transfer_to_temp
   
    path_out = ("/scale_wlg_persistent/filesets/project/nesi00263/bchow/"
                "results/seisflows")
     
    transfer_to = os.path.join(
     "/scale_wlg_persistent/filesets/project/nesi00263/bchow/results/seisflows",
     transfer_to)

    os.makedirs(transfer_to, exist_ok=True)
    transfer(transfer_from, transfer_to, slim) 
    
    os.symlink(src=transfer_to, dst=os.path.join(transfer_from, "saved_as"))



    
     
