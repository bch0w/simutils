"""
Only partially clean the Seisflows run folder, that is delete most of the 
unnecessary files that will not need to be used in the future, but retain
those files that may be useful for future work, or to restart a run if necessary

TO DO:
-get file paths from paths.py file?
"""
import os
import sys
import glob
import shutil


def delete_(path, dryrun=False):
    """
    Check if a path exists, then delete it according to the path type
    Allows wildcards to be entered as the path name but does not sort with glob
    
    :type path: str
    :param path: path to file or directory to be deleted
    :type dryrun: bool
    :param dryrun: if True, no deletions will occur, only print statements
    """
    # accepts wildcards and symlinks
    for i, fid in enumerate(glob.glob(path)):
        if i == 0:
            print(f"{os.path.basename(path)}", end=", ")
        if os.path.exists(fid) or os.path.islink(fid):
            if not dryrun:
                if os.path.isdir(fid):
                    shutil.rmtree(fid)
                elif os.path.islink(fid):
                    os.unlink(fid)
                else:
                    os.remove(fid)
        
        
def clean_solver(dryrun=False, save_mainsolver=True):
    """
    Cleans the current Seisflows solver directory, leaving the main solver

    :type dryrun: bool
    :param dryrun: if True, no deletions will occur, only print statements
    """
    cwd = os.getcwd()
    if not os.path.exists(os.path.join(cwd, "parameters.yaml")):
        raise Exception("No file parameters.yaml, not in a run folder")            

    solver = os.path.join(cwd, 'scratch' ,'solver')
    event_ids = glob.glob(os.path.join(solver, '*'))
    event_ids.sort()

    for i, event in enumerate(event_ids):
        if os.path.islink(event):
            continue
        print(f"scratch/solver/{os.path.basename(event)}")
        if save_mainsolver and i == 0:
            print("\tmainsolver, skipping")
            continue
        print("\tremoving...", end=" ")
        for del_path in [os.path.join(event, 'bin'),
                         os.path.join(event, 'DATA', 'tomo_files'),
                         # os.path.join(event, 'traces'), 
                         # os.path.join(event, 'SEM'),
                         os.path.join(event, 'OUTPUT_FILES', 'timestamp*'),
                         os.path.join(event, 'OUTPUT_FILES', 'DATABASES_MPI')]:
            delete_(del_path, dryrun)
        print("")

                
def clean_main(pyatoa, slurm, dryrun=False):
    """
    Removes figures from pyatoa.io directory
    """
    cwd = os.getcwd()
    if not os.path.exists(os.path.join(cwd, "parameters.py")):
        return

    if pyatoa:
        figs = os.path.join(cwd, 'pyatoa.io', 'figures')
        print("pyatoa.io/figures")
        print("\tremoving...", end=" ")
        delete_(os.path.join(figs, 'm??'), dryrun)
        print("")
    if slurm:
        outputslurm = os.path.join(cwd, 'output.slurm')
        print(os.path.basename(outputslurm))
        print("\tremoving...", end=" ")
        delete_(outputslurm, dryrun)
        print("")

if __name__ == "__main__":
    # Dry run so you don't delete things without checking
    try:
        dryrun = not bool(sys.argv[1])
    except IndexError:
        dryrun = True
        print("DRYRUN")
    
    # USER PARAMETERS
    save_mainsolver = True
    clean_pyatoa = False
    clean_slurm = False   
 
    clean_solver(dryrun, save_mainsolver)
    clean_main(pyatoa=clean_pyatoa, slurm=clean_slurm, dryrun=dryrun)
