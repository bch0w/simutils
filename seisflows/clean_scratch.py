"""
Partially clean a SeisFlows scratch directory while keeping the main components 
intact. This allows for maximum memory conservationw without removing the 
possibility of re-starting an inversion using the same directory with minimal 
required work.
"""
import os
from glob import glob
import sys
import shutil
import argparse


def _delete(path):
    """
    Check if a path exists, then delete it according to the path type
    Allows wildcards to be entered as the path name but does not sort with glob
    
    :type path: str
    :param path: path to file or directory to be deleted
    """
    # accepts wildcards and symlinks
    for i, fid in enumerate(glob(path)):
        if i == 0:
            print(f"\t{os.path.basename(path)} x {len(glob(path))}")
        if os.path.exists(fid) or os.path.islink(fid):
            if os.path.isdir(fid):
                shutil.rmtree(fid)
            elif os.path.islink(fid):
                os.unlink(fid)
            else:
                os.remove(fid)


def _dry_delete(path):
    """
    A version of delete that simply prints the filename that it will delete.
    Explicitely a different function to keep things distinct

    :type path: str
    :param path: path to file or directory to be dry deleted
    """
    # accepts wildcards and symlinks
    for i, fid in enumerate(glob(path)):
        if i == 0:
            print(f"\t{os.path.basename(path)} x {len(glob(path))}")

        
def clean(args):
    """
    Cleans the current Seisflows solver directory

    :type args: argparse.ArgumentParser
    :param args: command line arguments that control what parts of the directory
        to delete and whether or not this is a dry run
    """
    assert(os.path.exists(os.path.join(args.path, "parameters.yaml"))), \
        "Parameter file not found, not in a SeisFlows run folder"

    # Be explicit on whether or not this is a dry run
    if args.run:
        del_fx = _delete
    else:
        print(">>>>>>>>>> DRY RUN <<<<<<<<<<")
        del_fx = _dry_delete

    if args.evalfunc:
        evalfunc = os.path.join(args.path, "scratch", "evalfunc", "*")
        for fid in evalfunc:
            del_fx(fid) 

    if args.evalgrad:
        evalgrad = os.path.join(args.path, "scratch", "evalgrad", "*")
        for fid in evalgrad:
            del_fx(fid) 

    event_ids = sorted(glob(os.path.join(args.path, "scratch", "solver", "*")))

    for i, event in enumerate(event_ids):
        # Skip over the 'mainsolver' symlink
        if os.path.islink(event):
            continue
        print(f"EVENT ID: {os.path.basename(event)}")

        if i == 0 and not args.mainsolver:
            print("\tThis is the mainsolver... skipping")
            continue

        # Collect all the files to be deleted
        output_files = os.path.join(event, "OUTPUT_FILES", "{}")
        databases_mpi = os.path.join(event, "OUTPUT_FILES", 
                                     "DATABASES_MPI", "{}")
        to_delete = []

        if args.solver_logs:
            for wildcard in ["timestamp*", "output_list_s*.txt", "sr.vtk", 
                             "starttimeloop.txt"]:
                to_delete.append(output_files.format(wildcard))

        if args.databases_mpi:
            to_delete.append(databases_mpi.format("*"))
        else:
            if args.vtk_files:
                to_delete.append(databases_mpi.format("*.vt?"))

            if args.save_forward_arrays:
                for wildcard in ["proc*_save_forward_arrays.bin", 
                                 "proc*_absorb_field.bin"]:
                    to_delete.append(databases_mpi.format(wildcard))

        for d in to_delete:
            del_fx(d)
            
        

                
if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-p", "--path", type=str, help="path to the SeisFlows "
                        "working directory; defaults to cwd", 
                        default=os.getcwd())
    parser.add_argument("-R", "--run", action="store_true", 
                        help="actually remove files, if not used, then default "
                        "will be a dry run", default=False)
    parser.add_argument("-M", "--mainsolver", action="store_true",
                        help="include the mainsolver in the deletions, usually "
                        "not adviseable as the mainsolver can be used to "
                        "repopulate the other solver directories if a restart "
                        "is required")
    parser.add_argument("-f", "--evalfunc", action="store_true", help="remove "
                        "the contents of the evalfunc scratch directory",
                        default=False)
    parser.add_argument("-g", "--evalgrad", action="store_true", help="remove "
                        "the contents of the evalgrad scratch directory",
                        default=False)
    parser.add_argument("-l", "--solver_logs", action="store_true", 
                        help="remove any misc. log files inthe solver dirs "
                        "these include 'timestamp', 'sr.vtk', 'starttimeloop', "
                        "'output_list_sources/stations'. This is mostly to "
                        "keep file count down, not memory", default=False)
    parser.add_argument("-v", "--vtk_files", action="store_true",
                        help="remove all .vtk and .vtu files that have been "
                        "generated by SPECFEM and saved to DATABASES_MPI",
                        default=False)
    parser.add_argument("-d", "--databases_mpi", action="store_true",
                        help="remove all files inside the DATABASES_MPI "
                        "directory. This will probably prevent you from "
                        "restarting the inversion without considerable manual "
                        "effort, not advised")
    parser.add_argument("-s", "--save_forward_arrays", action="store_true",
                        help="remove the files generated when saving arrays "
                        "after a forward simulation, these include the "
                        "'proc???_save_forward_arrays.bin' and "
                        "'proc???_absorb_field.bin' files, these are probably "
                        "the largest memory sink")

    args = parser.parse_args()
        
    clean(args)    

