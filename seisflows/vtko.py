"""
An automation script to take the output models from Seisflows, convert the .bin
files into .vtk files, collect into approriate directories, create differences
and plot standard visualizations using Mayavi
"""
import os
import sys
import time
import subprocess
from glob import glob
from pprint import pprint
from pyatoa.utils.write import rcv_vtk_from_specfem, src_vtk_from_specfem


class VTKO:
    """
    VTK Outputter. A class used to automatically generate .vtk files from the
    SeisFlows output directory on a Slurm system. Organizes .bin files and
    symlinks them into the appropriate SPECFEM run directory,
    runs xcombine_vol_data_vtk, and then moves and organizes the output
    .vtk files in appropriate directories for easy access.
    """
    def __init__(self, bash_script, cwd=None, specfem=None, dirout=None):
        """
        Set up some internal path and parameter structure. Typical.
        """
        self.bash_script = bash_script

        # Determine the SeisFlows working directory
        if not cwd:
            self.cwd = os.getcwd()
        else:
            self.cwd = cwd

        # Find a suitable SPECFEM run directory
        if not specfem:
            self.specfem = os.path.join(self.cwd, "scratch", "solver",
                                        "mainsolver")
        else:
            self.specfem = specfem
        assert(os.path.exists(os.path.join(self.specfem, "bin"))), \
            "SPECFEM directory must have a `bin/` directory"
        assert(os.path.exists(os.path.join(self.specfem, "bin",
                                           "xcombine_vol_data_vtk"))), \
            "The SPECFEM bin/ dir must have the binary `xcombine_vol_data_vtk`"

        # Set the directory where the output files will be placed
        if not dirout:
            self.dirout = os.path.join(self.cwd, "scratch", "preprocess", "vtk")
        else:
            self.dirout = dirout
        if not os.path.exists(self.dirout):
            os.makedirs(self.dirout)

        # Attributes that should be populated by functions
        self.available = None
        self.pars = None

    def setup(self):
        """
        Setup a temporary 'SUM' directory that will hold symlinks to
        bin files. Ensure that if a 'SUM' directory already exists, that the
        user is notified incase files there should not be overwritten. Option
        to clean out all files inside 'SUM' directory for a clean start
        """
        sum_dir = os.path.join(self.specfem, "SUM")
        if not os.path.exists(sum_dir):
            os.makedirs(sum_dir)

        fids = glob(os.path.join(sum_dir, "*"))
        if fids:
            check = input("files found in SUM directory, clean? [y/(n)]: ")
            if check == "y":
                for f in fids:
                    os.remove(f)
            else:
                sys.exit("will not overwrite SUM directory, exiting")

    def discover(self):
        """
        Discover the available bin files within the SeisFlows output directory.
        The options are models, gradients and kernels
        """
        output = os.path.join(self.cwd,  "output")
        assert(os.path.exists(output)), f"output path {output} does not exist"

        def split_tag(t):
            """split the output seisflows tag between the proc_ and .bin"""
            t = os.path.basename(t)
            t = "_".join(t.split("_")[1:])  # e.g. mu_kernel.bin
            return "".join(t.split(".")[:-1])

        print("discovering available models/gradients/kernels")
        all_files = glob(os.path.join(output, "*_????"))

        # Figure out the labels of the directories, e.g. models, gradients etc.
        labels = set([os.path.basename(_).split("_")[0] for _ in all_files])
        available = {t: {} for t in labels}
        for lab in labels:
            # kernel directory must be treated differently
            if lab == "kernels":
                kernels = glob(os.path.join(output, "kernels_*"))  # e.g. kernels_0001
                kernel_tags = [os.path.basename(_) for _ in kernels]

                # !!! Assuming here that all kernels have the same events/params
                events = glob(os.path.join(kernels[0], "*"))
                event_tags = [os.path.basename(_) for _ in events]
                # summed kernels are already represented by the gradient
                for s in ["sum", "sum_nosmooth"]:
                    try:
                        event_tags.remove(s)
                    except ValueError:
                        continue

                # Determine the unique event and parameters used
                bin_files = glob(os.path.join(events[0], "*"))
                params = set([split_tag(_) for _ in bin_files])

                available[lab]["tags"] = kernel_tags
                available[lab]["events"] = event_tags
                available[lab]["params"] = list(params)

            else:
                directories = glob(os.path.join(output, f"{lab}_*"))
                tags = [os.path.basename(_) for _ in directories]

                bin_files = glob(os.path.join(directories[0], "*"))
                params = set([split_tag(_) for _ in bin_files])

                # Check that these gradients/models have not been made yet
                # Get a list of existing vtk files in the output directory
                check = [os.path.basename(_).split('.')[0] for _ in
                         glob(os.path.join(self.dirout, lab, '*'))]

                # Here we assume that if one parameter exists, then all exit
                for t in tags:
                    for p in params:
                        if f"{t}_{p}" in check:
                            print(f"\t{t}_{p} already exists, excluding...")
                            tags.remove(t)
                            break

                available[lab]["tags"] = tags
                available[lab]["params"] = list(params)

        if not available:
            sys.exit("no available models/gradients/kernels found")
        else:
            for lab in available.keys():
                print(f"{l}: {len(available[l]['tags'])}", end=" / ")
            print("")

        # Allow the user to manually edit the `available` before it is set as
        # a class attribute
        pprint(available)
        check = input("Manual remove some files? y/[n]: ")
        if check == "y":
            # File tags are located in `available`
            import ipdb; ipdb.set_trace()

        self.available = available

    def organize_models_gradients(self):
        """
        Models and gradients are organized by parameter and model number
        """
        c = 0
        for unique_id in self.pars["tags"]:  # e.g. model_0001
            directory = os.path.join("output", unique_id)
            for par in self.pars["params"]:
                srcs = glob(os.path.join(directory, f"proc*_{par}*.bin"))
                for src in srcs:
                    # e.g. proc000000_vs.bin
                    filename = os.path.basename(src)
                    # Add a unique tag to the filename, different if kernels
                    parts = filename.split("_")
                    # Rejoin end of string if e.g. 'vs_kernel' and not 'vs'
                    parts[1] = "_".join(parts[1:])
                    new_filename = "_".join([parts[0], unique_id, parts[1]])
                    # Set the source and destination for symlink
                    dst = os.path.join(self.specfem, "SUM", new_filename)
                    os.symlink(os.path.abspath(src), dst)
                    c += 1
        return c

    def organize_kernels(self, kernel):
        """
        Kernels are organized by event id, so require a further level of
        organization with respect to gradient and model directories

        :type kernel: str
        :param kernel: the kernel name
        """
        c = 0
        # Number of counts relates to time required for xcombine_vol_data_vtk
        events = glob(os.path.join(self.cwd, "output", kernel, "*"))
        for event in events:
            event_id = os.path.basename(event)
            for par in self.pars["params"]:
                srcs = glob(os.path.join(event, f"proc*_{par}*.bin"))
                for src in srcs:
                    # e.g. proc000000_vs.bin
                    filename = os.path.basename(src)
                    # Add a unique tag to the filename, different if kernels
                    parts = filename.split("_")
                    # Rejoin end of string if e.g. 'vs_kernel' and not 'vs'
                    parts[1] = "_".join(parts[1:])
                    new_filename = "_".join([parts[0], kernel, event_id,
                                             parts[1]])
                    # Set the source and destination for symlink
                    dst = os.path.join(self.specfem, "SUM", new_filename)
                    os.symlink(os.path.abspath(src), dst)
                    c += 1
        return c

    def run_xcomb(self):
        """
        Use subprocess to call sbatch of an existing script. A bit hacky but
        avoids having to generate an entire sbatch script in here.
        """
        os.chdir(self.specfem)
        stdout = subprocess.run(f"sbatch {self.bash_script}", shell=True,
                                capture_output=True, text=True).stdout

        # Maui specific command to check job status, similar to Seisflows
        for part in stdout.strip().split(" "):
            try:
                # integer ensures no trailing/leading whitespaces
                job_id = str(int(part))
                break
            except ValueError:
                continue
        else:
            sys.exit("subprocess did not return a job number")

        print(f"waiting for job number {job_id} to finish")
        time.sleep(20)
        while True:
            time.sleep(10)
            stdout = subprocess.run(
                            f"sacct -nL -o jobid,state -j {job_id}",
                            shell=True, capture_output=True, text=True).stdout

            # Choose the last occurrence of job id to check status from, as
            # earlier precursor jobs may be completed but main job still running
            # If job is still launching, else will wait another round
            status = None
            for lines in stdout.strip().split("\n"):
                if lines.split()[0] == job_id:
                    status = lines.split()[1]
            if status is None:
                print("waiting for status")
                continue

            if status == "RUNNING":
                print(".", end="")
                continue
            elif status == "PENDING":
                continue
            elif status == "COMPLETED":
                print(f"finished")
                return
            elif status == "FAILED":
                sys.exit(f"failed, exiting")
            else:
                print(f"unknown status {status}")
                continue

    def post_organize(self):
        """
        Remove symlinks and organize files in tagged directories, ensure that the
        directory containing symlinks is empty after organization
        """
        print("organizing output files")

        # Move files
        for i, src in glob(os.path.join(self.specfem, "SUM", "*.vtk")):
            # e.g. model, kernel, gradient
            tag = os.path.basename(src).split("_")[0]
            dst_dir = os.path.join(self.dirout, tag)
            if not os.path.exists(dst_dir):
                os.makedirs(dst_dir)
            dst = os.path.join(dst_dir, os.path.basename(src))
            os.rename(src, dst)
        print(f"\tmoved {i} .vtk files to output directory")

        # Remove symlinks
        c = 0
        for fid in glob(os.path.join(self.specfem, "SUM", "*.bin")):
            if os.path.islink(fid):
                os.remove(fid)
                c += 1
        print(f"\tremoved {c} symlinks from SUM directory")

        # Just make sure the directory is empty, otherwise something went wrong
        assert(len(glob(os.path.join(self.specfem, "SUM", "*"))) == 0)

    def create_srcrcv_vtk(self):
        """
        One-time creation of source and receiver VTK files from the current
        directory using the CMTSOLUTIONS and STATIONS files
        """
        if not os.path.exists(os.path.join(self.dirout, "srcs.vtk")):
            print("creating source vtk file")
            src_vtk_from_specfem(path_to_data=os.path.join(self.cwd, "sfDATA"),
                                 path_out=self.dirout)

        if not os.path.exists(os.path.join(self.dirout, "rcvs.vtk")):
            print("creating receiver vtk file")
            rcv_vtk_from_specfem(path_to_data=os.path.join(self.cwd, "sfDATA"),
                                 path_out=self.dirout)

    def run(self, **kwargs):
        """
        Main run function that calls upon the class to create .vtk from .bin
        """
        self.setup(**kwargs)
        self.discover()
        self.create_srcrcv_vtk()
        for tag in self.available.keys():
            # To be safe, set temp variables to None before they are set again
            self.pars = None

            self.pars = self.available[tag]
            if tag == "kernels":
                for kernel in self.pars["tags"]:
                    if self.organize_kernels(kernel):
                        self.run_xcomb()
                        self.post_organize()
            else:
                if self.organize_models_gradients():
                    self.run_xcomb()
                    self.post_organize()


if __name__ == "__main__":
    xcomb_bash_script = ("/scale_wlg_persistent/filesets/project/nesi00263/"
                         "PyPackages/simutils/seisflows/combvoldatavtk.sh")
    vtko = VTKO(bash_script=xcomb_bash_script)
    vtko.run()




