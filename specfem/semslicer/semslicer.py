"""
API to use the SPECFEM3D Fortran code sem_model_slice.f90. Provides utility
functions to generate regular grids, run the script on a SLURM HPC system,
and generate the input external tomography (.xyz) files required for a
subsequent SPECFEM run.

NOTE: The fortran binary file should be named 'xsemslicer'
"""
import os
import sys
import yaml
import datetime
import subprocess
import numpy as np
from glob import glob
from shutil import rmtree


class Pot(dict):
    """
    A pot to hold things: dictionary that allows access via attributes
    """
    def __setattr__(self, key, value):
        self[key] = value

    def __getattr__(self, key):
        return self[key]

    def mult(self, by):
        """
        Multiply all floats held inside the Pot by value 'by'
        :type by: float
        :param by: multiply by this value
        """
        for key, value in self.items():
            if isinstance(value, float):
                setattr(self, key, value * by)


class Semslicer:
    """
    A wrapper class for the model slice capabilities of Specfem, used to
    interpolate the unstructured Specfem grid onto a regular XYZ coordinate
    system which can then be used as a new external tomography file, or for
    plotting capabilities.

    .. note::
        All internal values are to be in units of meters
    """
    def __init__(self, cfg_fid):
        """
        Semslicer only requires an input configuration file, values will be
        stored as internal attributes
        """
        assert(os.path.exists(cfg_fid)), f"{cfg_fid} does not exist"
        config = yaml.safe_load(open(cfg_fid))
        for key in config:
            setattr(self, key, Pot(config[key]))

        # Convert values from km to m if required
        assert(self.domain.units.lower() in ["m", "km"]), \
            "units must be either 'm' or 'km'"
        if self.domain.units == "km":
            for part in vars(self).values():
                if isinstance(part, Pot):
                    part.mult(by=1E3)

        # Will be used to dynamically loop through all regions
        self.nregions = len(vars(self))

        self.nproc = None
        self.variables = None
        self._scratch = "./scratch"

        if not os.path.exists(self._scratch):
            os.makedirs(self._scratch)

    def clean(self):
        """
        Delete scratch directory
        """
        if os.path.exists(self._scratch):
            rmtree(self._scratch)
            os.makedirs(self._scratch)

    def check_depth_ranges(self):
        """
        Using `arange` to create the structured grids may sometimes leave gaps
        in the depth ranges if the start/end points and spacing are not chosen
        correctly. Quickly check that the regions cover the entire depth range
        """ 
        z_mins, z_maxs = [], []
        for i in range(1, self.nregions):
            reg = getattr(self, f"region_{i}")
            z_reg = np.arange(reg.zmin, reg.zmax + reg.dz, reg.dz)
            z_mins.append(z_reg.min())
            z_maxs.append(z_reg.max())
        
        z_mins = np.array(z_mins)
        z_maxs = np.array(z_maxs)       
 
        # Check that there are no gaps in the depths covered by each region
        assert ((z_maxs[1:] - z_mins[:-1] >= 0)).all(), f"Gaps in depth ranges"

    def make_grids(self):
        """
        Generate the input grid files that will be used for querying the model
        Store the number of points in the internal attributes. Optional write
        grid to text files.
        """
        self.check_depth_ranges()

        for i in range(1, self.nregions):
            reg = getattr(self, f"region_{i}")

            # Establish the regular grid (inclusive) based on each region
            x_reg = np.arange(self.domain.xmin, self.domain.xmax + reg.dx, 
                              reg.dx)
            y_reg = np.arange(self.domain.ymin, self.domain.ymax + reg.dy, 
                              reg.dy)
            z_reg = np.arange(reg.zmin, reg.zmax + reg.dz, reg.dz)

            x_grid, y_grid = np.meshgrid(x_reg, y_reg)
            x_out = x_grid.flatten()
            y_out = y_grid.flatten()

            npts = len(x_reg) * len(y_reg) * len(z_reg)
            print(f"{reg.tag:<8}: {npts} points")
            grid_tag = f"{reg.tag}_grid.xyz"

            reg.npts = npts
            reg.grid_tag = grid_tag
            # Add tag and npts attribute to the region
            setattr(self, f"region_{i}", reg)

            # Save these grid files to the scratch directory
            with open(os.path.join(self._scratch, grid_tag), "w") as f:
                for z in z_reg:
                    for x, y in zip(x_out, y_out):
                        f.write(f"{x:16.3f}\t{y:16.3f}\t{z:16.3f}\n")

    def get_model_values(self, model_dir):
        """
        Get the number of processors to be used in the job submission
        """
        binfiles = [os.path.basename(_) for _ in
                    glob(os.path.join(model_dir, "*"))]

        # list comprehension strips string parts, e.g. 'proc000001_vp.bin' -> 1
        self.nproc = max([int(_.split('_')[0][4:]) for _ in binfiles]) + 1

        # list comprehension strips all but var e.g. 'proc000001_vp.bin' -> 'vp'
        variables = [os.path.splitext(_.split('_')[1])[0] for _ in binfiles]
        self.variables = list(set(variables))

    def submit_jobs_maui(self, model_dir, data_names=None):
        """
        Submit an SBATCH job to the NZ cluster Maui. Auto generate the output
        file name based on the input grid file and data name.
        """
        if data_names is None:
            # These are the required values for an external tomography file
            data_names = ["vp", "vs", "rho", "qmu", "qkappa"]

        # Determine parameters from model directory and make the grid files
        assert(os.path.exists(model_dir)), f"{model_dir} does not exist"
        self.make_grids()
        if self.nproc is None:
            self.get_model_values(model_dir)

        # Scratch sbatch file to be overwritten multiple times
        script_fid = os.path.join(self._scratch, f"run_semslicer.sh")

        # Submit job for each region and data file
        for data_name in data_names:
            assert(data_name in self.variables), \
                f"{data_name} not in {self.variables}"
            for i in range(1, self.nregions):
                region = getattr(self, f"region_{i}")
                out_file = f"{data_name}_{region.grid_tag}"
                log_file = os.path.join(self._scratch, "log_semslicer_%j.out")

                # Job time is empirically derived with some conservative buffer
                job_time = str(datetime.timedelta(minutes=region.npts / 20000))

                # Strip microseconds from time
                job_time = job_time.split(".")[0]

                # Overwrite the run script each time, kinda nasty because 
                # were relative pathing back and forth. Maybe should fix this?
                with open(script_fid, "w") as f:
                    f.write(f"#!/bin/bash -e\n\n"
                            f"#SBATCH --account=nesi00263\n"
                            f"#SBATCH --job-name=semslicer\n"
                            f"#SBATCH --nodes={self.nproc // 40}\n"
                            f"#SBATCH --ntasks={self.nproc}\n"
                            f"#SBATCH --cpus-per-task=1\n"
                            f"#SBATCH --clusters=maui\n"
                            f"#SBATCH --partition=nesi_research\n"
                            f"#SBATCH --time={job_time}\n"
                            f"#SBATCH --output={log_file}\n"
                            )
                    f.write("\n")
                    f.write(f"srun -n {self.nproc:d} xsemslicer "
                            f"{os.path.join(self._scratch, region.grid_tag)} "
                            f"{os.path.join(model_dir, '')} "
                            f"{data_name} "
                            f"{out_file} "
                            )

                # Submit jobs with no job checking, since this is only run once
                # in a while, we can just manual check the outputs
                subprocess.run(["sbatch", script_fid])

    def write_xyz_files(self):
        """
        Write the header and the XYZ files based on the output values of the
        fortran file. Filenames based on 'model_tomography_bryant.f90'

        ..warning::
            Sometimes if adjacent values are the same, Fortran will output them
            formatted together: e.g. 0 0 0 => 3*0
            This will cause a ValueError in numpy.loadtxt(). Could write a 
            function to clean the files beforehand, but happens so rarely that
            manually changing may be enough.
        """
        # Get the grid attributes incase this is run in a separate instance
        self.make_grids()

        for i in range(1, self.nregions):
            # These are the required values for the external tomo files
            data = Pot(vp=[], vs=[], rho=[], qmu=[], qkappa=[])

            region = getattr(self, f"region_{i}")

            # Collect data from the output files
            for j, data_name in enumerate(data.keys()):
                fid = f"{data_name}_{region.grid_tag}"
                assert(os.path.exists(fid)), f"{fid} does not exist"

                # We only need x, y, and z once
                if j == 0:
                    x, y, z, v, _ = np.loadtxt(fid, dtype=float,
                                               delimiter=",").T
                else:
                    _, _, _, v, _ = np.loadtxt(fid, dtype=float,
                                               delimiter=",").T

                data[data_name] = v

            # Get the actualy grid point and spacing values
            xvals = np.unique(x)
            yvals = np.unique(y)
            zvals = np.unique(z)

            dx = xvals[1] - xvals[0]
            dy = yvals[1] - yvals[0]
            dz = zvals[1] - zvals[0]

            # External tomo file requires qp and qs but outputs qk and qm
            # Qp must be derived from Vp, Vs, Qmu, Qkappa;  Qs == Qmu
            # Conversion following Dahlen and Tromp Eq. 9.59
            f = 4 / 3 * (data.vs / data.vp) ** 2
            data["qp"] = 1 / ((1 - f) / data.qkappa + f / data.qmu) 

            # Use the output values to define the header, they may be different
            # from the input cfg values that are stored internally
            with open(f"tomography_model_{region.tag}.xyz", "w") as f:
                # Line 1, min and max spatial values
                f.write(f"{x.min():.1f} {y.min():.1f} {z.min():.1f} "
                        f"{x.max():.1f} {y.max():.1f} {z.max():.1f}\n")

                # Line 2, spacing values, assuming spacing is constant
                f.write(f"{dx:.1f} {dy:.1f} {dz:.1f}\n")

                # Line 3, number of unique points for each direction
                f.write(f"{len(xvals):d} {len(yvals):d} {len(zvals):d}\n")

                # Line 4, min and max values for each quantity
                f.write(f"{data.vp.min():.1f} {data.vp.max():.1f} ")
                f.write(f"{data.vs.min():.1f} {data.vs.max():.1f} ")
                f.write(f"{data.rho.min():.1f} {data.rho.max():.1f}\n")

                # Lines 5+: x, y, z, vp, vs, rho, qp, qs
                for k in range(len(x)):
                    f.write(f"{x[k]:.1f} {y[k]:.1f} {z[k]:.1f} ")
                    f.write(f"{data.vp[k]:.1f} {data.vs[k]:.1f} ")
                    f.write(f"{data.rho[k]:.1f} {data.qp[k]:.1f} ")
                    f.write(f"{data.qmu[k]:.1f}\n")


if __name__ == "__main__":
    assert(len(sys.argv) > 1), "Argument must be 'submit' or 'write'"
    
    cfg_file = "cfg.yaml"
    model = "model_notopo"

    if sys.argv[1] == "submit":
        ss = Semslicer(cfg_file)
        ss.submit_jobs_maui(model)
    elif sys.argv[1] == "write":
        ss = Semslicer(cfg_file)
        ss.write_xyz_files()
    else:
        "Argument must be 'submit' or 'write'"
           
