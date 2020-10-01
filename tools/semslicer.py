"""
API to use the SPECFEM3D Fortran code sem_model_slice.f90. Provides utility
functions to generate regular grids, run the script on a SLURM HPC system,
and generate the input external tomography (.xyz) files required for a
subsequent SPECFEM run.

NOTE: The fortran binary file should be named 'semslicer'
"""
import os
import yaml
import datetime
import subprocess
import numpy as np
from glob import glob


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

    def clean(self):
        """
        Delete scratch directory
        """
        if os.path.exists(self._scratch):
            os.remove(self._scratch)
            os.makedirs(self._scratch)

    def make_grids(self):
        """
        Generate the input grid files that will be used for querying the model
        Store the number of points in the internal attributes. Optional write
        grid to text files.
        """
        for i in range(1, self.nregions):
            reg = getattr(self, f"region_{i}")

            # Establish the regular grid based on the particular region
            x_reg = np.arange(self.domain.xmin, self.domain.xmax, reg.dx)
            y_reg = np.arange(self.domain.ymin, self.domain.ymax, reg.dy)
            z_reg = np.arange(reg.zmin, reg.zmax, reg.dz)

            x_grid, y_grid = np.meshgrid(x_reg, y_reg)
            x_out = x_grid.flatten()
            y_out = y_grid.flatten()

            npts = len(x_reg) * len(y_reg) * len(z_reg)
            print(f"{reg.tag}: {npts} points")
            grid_tag = f"{reg.tag}_grid.xyz"

            reg.npts = npts
            reg.grid_tag = grid_tag
            # Add tag and npts attribute to the region
            setattr(self, f"region_{i}", reg)

            # Save these grid files to the scratch directory
            if not os.path.exists(os.path.join(self._scratch, grid_tag)):
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
        self.nproc = max([int(_.split('_')[0][4:]) for _ in binfiles])

        # list comprehension strips all but var e.g. 'proc000001_vp.bin' -> 'vp'
        variables = [os.path.splitext(_.split('_')[1])[0] for _ in binfiles]
        self.variables = list(set(variables))

    def submit_jobs_maui(self, model_dir, data_names=None):
        """
        Submit an SBATCH job to the NZ cluster Maui. Auto generate the output
        file name based on the input grid file and data name.
        :return:
        """
        if data_names is None:
            data_names = ["vp", "vs", "rho", "qmu", "qkappa"]

        # Determine parameters from model directory and make the grid files
        assert(os.path.exists(model_dir)), f"{model_dir} does not exist"
        self.make_grids()
        if self.nproc is None:
            self.calc_nproc(model_dir)

        # Scratch name
        script_fid = os.path.join(self.scratch, f"run_semslicer.sh")

        # Submit job for each region and data file
        for data_name in data_names:
            assert(data_name in self.variables), \
                f"{data_name} not in {self.variables}"
            for i in range(1, self.nregions):
                region = getattr(self, "region_{i}")
                out_file = f"{data_name}_{region.grid_tag}"

                # Overwrite the run script each time
                with open(script_fid, "w") as f:
                    f.write("#!/bin/bash -e\n\n"
                            "srun -n {nproc:d} semslicer ")
                    f.write(f"{os.path.join(self._scratch, region.grid_tag)} ")
                    f.write(f"{os.path.join(model_dir, '')} ")
                    f.write(f"{data_name} ")
                    f.write(f"{out_file} ")

                # Submit the run script that was just generated
                jobtime = str(datetime.timedelta(minutes=region.npts / 25000))

                submit_call = " ".join([
                    "sbatch",
                    "--account=nesi00263",
                    "--job-name=semslicer",
                    f"--nodes={region.nproc // 40}",  # Maui has 40 cores/node
                    f"--ntasks={region.nproc}",
                    "--cpus-per-task=1",
                    "--clusters=maui",
                    "--partition=nesi_research",
                    f"--time={jobtime}",
                    "--output=semslicer_%j.out"
                    "script_fid"
                ])

                # Submit jobs with no job checking, since this is only run once
                # in a while, we can just manual check the outputs
                subprocess.run([submit_call], capture_output=True)

    def write_xyz_files(self):
        """
        Write the header and the XYZ files based on the output values of the
        fortran file
        :return:
        """
        # Get the grid attributes incase this is run in a separate instance
        self.make_grids()

        # These are the required values for the external tomo files
        data = Pot(vp=[], vs=[], rho=[], qmu=[], qkappa=[])

        for i in range(1, self.nregions):
            region = getattr(self, "region_{i}")

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

            dx = x[1] - x[0]
            dy = y[1] - y[0]
            dz = z[1] - z[0]

            # Check domain values against internal attributes
            assert (x.min() == self.domain.xmin), "xmin does not match"
            assert (x.max() == self.domain.xmax), "xmax does not match"
            assert (y.min() == self.domain.ymin), "ymin does not match"
            assert (y.min() == self.domain.ymax), "ymax does not match"

            # Check region values against internal attributes
            assert (z.min() == region.zmin()), \
                f"{region.tag} zmin does not match"
            assert (z.max() == region.zmax()), \
                f"{region.tag} zmax does not match"
            assert (dx == region.dx), f"dx does not match"
            assert (dy == region.dy), f"dy does not match"
            assert (dz == region.dz), f"dz does not match"

            # Use the output values to define the header
            with open(f"tomography_model_{region.tag}.xyz", "w") as f:
                # Line 1, min and max values
                f.write(f"{x.min():.1f} {y.min():.1f} {z.max():.1f} "
                        f"{x.max():.1f} {y.max():.1f} {z.min():.1f}\n")

                # Line 2, spacing values. Assuming spacing is constant
                f.write(f"{dx:.1f} {dy:.1f} {dz:.1f}\n")

                # Line 3, npts for each direction
                f.write(f"{len(x):d} {len(y):d} {len(z):d}\n")

                # Line 4, min and max values for each quantity
                f.write(f"{data.vp.min():.1f} {data.vp.max():.1f} ")
                f.write(f"{data.vs.min():.1f} {data.vs.max():.1f} ")
                f.write(f"{data.rho.min():.1f} {data.rho.max():.1f}\n")

                # Lines 5+: x, y, z, vp, vs, rho, qp, qs
                for k in range(len(x)):
                    f.write(f"{x[k]:.1f} {y[k]:.1f} {z[k]:.1f} ")
                    f.write(f"{data.vp[k]:.1f} {data.vs[k]:.1f} ")
                    f.write(f"{data.rho[k]:.1f} {data.qmu[k]:.1f} ")
                    f.write(f"{data.qkappa[k]:.1f}\n")




