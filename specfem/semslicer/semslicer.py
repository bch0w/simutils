"""
API to use the SPECFEM3D Fortran code sem_model_slice.f90. Provides utility
functions to generate regular grids, run the script on a SLURM HPC system,
and generate the input external tomography (.xyz) files required for a
subsequent SPECFEM run.

.. note:: The fortran binary file should be named 'xsemslicer'

.. note:: 
    Here are a few notes related to running semslicer from problems that I
    encountered while running it:

    * Make sure the spacing of your structured grid is smaller than the smallest
    gll point size in each direction. This ensures that you are actually 
    sampling each gll point. If you don't do this, the nearest neighbor 
    selection is at risk of pulling values that are quite far from expected 
    depending on how strong the parameter gradient is. 
    You made need to downsample afterwards
    * Beware of topography, both actual topography, and mesh topography
    introduced by doubling layers and high-skew elements. The nearest neighbor
    algorithm has trouble grabbing the correct values for regions with high skew
    which may lead to large differences between expected and interpolated vals.
    For me this was evident in a slice across the entire domain showing uniform
    difference between expected and interpolated values. Sampling my structured
    grid at very high sampling (< smallest GLL distance) remedied this
    * Breaking your model into regions (e.g. shallow, crust, mantle) is a good
    idea as usually the required grid spacing gets larger with depth. Even then
    very high sampling rate in the shallow may encounter segfaults due to the 
    extremely large number of points (10s of millions for me). Here you may need
    to break that region into multiple sub-regions (shallow_1, shallow_2, etc..)
    and then recombine them into a single region file.
    * When recombining region files, ensure that the point order is sorted in 
    the correct order, as SPECFEM (and semslicer) are stupid and dont actually
    read the grid point values, but assume it based on the start and end points.
    That means, grid points must go in order Z, Y, X, which looks like this:

    min_x,       min_y          min_z
    min_x + dx   min_y          min_z
    ...
    max_x        min_y          min_z
    min_x        min_y + dy     min_z
    ...
    min_x        max_y          min_z
    min_x        min_y          min_z + dz
    ...

    Hopefully you can figure out the rest from there. Took me three straight 
    weeks to figure out all this stuff, and what a struggle it was...
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
            if isinstance(value, float) or isinstance(value, int):
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
            try:
                setattr(self, key, Pot(config[key]))    
            except ValueError:
                setattr(self, key, config[key])

        # Convert values from km to m if required
        assert(self.domain.units.lower() in ["m", "km"]), \
            "units must be either 'm' or 'km'"
        if self.domain.units == "km":
            for part in vars(self).values():
                if isinstance(part, Pot):
                    part.mult(by=1E3)

        # Will be used to dynamically loop through all regions
        nreg = 1
        for key in vars(self).keys():
            if "region_" in key:
                nreg += 1
        self.nregions = nreg

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
        subregions = {}
        for i in range(1, self.nregions):
            reg = getattr(self, f"region_{i}")
            z_reg = np.arange(reg.zmin, reg.zmax + reg.dz, reg.dz)

            try:
                # Check if this region is split into subregions, these should 
                # NOT overlap
                int(reg.tag.split("_")[-1])
                continue
            except ValueError:
                # ValueError thrown in normal conditions, i.e. no subregion
                z_mins.append(z_reg.min())
                z_maxs.append(z_reg.max())
        
        z_mins = np.array(z_mins)
        z_maxs = np.array(z_maxs)       
 
        # Check that there are no gaps in the depths covered by each region
        # assert ((z_maxs[1:] - z_mins[:-1] >= 0)).all(), f"Gaps in depth ranges"

    def make_grids(self, write=True):
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
            if write:
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

    def submit_jobs_maui(self, model_dir):
        """
        Submit an SBATCH job to the NZ cluster Maui. Auto generate the output
        file name based on the input grid file and data name.
        """
        # Determine parameters from model directory and make the grid files
        assert(os.path.exists(model_dir)), f"{model_dir} does not exist"
        self.make_grids()
        if self.nproc is None:
            self.get_model_values(model_dir)

        # Scratch sbatch file to be overwritten multiple times
        script_fid = os.path.join(self._scratch, f"run_semslicer.sh")

        # Submit job for each region and data file
        for data_name in self.data:
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

    def check_xsemslicer_format(self):
        """
        Check output files for Fortran values that replace double entries 
        that are filled in with a 2* value with the actual values.
       
        .. note:: 
            Should be run by other functions, assumed that self.make_grids() has
            already been run

        e.g. A column that is meant to be 0., 0., 1., would be written as
             2*0., 1., ... We want to replace the 2* with two zeros.
        """
        self.make_grids(write=False)

        line_fmt = "{0:f}, {1:f}, {2:f}, {3:f}, {4:f}\n"

        for i in range(1, self.nregions):
            region = getattr(self, f"region_{i}")
            for j, data_name in enumerate(self.data):
                edited_lines = 0

                fid = f"{data_name}_{region.grid_tag}"
                print(fid)
                with open(fid, "r") as f:
                    lines = f.readlines()
                 
                for k, line in enumerate(lines[:]):
                    if "*" in line:
                        edited_lines += 1
                        new_array = []
                        for value in line.strip().split(","):
                            if "*" in value:
                                instances, value = value.split("*")
                                for _ in range(int(instances)):
                                    new_array.append(float(value))
                            else:
                                new_array.append(float(value))

                        new_array = line_fmt.format(*new_array)
                        print(f"\t{line}\n\t{new_array}")
                        lines[k] = new_array

                if edited_lines:
                    print(f"{fid} has {edited_lines} edited lines")
                    with open(fid, "w") as f:
                        f.writelines(lines)

    def write_xyz_files(self):
        """
        Write the header and the XYZ files based on the output values of the
        fortran file. Filenames based on 'model_tomography_bryant.f90'
        """
        # Get the grid attributes incase this is run in a separate instance
        self.make_grids(write=False)

        for i in range(1, self.nregions):
            # These are the required values for the external tomo files
            data = Pot({_:[] for _ in self.data})

            region = getattr(self, f"region_{i}")

            # Collect data from the output files
            for j, data_name in enumerate(data.keys()):
                fid = f"{data_name}_{region.grid_tag}"
                assert(os.path.exists(fid)), f"{fid} does not exist"

                # We only need x, y, and z once
                try:
                    if j == 0:
                        x, y, z, v, d = np.loadtxt(fid, dtype=float,
                                                   delimiter=",").T

                        # Simple check to see what the min and max distance from 
                        # grid point to nearest neighbor gll point. Useful for 
                        # catching unrealistically large values of d
                        print(f"{region.grid_tag}:\n"
                              f"MINIMUM DISTANCE = {d.min()}m\n"
                              f"MAXIMUM DISTANCE = {d.max()}m\n"
                              )

                    else:
                        _, _, _, v, _ = np.loadtxt(fid, dtype=float,
                                                   delimiter=",").T 
                except ValueError as e:
                    raise ValueError(f"{fid} needs manual review") from e

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

    def combine_subregions(self, sort=True):
        """
        If a region has very dense grid spacing, the cluster may segfault due to
        the temporary storage of a very large number of points. To work around 
        this, a region can be split into sub-regions with the same grid spacing,
        but covering a range of depths, effectively splitting the job using 
        embarassing-parallelization. This utility function recombines the 
        necessary subregions so that the function write_xyz_files() can be 
        called with the expected behavior.

        .. warning::
            This only handles the case where the zmin/zmax values of adjacent
            layers are touching, or not. If the layers overlap, then redundant
            points will be written into the resulting .xyz file.
            To avoid errors, ensure that your subregions touch, or are separated
            by 'dz'

        :type sort: bool
        :param sort: sort the output files by coordinate. required for SPECFEM
            to be able to work with the resulting .xyz files. But if the files
            are quite big then memory errors become a problem and this may have
            to be done externally
        """
        self.make_grids(write=False) 

        # This defines the line format passed to np.savetxt
        xyz_fmt = ["%8.1f", "%10.1f", "%8.1f", "%12.5f", "%12.5f"] 

        # Determine which regions are divided into sub-regions
        region_lists = {}
        for i in range(1, self.nregions):
            region_label = f"region_{i}"
            region = getattr(self, region_label)
            region_name = region.tag.split("_")[0]
            if region_name not in region_lists:
                region_lists[region_name] = [region_label]         
            else:
                region_lists[region_name].append(region_label)

        # For any regions divided into sub-regions, combine the .xyz files
        # into a single region file for each data parameter
        for data_name in self.data:
            for region_name, region_list in region_lists.items():
                # Don't care about the regions that aren't sub-divided
                if len(region_list) <= 1:
                    continue

                # Continuously write to disk to avoid storing data 
                fid_out = f"{data_name}_{region_name}_grid.xyz"
                print(fid_out)

                # Since we're opening the file in append mode, ensure that some
                # old file doesn't already exist
                if os.path.exists(fid_out):
                    os.remove(fid_out)

                with open(fid_out, "ab") as f:
                    # Used to check if the Z values overlap, which we don't want 
                    # because that is redundant
                    zmin_last, zmax_last = None, None

                    for i, region_label in enumerate(region_list):
                        region = getattr(self, region_label)
                        print(f"\t{region.tag}")

                        # Can only start checking depth after the first subregn
                        zmin = region.zmin
                        zmax = region.zmax
                        if i >= 1:
                            # Adjacent boundaries means have to exclude 
                            # redundant depth slices. Works if regions are being
                            # built upwards or downwards
                            if zmax == zmin_last:
                                zmax -= region.dz
                            elif zmin == zmax_last:
                                zmin += region.dz

                        fid = f"{data_name}_{region.grid_tag}"
                        d = np.loadtxt(fid, dtype=float, delimiter=',')

                        # Exclude adjacent boundaries, here we're assuming the 
                        # structure of the xyz file, i.e. Z values in 2nd column
                        _len_d_orgn = len(d)
                        d = d[np.where((d[:,2] >= zmin) & (d[:,2] <= zmax))]
                        print(f"\tboundaries npts: {_len_d_orgn} -> {len(d)}")
                        
                        np.savetxt(f, d, delimiter=",", fmt=xyz_fmt)
                                    
                        zmin_last = zmin
                        zmax_last = zmax
                
                # Free up memory before moving onto sorting
                if sort:
                    del d

                    # Sorta hacky: re-read written file and sort depth values in
                    # descending order since they won't have any natural order 
                    # and the xyz values need to match when writing into single 
                    # xyz file which means reverse sorted by Z then Y then X
                    d = np.genfromtxt(fid_out, delimiter=",", dtype=None,
                                      names=["x", "y", "z", "v", "d"])
                    d.sort(order=["z", "y", "x"])

                    print(f"\t{fid_out}: {len(d)} pts")
                    np.savetxt(fid_out, d, delimiter=",", fmt=xyz_fmt)
                

if __name__ == "__main__":
    cfg_file = "cfg.yaml"
    model = "birch_m11"
    ss = Semslicer(cfg_file)

    if len(sys.argv) > 1:
        if sys.argv[1] == "submit":
            # Step 1: submit xsemslicer jobs to Maui to interrogate 
            ss.submit_jobs_maui(model)
        elif sys.argv[1] == "check":
            # Step 1a: check that the xsemslicer outputs are formatted properly,
            # not mandatory but 'combine' and 'write' may fail if not run
            ss.check_xsemslicer_format()
        elif sys.argv[1] == "combine":
            # Step 2 (Optional): if region was split into segments, combine segs
            ss.combine_subregions(sort=False)
        elif sys.argv[1] == "write":
            # Step 3: write individual parameters into SPECFEM tomography files
            ss.write_xyz_files()
    else:
        # Define the grid structure that will be used for interrogation
        ss.make_grids(write=False)
           
