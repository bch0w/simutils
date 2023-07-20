# Project to FD Grid

SPECFEM3D\_Cartesian function `xproject_and_combine_vol_data_on_regular_grid` 
interpolates a GLL model onto a regular grid (FD = finite difference) to convert
an unstructured mesh to a regular XYZ grid.

Function must be called in a valid SPECFEM3D working directory with available
database files (requires `external_mesh` files).

## Workflow

General steps to convert your GLL model (.bin files) to an external tomography
model (.xyz files)

### 1. Make `xproject_and_combine_vol_data_on_regular_grid`

Important source code files to look through to see what's happening:
- auxiliaries/project\_and\_combine\_vol\_data\_on\_regular\_grid.f90
- inverse\_problem\_for\_model/projection\_on\_FD\_grid\_mod.f90

Configure specfem3D as normal and run the following:

```bash
cd specfem3D/
make project_and_combine_vol_data_on_regular_grid
```

### 2. Write `fd_proj_grid.txt`

This exectuable requires a file `fd_proj_grid.txt` in the main working 
directory which defines the finite difference grid you are trying to project 
your model on to.

Use the script `write_proj_grid_fd.py` to define your finite difference grid
and create the required grid file `fd_proj_grid.txt`. You must fill in the 
dictionary `fd_grid` in this script based on your mesh dimensions.

> __Note__: You may want to split your mesh up into multiple chunks 
    (shallow, crust, mantle) for convenience, and potentially due to memory 
    errors. See warning below for explanation.

The format of the output file `fd_proj_grid.txt` is:

```
ox oy oz 
hx hy hz
nx ny nz
```
Where:
- o? defines the origin point (should match your GLL points)
- h? defines the sampling rate 
- n? defines the number of points in each direction

> __Warning__: For large meshes, it is possible that you will run into a 
    SEGFAULT issue for very-finely sample FD grids. You may need to split your 
    FD grid into smaller chunks,  or decrease the sampling rate / number of 
    grid points. I did not find a suitable workaround for this issue.

### 3. Run `xproject_and_combine_vol_data_on_regular_grid`

You will need to run `xproject_and_combine_vol_data_on_regular_grid` for **EACH**
parameter you have in your velocity model (e.g., vp, vs, rho, qp, qs).

The system call looks like:

```bash
$ mpirun -n $NPROC bin/xproject_and_combine_vol_data_on_regular_grid \
    $DATA_FILENAME $INPUT_DIRECTORY $OUTPUT_DIRECTORY
```

Where $DATA\_FILENAME should match one of your model parameters (e.g., 'vs'), 
and $INPUT\_DIRECTORY should point to your DATABASES\_MPI/ (local path) directory.

The output file will be a single array (FORTRAN binary) called 
`*_projected.bin`, the length of the array will be npoints = nx * ny * nz.

### 4. Convert FD projection to XYZ

If you were able to run 1--3 successfully, you should have a number of FORTRAN
binary files (*_projected.bin). You can then use the script 
`convert_projection_xyz.py` to convert these .bin files to a .xyz file.

The Python script requires the file `fd_proj_grid.txt` and a path to your 
created binary files. It will read these files and create an XYZ tomography
file out of them. It also takes care of converting Qmu and Qkappa to Qp and Qs,
which is the required input format for SPECFEM3D.

Call this simply with
```python
python convert_projection_xyz.py
```

