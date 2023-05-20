# Project to FD Grid

Attempting to use SPECFEM3D\_Cartesian function 
`xproject_and_combine_vol_data_on_regular_grid` to interpolate a GLL mesh onto
a regular grid (FD = finite difference).

Function must be called in a valid SPECFEM3D working directory with available
database files (requires `external_mesh` files).

Call looks like

```bash
$ mpirun -n $NPROC bin/xproject_and_combine_vol_data_on_regular_grid \
    $DATA_FILENAME $INPUT_DIRECTORY $OUTPUT_DIRECTORY
```

Where $DATA\_FILENAME should match one of your model parameters (e.g., 'vs'), 
and $INPUT\_DIRECTORY should point to your DATABASES\_MPI/local path directory.

Important source code files to look through to see what's happening:
- auxiliaries/project\_and\_combine\_vol\_data\_on\_regular\_grid.f90
- inverse\_problem\_for\_model/projection\_on\_FD\_grid\_mod.f90

Requires a file `fd_proj_grid.txt` which defines the finite difference grid
you are trying to project onto.

Format of that file is:

```
ox oy oz 
hx hy hz
nx ny nz
```

- o? defines the origin point (should match your GLL points)
- h? defines the sampling rate 
- n? defines the number of points in each direction

The final file will be a single array (FORTRAN binary) with npoint nx * ny * nz



