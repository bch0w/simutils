# NZATOM r0.1 Workflow

Instruction set for how I fixed numerical artefacts in the NZATOM model

## Git Hashes
SPECFEM3D: 5275e1f

## Instructions
1. Generate a) shallow, b) crust, c) mantle models, regular grid. Mesh files
   stored in `simutils/specfem/reinterp_gll_model`
2. `xgenerate_databases` with IRIS EMC proper .xyz files 
3. smooth each model with trial-and-error determined Gaussian
4. `xproject_and_combine_vol_data_on_regular_grid` to generate binary array 
    files. Use `simutils/specfem/project_to_FD_grid/write_proj_grid_fd.py` to 
    create required file `fd_proj_grid.txt`. Grid files used sotred with mesh
    files in (1)
5. `simutils/specfem/project_to_FD_grid/convert_projection_xyz.py` to convert 
    data arrays to .xyz files, ensure that Qmu, Qkappa -> Qp, Qs
6. `simutils/tools/xyz_2_geocsv.py` (simutils) to convert .xyz files to GeoCSV 
    with appropriate headers
8. `iris-edu/emc-tools/GeoCSV_2_netCDF.py` to convert GeoCSV to netCDF files, 
    upload netCDF to IRIS
