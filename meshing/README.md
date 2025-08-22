# SPECFEM3D MESHING
The following utilities are used for meshing cartesian meshes for SPECFEM3D 
Cartesian. There are utilities to generate meshfem3D files, topography files
for your top interface, and to apply stochastic perturbations to an external
tomography model (.xyz). See each of the main functions for detailed 
instructions, but this README provides some basic intro, because inevitebly I 
will forget what I did

## Tool Overview

### `meshfem`
Contains utilities for generating mesh files compatible with SPECFEM3D. These scripts help define the geometry and discretization of your simulation domain.

**Example usage:**
```bash
python prepare_meshfem.py configs/<FID>.json
```

Meshfem3D files will be stored in `created/<FID>`

---

### `stochperts`
Provides tools to apply stochastic perturbations to a chosen 1D model. Future work will have this applied directly to a given .xyz model.

**Example usage:**
```bash
python stochpert_model.py configs/<FID>.json
```

.xyz tomography files will be stored in `created/<FID>`

---

### `topography`
Includes scripts for generating topography files from SRTM30P for the top interface of your mesh.

**Example usage:**
```bash
python generate_topo.py configs/<FID>.json
```

If you point the config file to the correct output of `meshfem`, it will automatically update the required files, otherwise you will need
to do some manual re-tooling and file movement.

---

Refer to each subdirectory's README or main script for more details and options.
