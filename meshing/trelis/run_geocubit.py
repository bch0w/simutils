"""
A "simple" script to run Geocubit, replaces a bash script that kept breaking
cause Bash is just a fuckin rock sometimes when I really need a ballpeen hammer.
"""
import os
import numpy as np
import subprocess

# Make sure environment is copied so it doesn't affect current python environ
myenv = os.environ.copy()
mypath = myenv["PATH"]

# Set export paths here
myenv["TRELISHOME"] = "/opt/Trelis-16.1"
myenv["CUBITLIB"] = "/opt/Trelis-16.1/bin:opt/Trelis-16.1/structure"
myenv["CUBITDIR"] = "/opt/Trelis-16.1"
myenv["CUBITHOME"] = "/opt/Trelis-16.1"
myenv["LD_LIBRARY_PATH"] = "{}/bin".format(myenv["CUBITDIR"])
myenv["PATH"] = "{path}:{cubit_dir}/bin".format(path=mypath,
                                                cubit_dir=myenv["CUBITDIR"]
                                                )
myenv["PYTHONPATH"] = ("/home/bchow/.conda/envs/mesher/bin/python:" 
                       "{cub}/bin:{cub}/structure".format(cub=myenv["CUBITDIR"]
                                                          )
                       )

# Set meshing parameters here
fid = "test_nz_north_68_3triples"
base = f"/seis/prj/fwi/bchow/tomo/meshing/trelis/new_zealand/{fid}"
config = f"{base}/{fid}.cfg"
python2 = "/home/bchow/.conda/envs/mesher/bin/python"
geocubit = "/seis/prj/fwi/bchow/packages/GEOCUBIT/GEOCUBIT.py"

# Turn on or off the various parts of the meshing process
run_mesh = True
run_merge = True
run_export = True
make_new_materials_file = True

# Starting workflow here. Get some prerequisite info from the config file
print(f"{fid}")
print("getting number of processors from config")
with open(config, "r") as f:
    lines = f.readlines()
    nproc_xi = nproc_eta = None
    for line in lines:
        if "output_dir" in line:
            output_dir = line.strip().split("=")[1].strip()
        elif "working_dir" in line:
            working_dir = line.strip().split("=")[1].strip()
        elif "number_processor_xi" in line:
            nproc_xi = int(line.strip().split("=")[1])
        elif "number_processor_eta" in line:
            nproc_eta = int(line.strip().split("=")[1])
            break
nproc = nproc_xi * nproc_eta
print(f"xi={nproc_xi}, eta={nproc_eta}, nproc={nproc}")

# Set some necessary directory paths
output_dir = os.path.join(base, output_dir)
export_dir = os.path.join(base, "export_mesh_specfem3d")
if not os.path.exists(export_dir):
    os.makedirs(export_dir)

# Run the mesher based on the number of processors given
if run_mesh:
    print("running geocubit to mesh")
    try:
        # Run these in parallel, redirect stdout to dev null to avoid output
        child_processes = []
        for i in range(nproc):
            if i in [0, nproc-1]:
                print("submitted meshing process {i}/{nproc}".format(i, nproc-1)

            mesh = [python2, geocubit, "--build_volume", "--mesh",
                    "--cfg=" + config, "--id_proc=" + f"{i}"]
            mesh_out = subprocess.Popen(mesh, env=myenv, 
                                        stdout=open(os.devnull, 'w'),
                                        stderr=open(os.devnull, 'w')
                                        )
            child_processes.append(mesh_out)
        
        # Wait for all processes to finish before continuing
        print(f"waiting..., check {working_dir} logs for status")
        for child in child_processes:
            child.wait()

    except CalledProcessError:
        print(mesh_out)
        sys.exit("error running mesh")

# Merge all the NPROC parts of the mesh together
if run_merge:
    print("running geocubit to merge")
    os.chdir(output_dir)
    merge = [python2, geocubit, "--collect", "--merge",
             "--meshfiles=mesh_vol_*.e", "--cpux=" + f"{nproc_xi}",
             "--cpuy=" + f"{nproc_eta}"]
    try:
        merge_out = subprocess.check_call(merge, env=myenv)
    except CalledProcessError:
        print(merge_out)
        sys.exit("error running merge")

# Export the mesh to the necessary files for Specfem3D
if run_export:
    print("running geocubit to export")
    os.chdir(output_dir)
    if os.path.exists(os.path.join(output_dir, "TOTALMESH_MERGED.e")):
        export = [python2, geocubit, "--export2SPECFEM3D",
                  "--meshfiles=TOTALMESH_MERGED.e"
                  ]
        try:
            export_out = subprocess.check_call(export, env=myenv)
        except CalledProcessError:
            print(export_out)
            sys.exit("error running export")
    else:
        sys.exit("TOTALMESH_MERGED.e does not exist")

    # Move the files to a single directory for easy export
    for export_qty in ["mesh_file", "materials_file", "nodes_coords_file",
                       "free_or_absorbing_surface_file_zmax",
                       "absorbing_surface_file_bottom",
                       "absorbing_surface_file_xmin",
                       "absorbing_surface_file_ymin",
                       "absorbing_surface_file_xmax",
                       "absorbing_surface_file_ymax",
                       "nummaterial_velocity_file"]:
        src = os.path.join(output_dir, export_qty)
        if os.path.exists(src):
            dst = os.path.join(export_dir, export_qty)
            os.rename(src, dst)

if make_new_materials_file:
    def make_new_materials_file(export_dir, layer=-8E3):
        """
        Make new materials file, adapted from Carl Tape's Matlab script
        External tomography files require different sets of material ids
        They get assigned here and shoved into "materials_file"
        """
        if os.path.exists(os.path.join(export_dir, "materials_file_original")):
            print("this function has already been run")
            return

        # Read in the created mesh files to get some info, skip header
        nodes_coords = np.loadtxt(
            os.path.join(export_dir, "nodes_coords_file"), skiprows=1
        )
        mesh_file = np.loadtxt(os.path.join(export_dir, "mesh_file"),
                               skiprows=1
                               )
        materials_file = np.loadtxt(os.path.join(export_dir, "materials_file"))

        # determine the total list of mesh elements that require ids
        n_z = nodes_coords[:, 3]  # element depth
        material_ids = materials_file[:, 0]
        material_indices = materials_file[:, 1]

        # delete the first column of the mesh file to get the indices for
        # each element
        node_indices = np.delete(mesh_file, 0, axis=1)

        # Rename the old materials file so we can take the name
        os.rename(src=os.path.join(export_dir, "materials_file"),
                  dst=os.path.join(export_dir, "materials_file_original")
                  )
        # Write to new file
        new_material_indices = []
        for i, ids in enumerate(material_ids):
            # Get the z value of all the nodes that define the element
            inodes = node_indices[i]
            inodes_as_indices = (inodes - 1).astype(int)
            zvals = n_z[inodes_as_indices]  # need int to use as indices
            min_zval = np.min(zvals)  # find the deepest node on the element

            mat_id = int(material_ids[i]) - 1
            if material_indices[mat_id] == 1:
                new_material_indices.append(-3)
            else:
                # Defines the shallow, separated at user parameter 'layer'
                if min_zval >= layer:
                    new_material_indices.append(-1)
                else:
                    new_material_indices.append(-2)

        # Join the two column-wise
        new_materials_file = np.column_stack((material_ids,
                                              np.array(new_material_indices))
                                             )
        # Write to file
        np.savetxt(os.path.join(export_dir, "materials_file"),
                   new_materials_file, "%d")

    print("making new materials file")
    make_new_materials_file(export_dir)







