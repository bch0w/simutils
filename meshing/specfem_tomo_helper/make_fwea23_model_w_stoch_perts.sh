#1/bin/bash -e
# Used to make perturbed FWEA23 model for Simblast work
# 2/10/26
# 1) Update CONFIGS/* for the shallow (perturbed) model and rest of mesh
# 2) Make sure CONFIGS/* match each other in all parameters except tomo resi
# 3) Make sure you have nothing in EXPORT/ that you want to keep
# 4) Run this script to create mesh and tomo files and perturb shallow tomo
# 5) Make manual edits of the `Mesh_Par_file` to set PML layer thickness etc.
# 	 - Thickness is given in unit meters
# 	 - Note down element sizes in whole model config
# 6) Manually upload to SPECFEM DATA/ directory
# 7) Modify SPECFEM Par_file to match `NPROC`

# Set up file system
rm -rf EXPORT
mkdir -p EXPORT/tomo_files
mkdir -p EXPORT/tomo_files_unperturbed
mkdir -p EXPORT/meshfem3D_files
rm -rf TMP/*
cp -r CONFIGS EXPORT

# Make tomography and mesh files
tomo-helper -c CONFIGS/whole_model_config_fwea23.yaml 
tomo-helper -c CONFIGS/shallow_only_config_fwea23.yaml

# Set up unperturbed model for ref. seismograms
cp -r TMP/WHOLE/tomo_files/tomography_model.xyz EXPORT/tomo_files_unperturbed/tomography_model_2.xyz
cp -r TMP/SHALLOW/tomo_files/tomography_model.xyz EXPORT/tomo_files_unperturbed/tomography_model_1.xyz

# Apply stochastic perturbations
python /Users/chow/Repos/simutils/meshing/stochperts/apply_stochpert_xyz.py \
	TMP/SHALLOW/tomo_files/tomography_model.xyz \
	EXPORT/tomo_files/tomography_model_1.xyz

# Move rest of files
cp -r TMP/WHOLE/meshfem3D_files/* EXPORT/meshfem3D_files
cp -r TMP/WHOLE/tomo_files/tomography_model.xyz EXPORT/tomo_files/tomography_model_2.xyz

# Echo reminders
echo "TOMO FILE 1 == PERTURBED SHALLOW MODEL"
echo "> TO DO: ADJUST THICKNESS PML LAYER"
echo "> TO DO: CHANGE MATERIAL_ID"
echo "vi EXPORT/meshfem3D_files/Mesh_Par_file"
echo "---"
echo "ADD_PML_AS_EXTRA_MESH_LAYERS    = .true."
echo "NUMBER_OF_PML_LAYERS_TO_ADD     = 3"
echo "---"
echo "NMATERIALS                      = 2"
echo "-1 tomography elastic tomography_model_1.xyz 0 2"
echo "-2 tomography elastic tomography_model_2.xyz 0 2"
echo "---"


