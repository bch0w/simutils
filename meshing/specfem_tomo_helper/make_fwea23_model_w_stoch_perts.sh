#!/bin/bash -e
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

SIMUTILS="/home/bhchow/REPOS/simutils"

# Set up file system
rm -rf EXPORT
mkdir -p EXPORT/tomo_files
mkdir -p EXPORT/tomo_files_unperturbed
mkdir -p EXPORT/meshfem3D_files
rm -rf TMP/*
cp -r CONFIGS EXPORT

# Make tomography and mesh files
tomo-helper -c CONFIGS/whole_model_config_FWEA23.yaml 
tomo-helper -c CONFIGS/shallow_only_config_FWEA23.yaml

# Set up unperturbed model for ref. seismograms
cp -r TMP/WHOLE/tomo_files/tomography_model.xyz EXPORT/tomo_files_unperturbed/tomography_model_whole.xyz
cp -r TMP/SHALLOW/tomo_files/tomography_model.xyz EXPORT/tomo_files_unperturbed/tomography_model_shallow.xyz

# Add attenuation to the unperturbed models
python ${SIMUTILS}/meshing/model_utils/append_attenuation.py EXPORT/tomo_files_unperturbed/tomography_model_whole.xyz
python ${SIMUTILS}/meshing/model_utils/append_attenuation.py EXPORT/tomo_files_unperturbed/tomography_model_shallow.xyz
 
# Apply stochastic perturbations. No depth cutoff
python ${SIMUTILS}/meshing/stochperts/apply_stochpert_xyz.py \
	TMP/SHALLOW/tomo_files/tomography_model.xyz \
	EXPORT/tomo_files/tomography_model_shallow.xyz  

# Move rest of files
cp -r TMP/WHOLE/meshfem3D_files/* EXPORT/meshfem3D_files
cp -r TMP/WHOLE/tomo_files/tomography_model.xyz EXPORT/tomo_files/tomography_model_whole.xyz

# Last minute fixes of the meshfem3D file
python ${SIMUTILS}/meshing/specfem_tomo_helper/modify_mesh_par_file.py

