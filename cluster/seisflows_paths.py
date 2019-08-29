export PYTHONUNBUFFERED=true  # quickly get Python output logs
export HDF5_USE_FILE_LOCKING='FALSE'  # trying to avoid (unable to lock file, h5py error)
export PATH=$PATH:/scale_wlg_persistent/filesets/project/nesi00263/PyPackages/seisflows/scripts
export PYTHONPATH=$PYTHONPATH:/scale_wlg_persistent/filesets/project/nesi00263/PyPackages/seisflows
export PYTHONPATH=$PYTHONPATH:/scale_wlg_persistent/filesets/project/nesi00263/PyPackages/seisflows-hpc
echo "seisfows paths exported"
