# Setting up a SPECFEM3D-SeisFlows example on Wisteria GPUs
Feb. 19, 2024

Notepad related to setting up an example problem on Wisteria GPUs using 
SeisFlows (v3.0.0-beta) and SPECFEM3D Cartesian (v4.1.0). 

## Setup

- Copied Samriddhi's DATA/ directory from `/work/gc62/share/Samriddhi/Koketsu_ALL/`
- Renamed CMTSOLUTION files to match SeisFlows format (i.e., CMTSOLUTION\_<ID>)
- Ran SPECFEM binaries to generate DATABASE and waveform files
- Setup parameter file following previous gpu example: `/home/r58003/work/adjtomo/testing/gpu_test`


## Make SPECFEM3D for GPU

```bash
module load cuda/12.2
module load gcc
module load ompi-cuda
./configure --with-cuda=cuda11 MPIFC=mpifort 'FCFLAGS=-O3 -march=native' 'CFLAGS=-O3 -march=native' 'CUDA_FLAGS=-O3 -Impi'
make all
```

## Run SeisFlows on Wisteria (GPU)

Run SeisFlows from an interactive environment to ensure no computation done on login nodes
Must load required modules and the correct Conda environment


```bash
pjsub -X --interact -g gr58 -L rscgrp=prepost
module load cuda/12.2
module load gcc
module load ompi-cuda
conda activate /work/01/gr58/share/adjtomo/conda/envs/adjtomo
seisflows submit
```

## Misc. Bugs
- Runscripts failed to locate relative paths for executables (./bin/xgenerate_databases)
  so I had to provide the full path. This was due to symlinks not evaluating
  the same for compute nodes, fix is to copy executables not symlink


### Problem with Cartopy unable to plot maps
- Due to different /home/ directories adopted on the compute nodes, some default
  Python pathing related to Matplotlib and Cartopy were not working
- Using the `environs` parameter in SeisFlows, I was able to overwrite the default
  paths in order to allow the compute nodes the required read/write access
- MPLCONFIGDIR needs to be a writeable directory but was pointing at some cache
  which was not writeable, we set this to /tmp/ to make it work
- Cartopy was trying to save shapefile data to /pjhome/.local on the compute node
  which was also not writeable. I needed to set the default path in the actual
  source code to a location that was writeable. This is pretty fragile so any
  updates to Cartopy may break this implementation.
- You will need to change path the _data_dir in `/home/r58003/work/adjtomo/conda/envs/adjtomo/lib/python3.11/site-packages/cartopy/__init__.py` to '/work/01/gr58/share/adjtomo/cartopy'


To deal with this I looked at the Cartopy source
  code to determine what was necessary to redirect the path. Using the env
  variable XDG_DATA_HOME I was able to redirect this to a local path. 
  Hopefully the sys admins don't complain that all the concurrent processes
  are trying to read from the filesystem at the same time
- !!! Nevermind XDG_DATA_HOME is required for other processes and leads to a 
  seg fault on SPECFEM


```python
# /home/r58003/work/adjtomo/conda/envs/adjtomo/lib/python3.11/site-packages/cartopy/__init__.py
#
# for the writable data directory (i.e. the one where new data goes), follow
# the XDG guidelines found at
# https://standards.freedesktop.org/basedir-spec/basedir-spec-latest.html
_writable_dir = os.path.join(os.path.expanduser('~'), '.local', 'share')
_data_dir = os.path.join(os.environ.get("XDG_DATA_HOME", _writable_dir),
                         'cartopy')
_cache_dir = os.path.join(tempfile.gettempdir(), 'cartopy_cache_dir')

config = {'pre_existing_data_dir': os.environ.get('CARTOPY_DATA_DIR', ''),
          'data_dir': _data_dir,
          'cache_dir': _cache_dir,
          'repo_data_dir': os.path.join(os.path.dirname(__file__), 'data'),
          'downloaders': {},
          }
```	


