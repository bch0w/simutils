# Setting up a SPECFEM3D-SeisFlows example on Wisteria GPUs
Feb. 19, 2024

Notepad related to setting up an example problem on Wisteria GPUs using 
SeisFlows (v2.3.0) and SPECFEM3D Cartesian (v4.1.0). 

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

## Run Example

Run SeisFlows from an interactive environment to ensure no computation done on login node

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




