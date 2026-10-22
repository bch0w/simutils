# SPECFEM on Chinook
Created 2025-10-21

## Important Links and Info
- UAF RCS Documentation for Chinook: https://uaf-rcs.gitbook.io/uaf-rcs-hpc-docs/hpc
- Bryant's Chinook bashrc with helpful alias's and ways to check queue and file: https://github.com/bch0w/simutils/blob/master/cluster/profiles/bashrc_chinook
- Bryant's Chinook SBATCH scripts for SPECFEM3D: https://github.com/bch0w/simutils/tree/master/cluster/runscripts/chinook/specfem3d
- Slurm Cheatsheet: https://slurm.schedmd.com/pdfs/summary.pdf
- `qmap`: A useful utility for checking cluster and queue usage


## To Do (on Chinook)
1. Load modules. The Intel compiler has been tested and works well but be aware there are other compiler suites we can use. Suggest putting the following in your `bashrc` so you don't have to remember to do it everytime.
```
module purge
module load slurm
module load intel/2022a
```
2. Clone SPECFEM somewhere in your `home` directory (either direct or from your own fork)
```
cd <path/to/your/repositories>
git clone https://github.com/SPECFEM/specfem2d.git --branch=devel
cd specfem2d
```
3. Configure and compile binaries
```
./configure F90=ifort FC=ifort MPIFC=mpiifort CC=icc MPICC=mpiicc
make all
```
4. Set up your SPECFEM working directory (`bin/`, `DATA/`, `OUTPUT_FILES`) somewhere on `center1` (`cd /center1/ERTHQUAK/$USER`).
   Here is an example of how I like to do it: `/import/c1/ERTHQUAK/bhchow/scratch/example_specfem3d_dir`
6. Set up your SBATCH scripts similar to: https://github.com/bch0w/simutils/tree/master/cluster/runscripts/chinook/specfem3d
7. Set up SPECFEM for your run and run your simulations with `sbatch <SBATCH_SCRIPT>.sh`, you can check the status of your jobs with `squeue -A $USER`
8. Rsync results back to your local machine for analysis. If you put the following command on your laptop's `bashrc`, it makes it easier:
```
alias rsh='rsync --ignore-existing -avzhLe ssh'
```
You can then run this from your laptop to quickly retrieve files from Chinook. Make sure you have `chinook` defined in `~/.ssh/config` and have SSH keys setup:
```
rsh chinook:/path/to/data/*
```
or, a more specific example for my case:
```
rsh chinook:/import/c1/ERTHQUAK/bhchow/scratch/example_specfem3d_dir/OUTPUT_FILES/*.ascii
```
