# adjTomo on Wisteria

## 1. Activate the Conda environment and load the correct environment

```bash
module purge
module load aquarius
module load impi
module load miniconda/py39_4.9.2
conda activate /work/01/gr58/share/adjtomo/conda/envs/adjtomo
```
OR 
```bash
source /work/01/gr58/share/adjtomo/activate_conda_aquarius.sh
```

> __NOTE__: If this is your first time running Conda, you will need to run the
  following commands after loading miniconda/before activating your environment

```bash
conda init
source ~/.bashrc
```

Test if this works by running the seisflows help message

```bash
seisflows -h
```

## 2. Run SeisFlows Example

Setup one of the SeisFlows examples on Wisteria-Aquarius 

```bash
cd /work/01/gr58/share/adjtomo/testing
mkdir test_seisflows_example  # or choose your own directory name
cd test_seisflows_example
seisflows examples setup 2 -r /work/01/gr58/share/adjtomo/REPOSITORIES/specfem2d
```

This will have set up your 2D example problem and created a parameter file. 
We need to make adjustments to the default parameter file to use the Wisteria
system:

```bash
seisflows swap system wisteria
seisflows par group gr58
seisflows par rscgrp debug-a  # running jobs on the Aquarius debug node
```

And now we can run the example

```bash
seisflows submit
```

## BRYANT'S INTERNAL NOTES

### 1. Creating Conda environment 

```bash
git clone --branch devel https://github.com/adjtomo/seisflows.git
cd seisflows
conda env create --prefix /work/01/gr58/share/adjtomo/conda/envs/adjtomo -f environment.yml
```

> __NOTE__: Use the --prefix command to install in a nonstandard location 
  SeisFlows will install all the other adjTomo packages

### 2. Setting up Wisteria

Q: Can't submit bulk jobs?

```bash
pjsub --bulk --sparam '1-3' test_submit.sh   
[ERR.] PJM 0070 pjsub No execute permission: pjsub-bulk.
```

A: no access for bulk jobs on Wisteria

Q: Can't access Conda/adjTomo on compute node
```bash
pjsub --interact -g gr58 -L rscgrp=interactive-a,node=1
module load aquarius
module load ompi/4.1.1
source /work/opt/local/x86_64/cores/miniconda/py39_4.9.2/bin/activate
conda activate /work/01/gr58/share/adjtomo/conda/envs/adjtomo
```

Difficulties:
- Can't submit batch job from a compute node
- Can't SSH to login node from compute node
- Can't provide command line arguments to run script when running pjsub


## Useful links
https://slurm.schedmd.com/rosetta.pdf
https://www.nrel.gov/hpc/assets/pdfs/pbs-to-slurm-translation-sheet.pdf
https://staff.cs.manchester.ac.uk/~fumie/internal/RikenDocuments/Job_Operation_Software_en.pdf
