# adjTomo on Wisteria

The following steps show you how to activate the shared adjTomo Conda
environment on Wisteria-Aquarius, and run an example problem.

## 0. Initialize the Conda environment

If this is your first time using Conda on Wisteria, you will need to 
perform a one-time setup that adds a conda initialization statement
to your bashrc file.

```bash
module load aquarius
module load miniconda/py39_4.9.2
conda init 
source ~/.bashrc 
```

## 1. Log on to the `prepost` node

To avoid running computation on the login nodes, we take advantage of the 
`prepost` node on Wisteria, which acts like the login node but does not burden 
shared resources (see Wisteria manual 5.2.3).

We submit an interactive job to `rscpgrp`=prepost in order to run the master 
job on the prepost node.

```bash
pjsub -X --interact -g gr58 -L rscgrp=prepost
```

## 2. Activate the Conda Environment

Now we can activate the adjTomo Conda environment. You will need to 
do this each time you log on to Wisteria and/or the prepost node.

```bash
module purge
module load aquarius
module load impi
module load miniconda/py39_4.9.2
conda activate /work/01/gr58/share/adjtomo/conda/envs/adjtomo
```

OR 

```bash
# This script does the same as the above commands
source /work/01/gr58/share/adjtomo/activate_conda_aquarius.sh
```

Test if this works by running the seisflows help message

```bash
seisflows -h
```

## 3. Run SeisFlows Example

Setup one of the SeisFlows examples on Wisteria-Aquarius. I have already compiled
SPECFEM2D so we can point to that directory for binary executables (`-r` flag).

```bash
cd /work/01/gr58/share/adjtomo/testing
mkdir test_seisflows_example  # or choose your own directory name
cd test_seisflows_example
seisflows examples setup 2 -r /work/01/gr58/share/adjtomo/REPOSITORIES/specfem2d
```

The above commands have set up your 2D example problem and created pre-defined
parameter file. 

We need to make adjustments to the default parameter file to use the Wisteria 
system module.

```bash
seisflows swap system wisteria
seisflows par group gr58
seisflows par rscgrp debug-a  # running jobs on the Aquarius debug node
```

And now we can run the example using:

```bash
seisflows submit
```

The master job will run directly on the login node and log messages will be 
printed to stdout. Simulation jobs will be submitted to the Aquarius debug
node.

You will know when the example problem completes successfully when you
see the following log message:

```
CLEANING WORKDIR FOR NEXT ITERATION
--------------------------------------------------------------------------------
2023-03-23 09:46:50 (I) | optimization has been restarted, defaulting to standard inversion workflow
2023-03-23 09:46:55 (I) | finished all 7 tasks in task list
2023-03-23 09:46:55 (I) | 
////////////////////////////////////////////////////////////////////////////////
                             COMPLETE ITERATION 02                              
////////////////////////////////////////////////////////////////////////////////
2023-03-23 09:46:55 (I) | setting current iteration to: 3
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

5/23/24 Updating AdjTomo
Old commit hashes are
SeisFlows: 7cfea606 (HEAD -> feature-wisteria_gpu, origin/feature-wisteria_gpu) okay now it works, does not like backslash newlines in bash commands
Pyatoa: 1a6b86a (HEAD -> devel, origin/devel) Remove abstraction from windowing procedure (#41)
788186d (HEAD -> devel, origin/devel) Merge branch 'master' into devel

## Useful links
https://slurm.schedmd.com/rosetta.pdf
https://www.nrel.gov/hpc/assets/pdfs/pbs-to-slurm-translation-sheet.pdf
https://staff.cs.manchester.ac.uk/~fumie/internal/RikenDocuments/Job_Operation_Software_en.pdf
