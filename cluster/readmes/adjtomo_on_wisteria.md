# Installing adjTomo on Wisteria

## To activate the Conda environment

```bash
module load aquarius
module load miniconda/py38_4.9.2
conda activate /work/01/gr58/share/adjtomo/conda/envs/adjtomo
```
OR 
```bash
source /home/r58003/work/adjtomo/activate_conda.sh
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

## Useful links
https://slurm.schedmd.com/rosetta.pdf
https://www.nrel.gov/hpc/assets/pdfs/pbs-to-slurm-translation-sheet.pdf
https://staff.cs.manchester.ac.uk/~fumie/internal/RikenDocuments/Job_Operation_Software_en.pdf
