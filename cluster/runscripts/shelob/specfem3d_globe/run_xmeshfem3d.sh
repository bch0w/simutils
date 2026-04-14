#!/bin/bash

#SBATCH --job-name=xmeshfem3D
#SBATCH --ntasks=48
#SBATCH --partition=defq
#SBATCH --time=00:10:00
#SBATCH --output=meshfem3D_%j.out


# Get the number of processors from the Par_file, ignore comments
BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `
NPROC_XI=`grep ^NPROC_XI DATA/Par_file | cut -d = -f 2 `
NPROC_ETA=`grep ^NPROC_ETA DATA/Par_file | cut -d = -f 2` 
NCHUNKS=`grep ^NCHUNKS DATA/Par_file | cut -d = -f 2 ` 
NPROC=$(( $NCHUNKS * $NPROC_XI * $NPROC_ETA )) 

# remake the output directories assuming we are starting fresh
rm -rf OUTPUT_FILES
mkdir -p OUTPUT_FILES

rm -rf ${BASEMPIDIR}
mkdir -p ${BASEMPIDIR}

echo "xmeshfem3D on ${NPROC} processors"
echo
echo "`date`"
time mpiexec -n ${NPROC} ./bin/xmeshfem3D
echo
echo "finished at: `date`"

