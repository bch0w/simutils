#!/bin/bash -e                                                                   
                                                                                 
#PJM -N test_submit                                                           
#PJM -L rscgrp=debug-a   
#PJM -L node=1                                                                   
#PJM -L elapse=00:01:00
#PJM -o output.log
#PJM -e output.err
#PJM -g gr58   

echo "1"
module purge
module load miniconda/py38_4.9.2
source /work/opt/local/x86_64/cores/miniconda/py38_4.9.2/bin/activate /work/01/gr58/share/adjtomo/conda/envs/adjtomo
python python_script.py
echo "2"
