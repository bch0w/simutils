#!/bin/bash -e                                                                   
                                                                                 
#PJM -N test_submit                                                           
#PJM -L rscgrp=regular-a                                                         
#PJM -L node=1                                                                   
#PJM -L elapse=00:00:05
#PJM -o output.log
#PJM -e output.err
#PJM -g gr58   

echo 'hello world ${PJM_BULKNUM}'
