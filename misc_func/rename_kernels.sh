#!/bin/bash
echo "	renaming smoothed kernels for use in summation"
if ! [ -d INPUT_KERNELS/ ]
then
	echo "INPUT_KERNELS/ does not exit"
	exit
fi
if ! [ -f kernels_list.txt ]
then 
	echo "kernels_list.txt does not exist"
	exit
fi
while read EVENT_DIR
do
	rename kernel.bin kernel_original.bin ${EVENT_DIR}*kernel.bin
	rename kernel_smooth.bin kernel.bin ${EVENT_DIR}*kernel_smooth.bin
done < kernels_list.txt

