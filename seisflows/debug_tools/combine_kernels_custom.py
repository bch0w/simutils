"""
3/21/24 For Ambient Noise Adjoint Tomography. Currently SeisFlows does not
sum the individual horizontal kernels (RR=ER+NR, TT=ET+NT) because it does 
this all in one step (G=ER+NR+ET+NT+ZZ) to save computational time. However
it is often nice to look at the individual RR, TT and ZZ kernels to understand
what is going on. This debug function sums all these individual kernels 
together
"""
# >>> BEGIN COPY-PASTE DIRECTLY INTO SEISFLOWS DEBUG MODE
import os
from glob import glob

system.partition = "debug"
system.tasktime = 10
workflow.setup()

def combine_all_horizontal_kernels_to_misfit_kernel():
    """
    For iteration i, create RR and TT misfit kernels by summing all sources
    """
    parameters = ["reg1_vsh_kernel"] or \
            [f"{par}_kernel" for par in solver._parameters]

    set_output_path = os.path.join(workflow.path.output, "TEST")
    kernels_path = glob(os.path.join(workflow.path.output, "KERNELS_08"))

    # e.g., output/KERNEL_01 ...
    for kernel_dir in kernels_path:
        for kernel in ["R", "T"]:
            kernel_number = os.path.basename(kernel_dir)  # e.g., KERNEL_01

            input_paths = []
            output_path = os.path.join(set_output_path, kernel_number)

            # Collect all relevant event kernels for the given comp R or T
            for src in glob(os.path.join(kernel_dir, "*")):  # e.g., AK_A21K
                input_paths.append(os.path.join(src, f"E{kernel}"))
                input_paths.append(os.path.join(src, f"N{kernel}"))

            # Combines ER+NR (or T) kernels for all sources in a given iteration
            solver.combine(input_paths, output_path, parameters)

system.run([combine_all_horizontal_kernels_to_misfit_kernel], single=True)
# >>> END COPY-PASTE DIRECTLY INTO SEISFLOWS DEBUG MODE


