# Copy-paste into SeisFlows debug to make VTK files for visualization

input_path = ""

# Sum individual event kernels
import time
kernels = ["ET", "NT", "ER", "NR"]
sources = ["AK_A21K"]
parameters = ["vsh_kernel", "vpv_kernel"]
output_path = "/import/c1/ERTHQUAK/bhchow/work/seisflows/SYNVERSION/output/VTK"

for kernel in kernels:
    #for srcname in solver.source_names:
    for srcname in sources:
        input_path = os.path.join(workflow.path.eval_grad, "kernels",
                                  srcname, kernel)
        solver.combine_vol_data_vtk(input_path=input_path, 
                                    output_path=output_path, 
                                    parameters=parameters)
        time.sleep(5)
        for par in parameters:
            src = os.path.join(output_path, f"reg_1_{par}.vtk")
            dst = os.path.join(output_path, f"{srcname}_{par}_{kernel}.vtk")
            os.rename(src, dst)

# =============================================================================
# Sum individual misfit kernels ZZ, RR, TT
import time

input_path = ""

kernels = ["ZZ", "RR", "TT"]
parameters = ["vsh_kernel", "vpv_kernel"]
output_path = "/import/c1/ERTHQUAK/bhchow/work/seisflows/SYNVERSION/output/VTK"

for kernel in kernels:
    input_path = os.path.join(workflow.path.eval_grad, "kernels", kernel)
    solver.combine_vol_data_vtk(input_path=input_path, 
                                output_path=output_path, 
                                parameters=parameters)
    time.sleep(5)
    for par in parameters:
        src = os.path.join(output_path, f"reg_1_{par}.vtk")
        dst = os.path.join(output_path, f"{kernel}_{par}.vtk")
        os.rename(src, dst)


# =============================================================================
# Sum all gradient or model values automatically
import time
import os
from glob import glob


# User-set input parameters 
set_input_path = "/import/c1/ERTHQUAK/bhchow/work/seisflows/NAKVERSION/output"
output_path = os.path.join(set_input_path, "VTK")
set_parameters = ["reg1_vsh"] or solver._parameters

# Used to get back after relative chdir'ing
cwd = os.getcwd()

if not os.path.exists(output_path):
    os.makedirs(output_path)

for name in ["GRADIENT", "MODEL"]:
    # Add _kernel to gradient/kernel parameter names to match naming schema
    if name == "GRADIENT":
        parameters = [f"{par}_kernel" for par in set_parameters]
    else:
        parameters = set_parameters

    # Strip 'reg1_' from parameter names because SeisFlows doesn't want that
    for parameter in parameters:
        if "reg" in parameter:
            parameters = [_[5:] for _ in parameters]  # strip reg1_ from param

    # Loop through all available gradients/models and generate VTK files
    for input_path in glob(os.path.join(set_input_path, f"{name}_??")):
        tag = os.path.basename(input_path)
        for par in parameters:
            dst = os.path.join(output_path, f"{tag}_{par}.vtk")
            if os.path.exists(dst):
                print(f"{os.path.basename(dst)} already exists, skipping")
                continue

            solver.combine_vol_data_vtk(input_path=input_path,
                                        output_path=output_path,
                                        parameters=parameters)

            time.sleep(5)
            src = os.path.join(output_path, f"reg_1_{par}.vtk")
            os.rename(src, dst)

os.chdir(cwd)
  

# =============================================================================
# Sum gradient, model or kernel by target directory

input_paths = {"MODEL_01": "/import/c1/ERTHQUAK/bhchow/work/seisflows/SYNVERSION/output/MODEL_01"}
output_path = "/import/c1/ERTHQUAK/bhchow/work/seisflows/SYNVERSION/output/VTKS"
set_parameters = ["reg1_vsh"] or solver._parameters

for name, input_path in inputs.items():
    if "KERNEL" in name or "GRADIENT" in name:
        parameters = [f"{par}_kernel" for par in set_parameters]
    else:
        parameters = set_parameters

    parameters = [_[5:] for _ in parameters]  # strip reg_1 from parameters
    solver.combine_vol_data_vtk(input_path=input_path, 
                                output_path=output_path, 
                                parameters=parameters)
    time.sleep(5)
    for par in parameters:
        src = os.path.join(output_path, f"reg_1_{par}.vtk")
        dst = os.path.join(output_path, f"{name}_{par}.vtk")
        os.rename(src, dst)
