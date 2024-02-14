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
# Sum gradient or model values
import time

# input_path = os.path.join(workflow.path.eval_grad, "gradient"); tag="kernel"
input_path = os.path.join(workflow.path.output, "GRADIENT_02"); tag="kernel"
# input_path = os.path.join(workflow.path.eval_grad, "kernels", "ZZ"); tag="kernel"

if tag in ["kernel", "gradient"]:
    parameters = [f"{par}_kernel" for par in solver._parameters]
elif tag in ["model"]:
    parameters = [par for par in solver._parameters]

parameters = [_[5:] for _ in parameters]  # strip reg_1 from parameters
output_path = "/import/c1/ERTHQUAK/bhchow/work/seisflows/SYNVERSION/output/VTK"
solver.combine_vol_data_vtk(input_path=input_path, 
                            output_path=output_path, 
                            parameters=parameters)
time.sleep(5)
for par in parameters:
    src = os.path.join(output_path, f"reg_1_{par}.vtk")
    dst = os.path.join(output_path, f"{tag}_{par}.vtk")
    os.rename(src, dst)
