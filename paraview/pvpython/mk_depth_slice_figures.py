import os
import sys
import math
from copy import copy
from glob import glob
from paraview.simple import *


def myround(x, base, choice=None):
    """
    Round values to the nearest base
    """
    if choice == "up":
        fx = math.ceil
    elif choice == "down":
        fx = math.floor
    else:
        fx = round
    return int(base * fx(float(x) / base))


# vUSER PARAMETERS HEREv
try:
    vtk_files = sys.argv[1:]
except IndexError:
    vtk_files = glob("vtk_files/MODEL_??_reg_1_vs?.vtk")
path_out = "./slice_figures"
template_state_file = "template.pvsm"
_lines_template = open(template_state_file).read()
temp_state_file = "temp_state_file.pvsm"
depths = [10, 30 , 50, 75, 100, 125, 150, 200, 250]
resolution = [1724, 1162]
auto_scale = True
auto_scale_base = 0.01

# Colorbar limits for each depth slice, determine manually
depth_clims = {
    "vsh": {"10km": [2.75, 4.0],
            "30km": [4.20, 4.50],
            "50km": [4.35, 4.85],
            "75km": [4.25, 5.65],
            "100km": [4.17, 4.50],
            "125km": [4.20, 4.50],
            "150km": [4.35, 4.45],
            "200km": [],
            "250km": []
            },
    "vsv": {"10km": [3.0, 4.0],
            "30km": [4.2, 5.3],
            "50km": [4.2, 4.75],
            "75km": [4.3, 4.55],
            "100km": [4.35, 4.55],
            "125km": [4.4, 4.5],
            "150km": [4.35, 4.45],
            "200km": [],
            "250km": []
            },
    "vpv": {"10km": [],
            "30km": [],
            "50km": [],
            "75km": [],
            "100km": [],
            "125km": [],
            "150km": [],
            "200km": [],
            "250km": []
            },
    "vpv": {"10km": [],
            "30km": [],
            "50km": [],
            "75km": [],
            "100km": [],
            "125km": [],
            "150km": [],
            "200km": [],
            "250km": []
            }
        }
# ^USER PARAMETERS HERE^  


# Loop through all VTK files chosen
for vtk_file in vtk_files:
    # Get some pertinent information about the file
    filename = os.path.basename(vtk_file)  # e.g., MODEL_14_reg_1_vsh.vtk
    fullname = filename.split(".")[0]  # e.g., MODEL_14_reg_1_vsh
    _, N, _, _, parameter = fullname.split("_")
    modelname = f"M{N}_{parameter.upper()}"

    print(filename)
        
    # Set up output directory based on the current model
    full_path_out = os.path.join(path_out, modelname)
    if not os.path.exists(full_path_out):
        os.makedirs(full_path_out)

    # Create the temporary state file to be used that has the correct 
    # naming and text scheme
    lines = copy(_lines_template).format(PAR=parameter, 
                                         FID=fullname)
    with open(temp_state_file, "w") as f:
        f.write(lines)

    # load state and set view
    LoadState(temp_state_file)
    renderView1 = FindViewOrCreate('RenderView1', viewtype='RenderView')
    SetActiveView(renderView1)

    # Set current camera placement for renderView1 and don't change
    renderView1.CameraPosition = [
        -0.45875962494824524, -0.23577305053017142, 1.284718733525635
        ]
    renderView1.CameraFocalPoint = [
        0.4576892617589983, 0.26903787019648195, -1.3025221649023777
        ]
    renderView1.CameraViewUp = [
        0.8329424388283602, 0.40735480689687115, 0.3745249723271819
        ]
    renderView1.CameraParallelScale = 0.4933477880175255

    # Find each depth slice
    for depth in depths:
        depth_str = f"{depth}km"
        print(f"\t{depth}")

        # Change the annotation text
        depth_text = FindSource("depth_text")
        depth_text.Text = f"{modelname} @ {depth_str.upper()}"

        # Set the depth slice in question
        src = FindSource(depth_str)
        SetActiveSource(src)

        # Variable is used to determine the name of the color object
        variable = src.PointData.GetArray(0).Name

        # Turn on view and colorbar
        Display = Show(src, renderView1, 'GeometryRepresentation')
        Display.SetScalarBarVisibility(renderView1, True)

        # toggle interactive widget visibility (only when running from the GUI)
        HideInteractiveWidgets(proxy=src.SliceType)

        # toggle interactive widget visibility (only when running from the GUI)
        HideInteractiveWidgets(proxy=src.HyperTreeGridSlicer)

        # Either manually set vmin vmax or round off the Data Range values
        if auto_scale:
            active = GetActiveSource()
            Display = GetDisplayProperties(src, view=renderView1)
            Display.RescaleTransferFunctionToDataRange(False, True)
            vmin, vmax = active.PointData.GetArray(0).GetRange(0)
            
            cvt = 1E2  # convert km -> m and back
            vmin = myround(vmin * cvt, auto_scale_base * cvt, "down") * 1 / cvt
            vmax = myround(vmax * cvt, auto_scale_base * cvt, "up") * 1 / cvt
        else:
            vmin, vmax = depth_clims[variable][depth_str]

        # Rescale the colorscale to a new vmin and vmax
        vshPWF = GetOpacityTransferFunction(variable)
        vshLUT = GetColorTransferFunction(variable)
        vshTF2D = GetTransferFunction2D(variable)

        vshPWF.RescaleTransferFunction(vmin, vmax)
        vshLUT.RescaleTransferFunction(vmin, vmax)
        vshTF2D.RescaleTransferFunction(vmin, vmax, 0.0, 1.0)

        # get layout
        layout1 = GetLayout()
        layout1.SetSize(resolution[0], resolution[1])

        # save screenshot
        fid_out = os.path.join(full_path_out, f"{variable}_{depth:0>3}km.png")
        SaveScreenshot(fid_out, renderView1, ImageResolution=resolution,
                       TransparentBackground=1)
        
        # Hide this 
        Hide(src, renderView1)

    if os.path.exists(temp_state_file):
        os.remove(temp_state_file)
