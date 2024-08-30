import os
import sys
from paraview.simple import *


# Input parameters
state = sys.argv[1] 
path_out = "/Users/chow/Work/inversions/nakversion/mkvtk/auto_slice_figures"
name = os.path.basename(state).split(".")[0]  # e.g., MODEL_11_reg_1_vsh
variable = name.split("_")[-1]  # e.g., vsh
depths = [10, 30, 50, 75, 100, 125, 150]# , 200, 250]
resolution = [1724, 1162]

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
            "100km": [4.35 4.55],
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



# Set up output directory
full_path_out = os.path.join(path_out, name)
if not os.path.exists(full_path_out):
    os.makedirs(full_path_out)

# load state and set view
LoadState(state)
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

# find source
for depth in depths:
    depth_str = f"{depth}km"
    vmin, vmax = depth_clims[depth_str]

    src = FindSource(depth_str)
    SetActiveSource(src)

    # toggle interactive widget visibility (only when running from the GUI)
    HideInteractiveWidgets(proxy=src.SliceType)

    # toggle interactive widget visibility (only when running from the GUI)
    HideInteractiveWidgets(proxy=src.HyperTreeGridSlicer)

    # get color transfer function/color map for 'vsh'
    vshLUT = GetColorTransferFunction('vsh')

    # get opacity transfer function/opacity map for 'vsh'
    vshPWF = GetOpacityTransferFunction('vsh')

    # get 2D transfer function for 'vsh'
    vshTF2D = GetTransferFunction2D('vsh')

    # Rescale transfer function
    vshLUT.RescaleTransferFunction(vmin, vmax)

    # Rescale transfer function
    vshPWF.RescaleTransferFunction(vmin, vmax)

    # Rescale 2D transfer function
    vshTF2D.RescaleTransferFunction(vmin, vmax, 0.0, 1.0)

    # get layout
    layout1 = GetLayout()

    # layout/tab size in pixels
    layout1.SetSize(resolution[0], resolution[1])

    # save screenshot
    fid_out = os.path.join(full_path_out, f"{variable}_{depth:0>3}km.png")
    SaveScreenshot(fid_out, renderView1, ImageResolution=resolution,
                   TransparentBackground=1)
    
    # Hide this 
    Hide(src, renderView1)
