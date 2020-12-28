"""
A script to control Paraview using Python, setting up a figure and taking 
screenshots at various depths or cross-sections through the model.
Can be run in the Paraview Python Shell, or using PvPython from the CLI
"""
import os
import sys
from glob import glob
from paraview.simple import *

"""
Preset parameters that should stay constant among certain file types
"""
PRESETS = {
        "model_vp": 
            {"cbar_title": "Vp [m/s]",
             "colormap": "Rainbow Desaturated",
             "invert_cmap": True,
             "center_cmap": False,
             "range_label_format": "%.1f",
             "round_base": 10,
             },
        "model_vs": 
            {"cbar_title": "Vs [m/s]",
             "colormap": "Rainbow Desaturated",
             "invert_cmap": True,
             "center_cmap": False,
             "range_label_format": "%.1f",
             "round_base": 10,
             },
        "gradient_vp_kernel": 
            {"cbar_title": "Vp Gradient [m^-2 s^2]",
             "colormap": "Cool to Warm (Extended)",
             "invert_cmap": False,
             "center_cmap": True,
             "range_label_format": "%.1E",
             "round_base": None,
             },
        "gradient_vs_kernel": 
            {"cbar_title": "Vs Gradient [m^-2 s^2]",
             "colormap": "Cool to Warm (Extended)",
             "invert_cmap": False,
             "center_cmap": True,
             "range_label_format": "%.1E",
             "round_base": None,
             },
        "update_vp": 
            {"cbar_title": "Vp Net Update [%]",
             "colormap": "Blue Orange (divergent)",
             "invert_cmap": False,
             "center_cmap": True,
             "range_label_format": "%.03f",
             "round_base": None,
             },
        "update_vs": 
            {"cbar_title": "Vs Net Update [%]",
             "colormap": "Blue Orange (divergent)",
             "invert_cmap": False,
             "center_cmap": True,
             "range_label_format": "%.03f",
             "round_base": None,
             },
        "ratio_poissons": 
            {"cbar_title": "Poisson's Ratio",
             "colormap": "Blue - Green - Orange",
             "invert_cmap": False,
             "center_cmap": False,
             "range_label_format": "%.1f",
             "round_base": None,
             },
        "ratio_vpvs": 
            {"cbar_title": "Vp/Vs Ratio",
             "colormap": "",
             "invert_cmap": False,
             "center_cmap": False,
             "range_label_format": "%.01f",
             "round_base": None,
             },
        }


def myround(x, base):
    """
    Round values to the nearest base
    """
    return int(base * round(float(x) / base))


def slice_plane(vtk, depth):
    """
    Slice a plane through a volume
    """
    # Generate the Z-normal slice and set at a specific depth
    slice_vtk = Slice(vtk)
    slice_vtk.SliceType = "Plane"
    slice_vtk.SliceType.Normal = [0., 0., 1.]
    origin = slice_vtk.SliceType.Origin
    origin = [origin[0], origin[1], abs(depth) * -1E3]
    slice_vtk.SliceType.Origin = origin

    # Apply the slice and render the new view
    renderView = GetActiveViewOrCreate("RenderView")
    sliceDisplay = Show(slice_vtk, renderView)
    RenameSource(f"{depth}km", slice_vtk)

    # Hide the plane widget
    Hide3DWidgets(proxy=slice_vtk)

    return slice_vtk

# ============================== USER PARAMETERS ===============================
# File ID's specifying the chosen files
state_fid = "/Users/Chow/Documents/academic/vuw/forest/forest.pvsm"

# Allow for user-input files on the command line
if len(sys.argv) > 1:
    data_fids = sys.argv[1:]
else:
    data_fids = glob("/Users/Chow/Documents/academic/vuw/forest/birch/vtk/"
                     "model/model*0007*vs*")

# Define where the slice the model in the XZ plane, surface implicit
slices = ["surface", 2, 4, 6, 8, 10, 15, 20, 25, 30, 50]

# Save path for screenshots
path_out = "/Users/Chow/Documents/academic/vuw/forest/birch/figures"
view_size = [1037, 813]
# ==============================================================================

# Constantly used to render views
renderView = GetActiveViewOrCreate("RenderView")
renderView.ViewSize = view_size

# ACTIVE STATE: Load an active state that has common traits, e.g. src-rcv glyphs
servermanager.LoadState(state_fid)
for source in ["receivers", "src_epicenter", "coast"]:
    vtk_ = FindSource(source)
    vtkDisplay_ = Show(vtk_, renderView, "GeometryRepresentation")
    vtkDisplay_.SetRepresentationType("Surface With Edges")
    vtkDisplay_.LineWidth = 2.0
    vtkDisplay_.Opacity = 0.4
    if source == "src_epicenter":
        vtkDisplay_.AmbientColor = [0., 1., 0.]  # Green
        vtkDisplay_.DiffuseColor = [0., 1., 0.]
    elif source == "coast":
        vtkDisplay_.PointSize = 1.5

for data_fid in data_fids:
    # Define the output directory to save figures
    save_fid = os.path.splitext(os.path.basename(data_fid))[0]
    full_path_out = os.path.join(path_out, save_fid)
    if not os.path.exists(full_path_out):
        os.makedirs(full_path_out)

    # Presets defined by labels assigned to the file names, ignoring model num
    preset = save_fid.split("_")
    preset = "_".join(preset[:1] + preset[2:])
    assert(preset in PRESETS), f"{preset} does not match preset keys"
    preset = PRESETS[preset]

    # Plot the entire model volume
    vtk = OpenDataFile(data_fid)
    RenameSource("surface", vtk)
    Show(vtk)
    ResetCamera()

    # ANNOTATION: Create annotation to be filled in at each depth slice
    text = Text()
    textDisplay = Show(text, renderView)
    textDisplay.Position = [0.65, 0.1]
    textDisplay.FontSize = 40
    Render()

    # CAMEAR: Move camera a bit closer
    camera = GetActiveCamera()
    camera.Dolly(1.2)

    # COLORMAP: Change the colormap to the desired preset value
    quantity = vtk.PointData.GetArray(0).Name  # e.g. model_init_vp
    vsLUT = GetColorTransferFunction(quantity)
    vsLUT.ApplyPreset(preset["colormap"], True)
    vsLUT.UseAboveRangeColor = 0
    vsLUT.UseBelowRangeColor = 0
    if preset["invert_cmap"]:
        vsLUT.InvertTransferFunction()

    # COLORBAR: Create the colorbar and set a common look
    cbar = GetScalarBar(vsLUT)
    cbar.Title = preset["cbar_title"]
    cbar.ComponentTitle = ""
    cbar.AutoOrient = 0
    cbar.Orientation = "Vertical"
    cbar.Position = [0.249, 0.78]
    cbar.RangeLabelFormat = preset["range_label_format"]
    cbar.AddRangeLabels = 1
    cbar.ScalarBarThickness = 35
    cbar.ScalarBarLength = 0.15
    cbar.TitleFontSize = 15
    cbar.LabelFontSize = 15

    # SLICE: Generate slices through the volume at desired depth levels
    slice_vtk = None
    for slice_ in slices:
        if slice_ == "surface":
            tag = "map"

            # Hacky method to get data range by selecting surface points
            v = GetActiveView()
            source = FindSource("surface")
            SelectSurfacePoints(Rectangle=[0, 0, v.ViewSize[0], v.ViewSize[1]])
            extract = ExtractSelection(Input=source)
            vmin, vmax = extract.PointData.GetArray(0).GetRange(0)

            # Get rid of the extraction surface
            Delete(extract)
            del extract
            SetActiveSource(source)
            ClearSelection()
        else:
            slice_vtk = slice_plane(vtk, slice_)
            tag = f"{slice_}km"

            # Determine colorbar bounds of slice to rescale colorbar
            vmin, vmax = slice_vtk.PointData.GetArray(0).GetRange(0)

        # Some fields should be centered on 0 despite the actual data bounds
        if preset["center_cmap"]:
            vabsmax = max(abs(vmin), abs(vmax))
            vmin, vmax = -1 * vabsmax, vabsmax
        # If desired, round the colobar bounds to some base value
        if preset["round_base"]:
            vmin = myround(vmin, preset["round_base"])
            vmax = myround(vmax, preset["round_base"])

        vsLUT.RescaleTransferFunction(vmin, vmax)

        # Annotate the depth level
        text.Text = tag

        # Ensure that the colorbar shows for each new slice
        active = GetActiveSource()
        display = GetDisplayProperties(active, view=renderView)
        display.SetScalarBarVisibility(renderView, True)

        # SCREENSHOT: Save screenshot, hide slice and move on
        renderView = GetActiveViewOrCreate("RenderView")
        SaveScreenshot(os.path.join(full_path_out, f"{tag}.png"), 
                       renderView, ImageResolution=view_size, 
                       TransparentBackground=1)
        if slice_vtk:
            Hide(slice_vtk, renderView)
        else:
            Hide(vtk, renderView)

