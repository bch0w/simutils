"""
A script to control Paraview using Python, setting up a figure and taking 
screenshots at various depths or cross-sections through the model.
Can be run in the Paraview Python Shell, or using PvPython from the CLI
"""
import os
import sys
from glob import glob
from paraview.simple import *


# 'Global' constants set here
COLOR_TABLE = {"k": [0., 0., 0.], "w": [1., 1., 1.], "r": [1., 0., 0.],
               "b": [0., 0., 1.], "g": [0., 1., 0.], "y": [1., 1., 0.],
               "o": [1., .5, 0.], "c": [0., 1., 1.], "gray": [.5, .5, .5]}

FONTSIZE = 15
COLOR = COLOR_TABLE["k"]
"""
Preset parameters that should stay constant among certain file types

:cbar_title (str): Label for the colorbar
:colormap (str): Colormap to define LUT
:invert_cmap (bool): Invert the original bounds of the colormap
:center_cmap (bool): Ensure that the middle value of the colormap is 0
:range_label_format (str): String formatter for the range bounds on colorbar
:round_base (int): Round range labels to a base value, if None, no rounding
:constant_bounds (bool): Keep the colorbar bounds constant for each screenshot
"""
PRESETS = {
        "model_vp": 
            {"cbar_title": "Vp [m/s]",
             "colormap": "Rainbow Desaturated",
             "invert_cmap": True,
             "center_cmap": False,
             "range_label_format": "%.1f",
             "round_base": 10,
             "constant_bounds": False,
             },
        "model_vs": 
            {"cbar_title": "Vs [m/s]",
             "colormap": "Rainbow Desaturated",
             "invert_cmap": True,
             "center_cmap": False,
             "range_label_format": "%.1f",
             "round_base": 10,
             "constant_bounds": False,
             },
        "gradient_vp_kernel": 
            {"cbar_title": "Vp Gradient [m^-2 s^2]",
             "colormap": "Cool to Warm (Extended)",
             "invert_cmap": True,
             "center_cmap": True,
             "range_label_format": "%.1E",
             "round_base": None,
             "constant_bounds": True,
             },
        "gradient_vs_kernel": 
            {"cbar_title": "Vs Gradient [m^-2 s^2]",
             "colormap": "Cool to Warm (Extended)",
             "invert_cmap": True,
             "center_cmap": True,
             "range_label_format": "%.1E",
             "round_base": None,
             "constant_bounds": True,
             },
        "update_vp": 
            {"cbar_title": "Vp Net Update [ln(m/m00)]",
             "colormap": "Blue Orange (divergent)",
             "invert_cmap": True,
             "center_cmap": True,
             "range_label_format": "%.02f",
             "round_base": None,
             "constant_bounds": True,
             },
        "update_vs": 
            {"cbar_title": "Vs Net Update [ln(m/m00)]",
             "colormap": "Blue Orange (divergent)",
             "invert_cmap": True,
             "center_cmap": True,
             "range_label_format": "%.02f",
             "round_base": None,
             "constant_bounds": True,
             },
        "ratio_poissons": 
            {"cbar_title": "Poisson's Ratio",
             "colormap": "Blue - Green - Orange",
             "invert_cmap": False,
             "center_cmap": False,
             "range_label_format": "%.1f",
             "round_base": None,
             "constant_bounds": True,
             },
        "ratio_vpvs": 
            {"cbar_title": "Vp/Vs Ratio",
             "colormap": "Cool to Warm (Extended)",
             "invert_cmap": False,
             "center_cmap": False,
             "range_label_format": "%.1f",
             "round_base": None,
             "constant_bounds": True,
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
    RenameSource(f"z_{depth}km", slice_vtk)

    # Hide the plane widget
    Hide3DWidgets(proxy=slice_vtk)

    return slice_vtk

def make_ruler(point1, point2, label=""):
    """
    Generate a ruler to be used as a scalebar
    """
    renderView = GetActiveViewOrCreate("RenderView")

    ruler = Ruler()
    ruler.Point1 = point1
    ruler.Point2 = point2
    rulerDisplay = Show(ruler, renderView, "RulerSourceRepresentation")
    rulerDisplay.LabelFormat = label
    rulerDisplay.AxisColor = COLOR
    rulerDisplay.Color = COLOR
    rulerDisplay.FontSize = FONTSIZE

    Hide3DWidgets(proxy=ruler)

    return ruler


# ============================== USER PARAMETERS ===============================

# File ID's specifying the chosen files
state_fid = "/Users/Chow/Documents/academic/vuw/forest/forest.pvsm"

# Common path for finding and saving files
# base_path = "/Users/Chow/Documents/academic/vuw/forest/birch/"
base_path = os.path.abspath("../..")  # assuming directory structure here

# Allow for user-input files on the command line
if len(sys.argv) > 1:
    data_fids = sys.argv[1:]
else:
    data_fids = glob(os.path.join(base_path, "vtk", "model", "ratio*0007*vs*"))

# Define where the slice the model in the XZ plane, surface implicit
slices = ["surface", 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 25, 30, 50]
slices = ["surface"]

# Save path for screenshots
path_out = os.path.join(base_path, "figures")
view_size = [1037, 813]
# ==============================================================================

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

    # Constantly used to render views
    renderView = GetActiveViewOrCreate("RenderView")
    renderView.ViewSize = view_size

    # ACTIVE STATE: Load an active state with glyphs, coastline etc.
    servermanager.LoadState(state_fid)
    for source in ["receivers", "src_epicenter", "coast"]:
        vtk_ = FindSource(source)
        vtkDisplay_ = Show(vtk_, renderView, "GeometryRepresentation")
        vtkDisplay_.SetRepresentationType("Surface With Edges")
        vtkDisplay_.LineWidth = 2.0
        vtkDisplay_.Opacity = 0.4
        if source == "src_epicenter":
            vtkDisplay_.AmbientColor = COLOR_TABLE["g"]
            vtkDisplay_.DiffuseColor = COLOR_TABLE["g"]
        elif source == "coast":
            vtkDisplay_.PointSize = 1.5

    # Plot the entire model volume
    vtk = OpenDataFile(data_fid)
    vol_vmin, vol_vmax = vtk.PointData.GetArray(0).GetRange(0)
    RenameSource("surface", vtk)
    Show(vtk)
    ResetCamera()

    # ANNOTATIONS: Depth slice annotation
    # Depth slice annotation
    text = Text()
    textDisplay = Show(text, renderView)
    textDisplay.Position = [0.625, 0.1]
    textDisplay.FontSize = FONTSIZE * 2
    textDisplay.Color = COLOR

    # File name annotation
    text_tag = Text()
    text_tag.Text = save_fid
    textTagDisplay = Show(text_tag, renderView)
    textTagDisplay.Position = [0.65, 0.95]
    textTagDisplay.FontSize = FONTSIZE
    textTagDisplay.Color = COLOR

    Render()

    # RULER: Create rulers as scalebar, points determined manually
    ruler_origin = [188068, 5722936, 0] 
    ruler_x = [288068, 5722936, 0]
    ruler_y = [188068, 5622936, 0]

    make_ruler(ruler_origin, ruler_x, "100km")
    make_ruler(ruler_origin, ruler_y)

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
    cbar.TitleFontSize = FONTSIZE
    cbar.LabelFontSize = FONTSIZE
    cbar.TitleColor = COLOR
    cbar.LabelColor = COLOR

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
            tag = f"z_{slice_}km"

            # Determine colorbar bounds of slice to rescale colorbar
            vmin, vmax = slice_vtk.PointData.GetArray(0).GetRange(0)

        # Make sure that global bounds are set if required
        if preset["constant_bounds"]:
            vmin, vmax = vol_vmin, vol_vmax
    
        # Some fields should be centered on 0 despite the actual data bounds
        if preset["center_cmap"]:
            vabsmax = max(abs(vmin), abs(vmax))
            vmin, vmax = -1 * vabsmax, vabsmax

        # If desired, round the colobar bounds to some base value
        if preset["round_base"]:
            vmin = myround(vmin, preset["round_base"])
            vmax = myround(vmax, preset["round_base"])
       
        # Apply depth specific values
        vsLUT.RescaleTransferFunction(vmin, vmax)
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

    # CLEAR: ResetSession() wasnt getting rid of the text so delete manually
    Delete(text)
    del text
    ResetSession()

