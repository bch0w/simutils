"""
A script to control Paraview using Python, to be run in Paraview Python Shell
"""
import os
import sys
import math
import string
from glob import glob
from paraview.simple import *



# 'Global' constants set here
COLOR_TABLE = {"k": [0., 0., 0.], "w": [1., 1., 1.], "r": [1., 0., 0.],
               "b": [0., 0., 1.], "g": [0., 1., 0.], "y": [1., 1., 0.],
               "o": [1., .5, 0.], "c": [0., 1., 1.], "gray": [.5, .5, .5]}

FONTSIZE = 15
COLOR = COLOR_TABLE["k"]
VIEW_SIZE = [1037, 813]

"""
Preset parameters that should stay constant among certain file types

:cbar_title (str): Label for the colorbar
:colormap (str): Colormap to define LUT
:invert_cmap (bool): Invert the original bounds of the colormap
:center_cmap (bool): Ensure that the middle value of the colormap is 0
:range_label_format (str): String formatter for the range bounds on colorbar
:round_base (int): Round range labels to a base value, if None, no rounding
:constant_bounds (bool): Keep the colorbar bounds constant for each screenshot
:num_table_values (int): Number of segmentations in the colorbar
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
             "num_table_values": 64,
             },
        "model_vs": 
            {"cbar_title": "Vs [m/s]",
             "colormap": "Rainbow Desaturated",
             "invert_cmap": True,
             "center_cmap": False,
             "range_label_format": "%.1f",
             "round_base": 10,
             "constant_bounds": False,
             "num_table_values": 64,
             },
        "gradient_vp_kernel": 
            {"cbar_title": "Vp Gradient [m^-2 s^2]",
             "colormap": "Cool to Warm (Extended)",
             "invert_cmap": True,
             "center_cmap": True,
             "range_label_format": "%.1E",
             "round_base": None,
             "constant_bounds": True,
             "num_table_values": 64,
             },
        "gradient_vs_kernel": 
            {"cbar_title": "Vs Gradient [m^-2 s^2]",
             "colormap": "Cool to Warm (Extended)",
             "invert_cmap": True,
             "center_cmap": True,
             "range_label_format": "%.1E",
             "round_base": None,
             "constant_bounds": True,
             "num_table_values": 64,
             },
        "update_vp": 
            {"cbar_title": "Vp Net Update [ln(m/m00)]",
             "colormap": "Blue Orange (divergent)",
             "invert_cmap": True,
             "center_cmap": True,
             "range_label_format": "%.02f",
             "round_base": None,
             "constant_bounds": True,
             "num_table_values": 64,
             },
        "update_vs": 
            {"cbar_title": "Vs Net Update [ln(m/m00)]",
             "colormap": "Blue Orange (divergent)",
             "invert_cmap": True,
             "center_cmap": True,
             "range_label_format": "%.02f",
             "round_base": None,
             "constant_bounds": True,
             "num_table_values": 64,
             },
        "ratio_poissons": 
            {"cbar_title": "Poisson's Ratio",
             "colormap": "Blue - Green - Orange",
             "invert_cmap": False,
             "center_cmap": False,
             "range_label_format": "%.1f",
             "round_base": None,
             "constant_bounds": True,
             "num_table_values": 64,
             },
        "ratio_vpvs": 
            {"cbar_title": "Vp/Vs Ratio",
             "colormap": "Cool to Warm (Extended)",
             # "colormap": "Rainbow Desaturated",
             "invert_cmap": False,
             "center_cmap": False,
             "range_label_format": "%.1f",
             "round_base": None,
             "constant_bounds": True,
             "num_table_values": 28,
             },
        }


def myround(x, base):
    """
    Round values to the nearest base
    """
    return int(base * round(float(x) / base))


def depth_slice(vtk, depth):
    """
    Slice a horizontal plane through a volume normal to the Z axis

    :type vtk: 
    :param vtk: opened data file
    :type depth: float
    :param depth: depth at which to slice through the model
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

def cross_section(vtk, normal, origin, name):
    """
    Cut a vertical cross section through a volume, with the plane parallel to
    the Z axis.

    :type vtk: 
    :param vtk: opened data file
    :type normal: list of float
    :param normal: normal vector of the plane
    :type origin: list of float
    :param origin: origin point for the 
    """
    slice_vtk = Slice(vtk)
    slice_vtk.SliceType = "Plane"
    slice_vtk.SliceType.Normal = normal
    slice_vtk.SliceType.Origin = origin
    
    # Apply the slice and render the new view
    renderView = GetActiveViewOrCreate("RenderView")
    sliceDisplay = Show(slice_vtk, renderView)
    RenameSource(name, slice_vtk)

    # Hide the plane widget
    Hide3DWidgets(proxy=slice_vtk)

    return slice_vtk


def create_ruler(point1, point2, label="", ticknum=5, color=None):
    """
    Generate a ruler to be used as a scalebar
    """
    renderView = GetActiveViewOrCreate("RenderView")

    ruler = Ruler()
    ruler.Point1 = point1
    ruler.Point2 = point2
    rulerDisplay = Show(ruler, renderView, "RulerSourceRepresentation")
    rulerDisplay.LabelFormat = label
    rulerDisplay.NumberOfTicks = ticknum
    rulerDisplay.AxisColor = color or COLOR
    rulerDisplay.Color = color or COLOR
    rulerDisplay.FontSize = FONTSIZE

    Hide3DWidgets(proxy=ruler)

    return ruler


def make_depth_slices(vtk, slices, save_path=os.getcwd()):
    """
    Main function for creating and screenshotting depth slices (slice plane 
    normal/perpendicular to Z axis). Creates a standardized look for the 
    depth slices with scale bar, colorbar and annotations
    """
    renderView = GetActiveViewOrCreate("RenderView")

    # Text needs to be deleted manually at the end
    text_objects = []

    # ANNOTATIONS: Depth slice annotations sit at the bottom-right corner
    text = Text()
    text_objects.append(text)
    textDisplay = Show(text, renderView)
    textDisplay.Position = [0.625, 0.1]
    textDisplay.FontSize = FONTSIZE * 2
    textDisplay.Color = COLOR

    # File name annotations sit at the top right corner
    text_tag = Text()
    text_objects.append(text_tag)
    text_tag.Text = save_fid
    textTagDisplay = Show(text_tag, renderView)
    textTagDisplay.Position = [0.55, 0.95]
    textTagDisplay.FontSize = FONTSIZE
    textTagDisplay.Color = COLOR

    Render()

    # RULER: Create rulers as scalebar, sits under the colorbar, mid left
    ruler_origin = [188068, 5722936, 0] 
    ruler_x = [288068, 5722936, 0]
    ruler_y = [188068, 5622936, 0]

    create_ruler(ruler_origin, ruler_x, "100km")
    create_ruler(ruler_origin, ruler_y)

    # CAMERA: Move camera a bit closer
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
            slice_vtk = depth_slice(vtk, slice_)
            tag = f"z_{slice_}km"

            # Determine colorbar bounds of slice to rescale colorbar
            vmin, vmax = slice_vtk.PointData.GetArray(0).GetRange(0)

        # Make sure that global bounds are set if required
        if preset["constant_bounds"]:
            vol_vmin, vol_vmax = vtk.PointData.GetArray(0).GetRange(0)
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
        SaveScreenshot(os.path.join(save_path, f"{tag}.png"), renderView,
                       ImageResolution=VIEW_SIZE, TransparentBackground=1)
        if slice_vtk:
            Hide(slice_vtk, renderView)
        else:
            Hide(vtk, renderView)

    # CLEAR: ResetSession() wasnt getting rid of the text so delete manually
    for txt in text_objects:
        Delete(txt)
        del txt
    ResetSession()


def make_trench_cross_sections(vtk, origins, save_path=os.getcwd()):
    """
    Create cross sections of a volume perpendicular to the Hikurangi trench 
    based on user-defined origin locations. Mark the origin locations 
    on the cross section for reference. Take screenshots
    """
    normal = [-0.64, -0.76, 0.]  # 40deg to the x-axis, from Donna's 2015 paper
    depth_cutoff_km = 100
    vol_vmin, vol_vmax = vtk.PointData.GetArray(0).GetRange(0)

    renderView = GetActiveViewOrCreate("RenderView")
    Hide(vtk, renderView)

    for i, (name, origin) in enumerate(origins.items()):
        ab = string.ascii_uppercase[i]
        tag = f"{ab}_{name.lower()}_trench_xsection"

        slice_vtk = cross_section(vtk=vtk, normal=normal, origin=origin,
                                  name=name)
        Hide(slice_vtk, renderView)

        # Cut off depths below a certain range as we are not interested in deep
        clip_vtk = Clip(slice_vtk)
        clip_vtk.ClipType = "Plane"
        clip_vtk.ClipType.Origin = [0., 0., abs(depth_cutoff_km) * -1E3]
        clip_vtk.ClipType.Normal = [0., 0., -1.]
        Show(clip_vtk, renderView, "UnstructuredGridRepresentation")

        # Reset the colorbounds to data range
        active = GetActiveSource()
        display = GetDisplayProperties(active, view=renderView)
        display.RescaleTransferFunctionToDataRange(False, True)

        # In order to set the camera, annotations, etc. in a general fashion,
        # we need to determine the corner grid locations of the slice
        dataPlane = servermanager.Fetch(clip_vtk)
        x, y, z = [], [], []
        for j in range(dataPlane.GetNumberOfPoints()):
            x_, y_, z_ = dataPlane.GetPoint(j)
            x.append(x_)
            y.append(y_)
            z.append(z_)

        # RULERS: Create rulers for use as scalebars
        dist_m = 100E3
        angle_rad = math.atan(normal[0] / normal[1])
        ruler_origin = [min(x), max(y), min(z)]  # bottom left

        # Horizontal ruler dimensions must be found with the power of geometry!
        # Took me way too long to figure out how to do this properly
        ruler_h = [min(x) + dist_m * math.cos(angle_rad),
                   max(y) - dist_m * math.sin(angle_rad),
                   min(z)]
        ruler_v = [min(x), max(y), min(z) + dist_m]

        rulerhDisplay = create_ruler(ruler_origin, ruler_h,
                                     f"{dist_m/1E3:.0f}km")
        rulervDisplay = create_ruler(ruler_origin, ruler_v)

        # GLPYH: Create a reference point based on the landmark location
        point = PointSource()
        point.Center = origin
        glyph = Glyph(Input=point, GlyphType="Cone")
        glyph.GlyphType.Direction = [0., 0., -1.]
        glyph.ScaleFactor = 10000
        Show(glyph, renderView, "GeometryRepresentation")

        # Annotate text to match glyph
        text = Text()
        text.Text = f"{ab}. {name}"
        textDisplay = Show(text, renderView)
        textDisplay.Position = [0.25, 0.65]
        textDisplay.FontSize = FONTSIZE * 2
        textDisplay.Color = COLOR

        # COLORMAP:
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
        cbar.Orientation = "Horizontal"
        cbar.TextPosition = "Ticks left/bottom, annotations right/top"
        cbar.Position = [0.25, 0.3]
        cbar.RangeLabelFormat = preset["range_label_format"]
        cbar.AddRangeLabels = 1
        cbar.ScalarBarThickness = 35
        cbar.ScalarBarLength = 0.15
        cbar.TitleFontSize = FONTSIZE
        cbar.LabelFontSize = FONTSIZE
        cbar.TitleColor = COLOR
        cbar.LabelColor = COLOR

        vmin, vmax = clip_vtk.PointData.GetArray(0).GetRange(0)

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

        vsLUT.NumberOfTableValues = preset["num_table_values"]
        vsLUT.RescaleTransferFunction(1.5, 2.1)
        # vsLUT.RescaleTransferFunction(vmin, vmax)

        # Ensure that the colorbar shows for each new slice
        display = GetDisplayProperties(clip_vtk, view=renderView)
        display.SetScalarBarVisibility(renderView, True)

        # Reset camera view to be normal to the plane. Specific to this plane
        renderView.CameraPosition = [-249942., 4837978., -57594.]
        renderView.CameraFocalPoint = [1145158., 6494661., -57594.]
        renderView.CameraViewUp = [0, 0, 1]


        SaveScreenshot(os.path.join(save_path, f"{tag}.png"), renderView,
                       ImageResolution=VIEW_SIZE, TransparentBackground=1)
        for object in [clip_vtk, rulerhDisplay, rulervDisplay, glyph]:
                       # textDisplay]:
            Hide(object, renderView)
        Delete(text)
        del text


def load_forest_state(fid):
    """
    Load a Paraview State file that may contain, e.g. source receiver
    glyphs, a coastline outline, etc. Specific to the Forest inversion which
    contains specific object names

    :type fid: str
    :param fid: file identifier for the state file
    """
    renderView = GetActiveViewOrCreate("RenderView")
    servermanager.LoadState(fid)

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


# =========================== USER PARAMETERS ==============================

# File ID's specifying the chosen files
state_fid = "/Users/Chow/Documents/academic/vuw/forest/forest.pvsm"

# Common path for finding and saving files
# base_path = "/Users/Chow/Documents/academic/vuw/forest/birch/"
base_path = "/Users/Chow/Documents/academic/vuw/forest/birch/"

# Allow for user-input files on the command line
data_fids = glob(os.path.join(base_path, "vtk", "model",
                              "ratio*0007*vs*"))

# Define where the slice the model in the XZ plane, surface implicit
slices = ["surface", 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 25, 30, 50]
slices = []

# Define trench parallel cross sections w/ origins based on landmark
# locations in UTM -60. Normal is defined as 40deg from the X-axis
origins = {"Wellington": [314007., 5426403., 0.,],
           "Flatpoint": [413092., 5433559., 0.,],
           "Castlepoint": [434724., 5471821., 0.,],
           "Akitio": [436980., 5511241., 0.],
           "Porangahau": [467051., 5538717., 0.,],
           "Elsthorpe": [484394., 5581561., 0.,],
           "Napier": [489374., 5626518., 0.,],
           "Mohaka": [507922., 5670909. ,0.,],
           "Mahia": [575567., 5665558., 0.,],
           "Gisborne": [588984., 5720001., 0.,],
           }

# Save path for screenshots
path_out = os.path.join(base_path, "test")
# ==========================================================================

for data_fid in data_fids:
    # Define the output directory to save figures
    save_fid = os.path.splitext(os.path.basename(data_fid))[0]
    save_path = os.path.join(path_out, save_fid)
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    # Presets defined by labels assigned to the file names, ignore model num
    preset = save_fid.split("_")
    preset = "_".join(preset[:1] + preset[2:])
    assert(preset in PRESETS), f"{preset} does not match preset keys"
    preset = PRESETS[preset]

    # Common viewing size so image proportions are correct in screenshots
    renderView = GetActiveViewOrCreate("RenderView")
    renderView.ViewSize = VIEW_SIZE

    # Open the model volume
    vtk = OpenDataFile(data_fid)
    RenameSource("surface", vtk)
    Show(vtk)
    ResetCamera()

    # Create depth slices
    if slices:
        load_forest_state(state_fid)
        make_depth_slices(vtk, slices, save_path=save_path)

    if origins:
        make_trench_cross_sections(vtk, origins, save_path=save_path)


