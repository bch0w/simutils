"""
A script to control Paraview using Python, setting up a figure and taking 
screenshots at various depths or cross-sections through the model.
Speficially designed for making screenshots of the NZNorth model during the 
Forest inversion

Must be run using PvPython from the CLI, which should be packaged with Paraview
"""
import os
import math
import string
import argparse
from glob import glob
from paraview.simple import (Slice, GetActiveViewOrCreate, RenameSource,
                             Hide3DWidgets, Ruler, Show, Text, Render,
                             GetActiveCamera, GetScalarBar,
                             GetColorTransferFunction, GetActiveView,
                             FindSource, SelectSurfacePoints,
                             ExtractSelection, Delete, SetActiveSource,
                             ClearSelection, GetActiveSource,
                             GetDisplayProperties, SaveScreenshot,
                             ResetSession, Hide, servermanager, OpenDataFile,
                             ResetCamera, Clip, PointSource, Glyph,
                             CreateRenderView, SetActiveView, GetSources)



# Constants
# An RGB color table to avoid looking up values, correspond to Python names
COLOR_TABLE = {"k": [0., 0., 0.], "w": [1., 1., 1.], "r": [1., 0., 0.],
               "b": [0., 0., 1.], "g": [0., 1., 0.], "y": [1., 1., 0.],
               "o": [1., .5, 0.], "c": [0., 1., 1.], "gray": [.5, .5, .5]}
# General fontsize and color for standard objects, e.g. text, rulers
FONTSIZE = 15
COLOR = COLOR_TABLE["k"]
# Size of the output screenshots
VIEW_SIZE = [1037, 813]

"""
Preset parameters that should stay constant among certain file types

:cbar_title (str): Label for the colorbar
:colormap (str): Colormap to define LUT
:invert_cmap (bool): Invert the original bounds of the colormap
:center_cmap (bool): Ensure that the middle value of the colormap is 0
:range_label_format (str): String formatter for the range bounds on colorbar
:round_base (int): Round range labels to a base value, if None, no rounding
:constant_bounds (bool or tuple): Keep the colorbar bounds constant for every
    screenshot. If 'true', uses the min/max values of the entire volume. If 
    tuple, then tuple must be in the form (min, max).
:num_table_values (int): Number of segmentations in the colorbar
"""
PRESETS = {
        "model_rho":
            {"cbar_title": "Density [kg m^-3]",
             "colormap": "Rainbow Desaturated",
             "invert_cmap": True,
             "center_cmap": False,
             "range_label_format": "%.1f",
             "round_base": 10,
             "constant_bounds": False,
             "num_table_values": 64,
             },
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
             "constant_bounds": (1.55, 1.9),
             "num_table_values": 28,
             },
        }


# Pre-defined trench parallel cross sections w/ origins based on landmark
# locations in UTM -60. Trench normal is defined as 40deg from the X-axis
TRENCH_XSECTIONS = {"Wellington": [314007., 5426403., 0.,],
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


def myround(x, base):
    """
    Round values to the nearest base
    """
    return int(base * round(float(x) / base))


def parse_slice_list(slice_list):
    """
    Slices are provided by the user and may consist of a list of mixed types, 

    1) "surface" indicates a surface projection and should be zeroth index
    2) depths may be provided as ints or floats, e.g. 5, 8.5. 
    3) depths can also be provided as ranges e.g. 10-15,1 would indicate a range
       of 10 to 15 by incriments of 1, inclusive, so 10, 11, 12, 13, 14, 15

    This function will parse this mixed list into an absolute list that can
    be iterated over by the main plotting function

    :type slice_list: list
    :param slice_list: list of input values that need to be parsed
    :rtype: list
    :return: parsed list that contains only integers and 'surface' (optional)
    """
    new_list = []
    include_surface = False

    for value in slice_list:
        if value == "surface":
            # Do this at the end so we can sort the list first
            include_surface = True
            continue
        else:
            try:
                new_list.append(int(value))
            except ValueError:
                try:
                    bounds, dx = value.split(",")
                    min_bound, max_bound = bounds.split("-")
                    for val in range(int(min_bound), 
                                     int(max_bound) + int(dx), int(dx)):
                        new_list.append(val)
                except Exception as e:
                    raise ValueError(f"{value} not in expected format: "
                                     "'surface', float, int or 'a-b,c'")

    # Remove any duplicates from the list
    new_list = sorted(list(set(new_list)))

    if include_surface:
        new_list.insert(0, "surface")

    return new_list


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
    renderView = GetActiveView()
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
    renderView = GetActiveView()
    Show(slice_vtk, renderView)
    RenameSource(name, slice_vtk)

    # Hide the plane widget
    Hide3DWidgets(proxy=slice_vtk)

    return slice_vtk


def create_ruler(point1, point2, label="", ticknum=5, axis_color=None,
                 font_color=None, reg_name="ruler"):
    """
    Generate a ruler to be used as a scalebar
    """
    renderView = GetActiveView()

    ruler = Ruler(registrationName=reg_name)
    ruler.Point1 = point1
    ruler.Point2 = point2
    rulerDisplay = Show(ruler, renderView, "RulerSourceRepresentation")
    rulerDisplay.LabelFormat = label
    rulerDisplay.NumberOfTicks = ticknum
    rulerDisplay.AxisColor = axis_color or COLOR
    rulerDisplay.Color = font_color or COLOR
    rulerDisplay.FontSize = FONTSIZE

    Hide3DWidgets(proxy=ruler)

    return ruler

def create_text(s, position, fontsize=None, color=None, bold=False,
                reg_name="text"):
    """
    Create and show a text object at a given position
    """
    renderView = GetActiveView()

    text = Text(registrationName=reg_name)
    text.Text = s
    textDisplay = Show(text, renderView)
    textDisplay.Position = position
    textDisplay.Bold = int(bold)
    textDisplay.FontSize = fontsize or FONTSIZE
    textDisplay.Color = color or COLOR

    return text, textDisplay


def create_xsection_data_axis_grid(source, min_depth_km=0, max_depth_km=100,
                                   dz_km=25, camera_parallel=False):
    """
    Create a standard looking axis grid for the trench cross sections which 
    puts horizontal lines at pre-determined depth levels. Unfornuately I havent
    found a way to do this for the vertical axis because we are plotting at some
    angle to the XY plane and Paraview doesnt have an easy way to shift the
    perspective
    """
    renderView = GetActiveView()
    renderView.CameraParallelProjection = int(camera_parallel)

    display = GetDisplayProperties(source, view=renderView)
    display.DataAxesGrid.GridAxesVisibility = 1
    display.DataAxesGrid.ShowGrid = 1

    display.DataAxesGrid.XAxisUseCustomLabels = 1
    display.DataAxesGrid.XAxisLabels = []

    display.DataAxesGrid.YAxisUseCustomLabels = 1
    display.DataAxesGrid.YAxisLabels = []

    display.DataAxesGrid.ZAxisUseCustomLabels = 1

    # Ensuring that all the depth values are formatted properly before ranging
    dz_m = -1 * int(abs(dz_km * 1E3))
    min_depth_m = int(abs(min_depth_km * 1E3))
    max_depth_m = -1 * int(abs(max_depth_km * 1E3)) + dz_m
    zrange = range(min_depth_m, max_depth_m, dz_m)

    display.DataAxesGrid.ZAxisLabels = list(zrange)
    display.DataAxesGrid.ZLabelFontSize = 1



def set_colormap_create_colorbar(vtk, position, orientation, colormap=None,
                                 title=None, fontsize=None, color=None):
    """
    Set the color transfer function based on preset values. Create a colorbar to
    match the colormap
    """
    # COLORMAP: Change the colormap to the desired preset value
    quantity = vtk.PointData.GetArray(0).Name  # e.g. model_init_vp
    vsLUT = GetColorTransferFunction(quantity)
    vsLUT.ApplyPreset(colormap or preset["colormap"], True)
    vsLUT.UseAboveRangeColor = 0
    vsLUT.UseBelowRangeColor = 0
    vsLUT.NumberOfTableValues = preset["num_table_values"]
    if preset["invert_cmap"]:
        vsLUT.InvertTransferFunction()

    # COLORBAR: Create the colorbar and set a common look
    cbar = GetScalarBar(vsLUT)
    cbar.Title = title or preset["cbar_title"]
    cbar.ComponentTitle = ""
    cbar.AutoOrient = 0
    cbar.Orientation = orientation
    cbar.Position = position
    cbar.RangeLabelFormat = preset["range_label_format"]
    cbar.AddRangeLabels = 1
    cbar.ScalarBarThickness = 35
    cbar.ScalarBarLength = 0.15
    cbar.TitleFontSize = fontsize or FONTSIZE
    cbar.LabelFontSize = fontsize or FONTSIZE
    cbar.TitleColor = color or COLOR
    cbar.LabelColor = color or COLOR

    return vsLUT, cbar


def rescale_colorscale(vsLUT, src, vtk, preset):
    """
    Rescale the colorbar based on user preference. Either centering colorscale
    on 0, setting global bounds, or rounding bounds to some base value

    :type vsLUT: ColorTransferFunction
    :param vsLUT: look up table that controls the colormap
    :type src: vtk object
    :param src: specific slice that is being rescaled, Slice, Clip, Extract
    :type vtk: vtk object
    :param vtk: The total volume used to find global min max values
    :type preset: dict:
    :param preset: the preset choices for how to deal with the colorscale/ map
    """
    # Set global bounds based on the min and max of the entire model
    # or set based on the given source file
    if preset["constant_bounds"]:
        if isinstance(preset["constant_bounds"], tuple):
            vmin, vmax = preset["constant_bounds"]
        else:
            vmin, vmax = vtk.PointData.GetArray(0).GetRange(0)
    else:
        vmin, vmax = src.PointData.GetArray(0).GetRange(0)

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

    return vsLUT


def get_coordinates(src):
    """
    Return the X, Y, Z coordinates of a given object
    :return:
    """
    dataPlane = servermanager.Fetch(src)
    x, y, z = [], [], []
    for j in range(dataPlane.GetNumberOfPoints()):
        x_, y_, z_ = dataPlane.GetPoint(j)
        x.append(x_)
        y.append(y_)
        z.append(z_)

    return x, y, z


def show_colorbar():
    """
    Sometimes colorbar is not set to visible, this function will force it out
    """
    renderView = GetActiveView()
    active = GetActiveSource()
    display = GetDisplayProperties(active, view=renderView)
    display.SetScalarBarVisibility(renderView, True)

def reset():
    """
    Delete any active sources (e.g. text, ruler) and reset the session,
    returning the pipeline to a blank slate
    """
    for x in GetSources().values():
        Delete(x[0])
    ResetSession()


def load_forest_state():
    """
    Load a Paraview State file that may contain, e.g. source receiver
    glyphs, a coastline outline, etc. Specific to the Forest inversion which
    contains specific object names
    """
    renderView = GetActiveView()
    servermanager.LoadState(
        "/Users/Chow/Documents/academic/vuw/forest/forest.pvsm"
    )

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


def make_depth_slices(fid, slices, preset, save_path=os.getcwd()):
    """
    Main function for creating and screenshotting depth slices (slice plane
    normal/perpendicular to Z axis). Creates a standardized look for the
    depth slices with scale bar, colorbar and annotations
    """
    # Open the model volume
    vtk = OpenDataFile(fid)
    RenameSource("surface", vtk)
    if "surface" in slices:
        Show(vtk)

    # ANNOTATIONS: Depth slice annotations sit at the bottom-right corner
    text, _ = create_text(s="", position=[0.625, 0.1], reg_name="text1",
                fontsize=FONTSIZE * 2)
    # Text showing the file name for easy id of data
    create_text(s=os.path.splitext(os.path.basename(data_fid))[0],
                position=[0.45, 0.95], reg_name="text2", fontsize=FONTSIZE * 2)

    # RULER: Create rulers as scalebar, sits under the colorbar, mid left
    ruler_origin = [198000, 5722936, 0]
    ruler_x = [298000, 5722936, 0]
    ruler_y = [198000, 5622936, 0]

    create_ruler(point1=ruler_origin, point2=ruler_x, label="100km [X]",
                 reg_name="ruler1")
    create_ruler(point1=ruler_origin, point2=ruler_y, label="[Y]",
                 reg_name="ruler2")

    # CAMERA: Hard set the camera as a top-down view over model, manually set
    ResetCamera()
    renderView = GetActiveView()
    renderView.CameraPosition = [401399., 5567959., 1276057.]
    renderView.CameraFocalPoint = [401399., 5567959., -197500.0]
    renderView.CameraViewUp = [0, 1, 0]
    renderView.OrientationAxesVisibility = 1
    Render()

    vsLUT, cbar = set_colormap_create_colorbar(vtk, position=[0.249, 0.78],
                                               orientation="Vertical",)

    # SLICE: Generate slices through the volume at desired depth levels
    for slice_ in slices:
        if slice_ == "surface":
            tag = "map"
            slice_vtk = vtk

            # 'Hacky' method to get data range by selecting surface points
            SelectSurfacePoints(Rectangle=[0, 0, renderView.ViewSize[0],
                                           renderView.ViewSize[1]])
            extract = ExtractSelection(Input=slice_vtk)
            rescale_colorscale(vsLUT, src=extract, vtk=vtk, preset=preset)

            # Get rid of the extraction surface
            Delete(extract)
            del extract
            SetActiveSource(slice_vtk)
            ClearSelection()
        else:
            slice_vtk = depth_slice(vtk, float(slice_))
            tag = f"z_{slice_}km"
            rescale_colorscale(vsLUT, src=slice_vtk, vtk=vtk, preset=preset)

        text.Text = tag
        show_colorbar()

        # SCREENSHOT: Save screenshot, hide slice and move on
        SaveScreenshot(os.path.join(save_path, f"{tag}.png"), renderView,
                       ImageResolution=VIEW_SIZE, TransparentBackground=1)
        Hide(slice_vtk, renderView)

    # CLEAR: Remove any additional objects created during plotting
    for reg_name in ["text1", "text2", "ruler1", "ruler2"]:
        src = FindSource(reg_name)
        Delete(src)


def make_cross_sections(fid, percentages, normal, preset, depth_cutoff_km=100,
                       save_path=os.getcwd()):
    """
    Create cross sections normal to the X or Y axes, create screenshots

    .. note::
        If you get an error with the min values of x, y or z, it means that the
        slice or clip operation has gone off the edge of the model

    :type fid: str
    :param fid: file id to plot
    :type sections: list of float
    :param sections: where to cut the cross section in the axis normal to the
        camera. given as percentages of the volume [0, 1].
    :type preset: dict
    :param preset: the preset dictionary that controls the colormap
    :type axis: str
    :param axis: the axis to plot, available: 'x', 'y'. If given as 'x', the
        x-axis of the plot will be the X axis of the model. the y-axis will
        always by the Z axis of the model.
    """
    normal = normal.lower()

    # Open the model volume
    vtk = OpenDataFile(fid)
    RenameSource("surface", vtk)
    Show(vtk)
    ResetCamera()

    # Get global min and max of array for colorscale before slicing
    renderView = GetActiveView()
    Hide(vtk, renderView)

    # In order to set the camera, annotations, etc. in a general fashion,
    # we need to determine the absolute extent of the volume
    x, y, z = get_coordinates(vtk)

    # Generally assign the horizontal axis based on the chosen normal
    if normal == "x":
        naxis = x
        nvector = [1., 0., 0.]
    elif normal == "y":
        naxis = y
        nvector = [0., 1., 0.]
    else:
        raise ValueError("normal must be 'x' or 'y'")

    for i, pct in enumerate(percentages):
        tag = f"{normal}_{pct}pct"

        # Determine what the actual distance is in relation to the given pct
        axis_dist = (pct * 1E-2) * (max(naxis) - min(naxis)) + min(naxis)
        origin = [_ * axis_dist for _ in nvector]
        slice_vtk = cross_section(vtk=vtk, normal=nvector,
                                  origin=origin, name=tag)
        Hide(slice_vtk, renderView)

        # Cut off depths below a certain range as we are not interested in deep
        clip_vtk = Clip(slice_vtk)
        clip_vtk.ClipType = "Plane"
        clip_vtk.ClipType.Origin = [0., 0., abs(depth_cutoff_km) * -1E3]
        clip_vtk.ClipType.Normal = [0., 0., -1.]
        Show(clip_vtk, renderView, "UnstructuredGridRepresentation")

        create_xsection_data_axis_grid(clip_vtk, camera_parallel=False)

        # Reset the colorbounds to data range
        active = GetActiveSource()
        display = GetDisplayProperties(active, view=renderView)
        display.RescaleTransferFunctionToDataRange(False, True)

        # RULERS: Create rulers for use as scalebars
        x, y, z = get_coordinates(clip_vtk)
        dist_m = 100E3
        ruler_origin = [min(x), max(y), min(z)]  # bottom left
        if normal == "x":
            ruler_h = [min(x), max(y) - dist_m, min(z)]
        elif normal == "y":
            ruler_h = [min(x) + dist_m, max(y), min(z)]
        ruler_v = [min(x), max(y), min(z) + dist_m]

        create_ruler(point1=ruler_origin, point2=ruler_h,
                     label=f"{dist_m / 1E3:.0f}km",
                     reg_name="ruler1")
        create_ruler(point1=ruler_origin, point2=ruler_v, label="[Z]",
                     reg_name="ruler2", font_color=COLOR_TABLE["w"])

        # Annotate distance
        create_text(f"{normal}={(axis_dist-min(naxis))*1E-3:.2f}km",
                    position=[0.25, 0.65], reg_name="text1")

        # Text showing the file name for easy id of data
        create_text(s=os.path.splitext(os.path.basename(data_fid))[0],
                    position=[0.55, 0.65], reg_name="text2",
                    fontsize=FONTSIZE * 2)

        # Generate and rescale the colorbar/ colormap
        vsLUT, cbar = set_colormap_create_colorbar(vtk, position=[0.25, 0.3],
                                                   orientation="Horizontal", )
        cbar.TextPosition = "Ticks left/bottom, annotations right/top"
        rescale_colorscale(vsLUT, src=clip_vtk, vtk=vtk, preset=preset)
        display = GetDisplayProperties(clip_vtk, view=renderView)
        display.SetScalarBarVisibility(renderView, True)

        # Reset camera view to be normal to the plane. Specific to this plane
        if normal == "x":
            renderView.CameraPosition = [-796319., 5599778., -75011.]
            renderView.CameraFocalPoint = [1347584., 5599778., -75011.]
            renderView.CameraViewUp = [0.0, 0.0, 1.0]
            renderView.CameraParallelScale = 554883.
        elif normal == "y":
            renderView.CameraPosition = [403183., 4578846., -62023.]
            renderView.CameraFocalPoint = [403183., 6971952., -62023.]
            renderView.CameraViewUp = [0.0, 0.0, 1.0]
            renderView.CameraParallelScale = 619381.
        renderView.OrientationAxesVisibility = 1

        Render()

        SaveScreenshot(os.path.join(save_path, f"{tag}.png"), renderView,
                       ImageResolution=VIEW_SIZE, TransparentBackground=1)

        # Clean up for next plot
        Hide(clip_vtk, renderView)
        for reg_name in ["ruler1", "ruler2", "text1", "text2"]:
            src = FindSource(reg_name)
            Delete(src)


def make_trench_cross_sections(fid, preset, depth_cutoff_km=100.,
                               save_path=os.getcwd()):
    """
    Create cross sections of a volume perpendicular to the Hikurangi trench
    based on user-defined origin locations. Mark the origin locations
    on the cross section for reference. Take screenshots

    :type depth_cutoff_km: float
    :param depth_cutoff_km: define where the bottom edge of the cross section
        will be. defaults to 100km depth
    """
    # Open the model volume
    vtk = OpenDataFile(fid)
    RenameSource("surface", vtk)
    Show(vtk)
    ResetCamera()

    # Pre-defined values for uniform look
    normal = [-0.64, -0.76, 0.]  # 40deg to the x-axis, from Donna's 2015 paper

    # Get global min and max of array for colorscale before slicing
    renderView = GetActiveView()
    Hide(vtk, renderView)

    for i, (name, origin) in enumerate(TRENCH_XSECTIONS.items()):
        ab = string.ascii_uppercase[i]  # alphabetize for easier identification
        tag = f"t_{ab.lower()}_{name.lower()}"

        slice_vtk = cross_section(vtk=vtk, normal=normal, origin=origin,
                                  name=name)
        Hide(slice_vtk, renderView)

        # Cut off depths below a certain range as we are not interested in deep
        clip_vtk = Clip(slice_vtk)
        clip_vtk.ClipType = "Plane"
        clip_vtk.ClipType.Origin = [0., 0., abs(depth_cutoff_km) * -1E3]
        clip_vtk.ClipType.Normal = [0., 0., -1.]
        Show(clip_vtk, renderView, "UnstructuredGridRepresentation")

        create_xsection_data_axis_grid(clip_vtk, camera_parallel=True)

        # Reset the colorbounds to data range
        active = GetActiveSource()
        display = GetDisplayProperties(active, view=renderView)
        display.RescaleTransferFunctionToDataRange(False, True)

        # In order to set the camera, annotations, etc. in a general fashion,
        # we need to determine the corner grid locations of the slice
        x, y, z = get_coordinates(clip_vtk)

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

        create_ruler(point1=ruler_origin, point2=ruler_h,
                     label=f"{dist_m/1E3:.0f}km", reg_name="ruler1")
        create_ruler(point1=ruler_origin, point2=ruler_v, label="[Z]",
                     reg_name="ruler2", font_color=COLOR_TABLE["w"])

        # GLPYH: Create a reference point based on the landmark location
        point = PointSource(registrationName="point1")
        point.Center = origin
        glyph = Glyph(Input=point, GlyphType="Cone", registrationName="glyph1")
        glyph.GlyphType.Direction = [0., 0., -1.]
        glyph.ScaleFactor = 10000
        Show(glyph, renderView, "GeometryRepresentation")

        # Annotate text to match glyph
        create_text(f"{ab}. {name}", [0.25, 0.65], reg_name="text1")

        # Text showing the file name for easy id of data
        create_text(s=os.path.splitext(os.path.basename(data_fid))[0],
                    position=[0.55, 0.65], reg_name="text2",
                    fontsize=FONTSIZE * 2)

        # Generate and rescale the colorbar/ colormap
        vsLUT, cbar = set_colormap_create_colorbar(vtk, position=[0.25, 0.3],
                                                   orientation="Horizontal",)
        cbar.TextPosition = "Ticks left/bottom, annotations right/top"
        rescale_colorscale(vsLUT, src=clip_vtk, vtk=vtk, preset=preset)
        display = GetDisplayProperties(clip_vtk, view=renderView)
        display.SetScalarBarVisibility(renderView, True)

        # Reset camera view to be normal to the plane. Specific to this plane
        renderView.CameraPosition = [-400553., 4689520., -49366.]
        renderView.CameraFocalPoint = [402390.0, 5595515., -49366.]
        renderView.CameraParallelScale = 300000.
        renderView.CameraViewUp = [0, 0, 1]
        renderView.OrientationAxesVisibility = 1
        Render()

        SaveScreenshot(os.path.join(save_path, f"{tag}.png"), renderView,
                       ImageResolution=VIEW_SIZE, TransparentBackground=1)

        # Clean up for next plot
        Hide(clip_vtk, renderView)
        for reg_name in ["ruler1", "ruler2", "text1", "point1",
                         "glyph1", "text2"]:
            src = FindSource(reg_name)
            Delete(src)


if __name__ == "__main__":
    # Command line arguments to define behavior of script
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="+", help="file to plot with paraview")
    parser.add_argument("-o", "--output", type=str, help="output path",
                        default=os.getcwd())
    parser.add_argument("-p", "--preset", type=str, 
                        help="preset colormap and labels given file type")
    parser.add_argument("-f", "--forest", action="store_true", 
                        help="include forest state file", default=False)
    parser.add_argument("-z", "--default_zslices", action="store_true", 
                        default=False, 
                        help="create default depth slices at the surface, "
                             "2-20km by increments of 2km, and 25-50km by "
                             "increments of 5km. Can be used in conjunciton "
                             "with -Z")
    parser.add_argument("-Z", "--zslices", nargs="+", 
                        help="manually set depths to slice at, 'surface' equals "
                             "surface projection. entries can be given as "
                             "ranges, e.g. 5-10,1 will produce 5,6,7,8,9,10; "
                             "depths and increments must be given in integers; "
                             "can be used in conjunction with -z",
                        default=[])
    parser.add_argument("-x", "--xslices", action="store_true", default=False,
                        help="create vertical cross sections with the viewing"
                             "plane normal to the X-axis.")
    parser.add_argument("-y", "--yslices", action="store_true", default=False,
                        help="create vertical cross sections with the viewing"
                             "plane normal to the Y-axis.")
    parser.add_argument("-t", "--trench", action="store_true",
                        help="plot pre-defined cross-sections normal to trench",
                        default=False)
    parser.add_argument("-c", "--cutoff_km", type=float, default=100,
                        help="For any vertical cross sections (Y-axis of figure "
                             "normal to Z axis of volume), define the depth"
                             "cutoff of the screenshot as usually were not "
                             "interested in looking at the entire volume. "
                             "Units of km, positive values only.")
    parser.add_argument("-v", "--verbose", action="store_true", default=False,
                        help="print output messages during the plotting")
    args = parser.parse_args()

    # Set up the active render view
    renderView = GetActiveViewOrCreate("RenderView")
    renderView.ViewSize = VIEW_SIZE
    SetActiveView(renderView)

    for data_fid in args.files:
        if args.verbose:
            print(f"Plotting file {data_fid}")

        assert(os.path.exists(data_fid)), f"{data_fid} not found"

        # Pre-define the output directory to save figures
        save_fid = os.path.splitext(os.path.basename(data_fid))[0]
        save_path = os.path.join(args.output, save_fid)

        if not os.path.exists(save_path):
            os.makedirs(save_path)

        # Presets defined by labels assigned to the file names or by user input
        if args.preset:
            preset_key = args.preset
        else:
            preset_key = save_fid.split("_")
            preset_key = "_".join(preset_key[:1] + preset_key[2:])
        assert(preset_key in PRESETS), (f"'{preset_key}' does not match "
                                        f"available presets {PRESETS.keys()}")
        if args.verbose:
            print(f"\tPreset is set to '{preset_key}'")
        preset = PRESETS[preset_key]

        # Create depth slices and generate screenshots
        zslices = []
        if args.default_zslices:
            zslices += ["surface", "2-20,2", "25-50,5"]
        if args.zslices:
            zslices += args.zslices

        if zslices:
            if args.forest:
                load_forest_state()

            zslices = parse_slice_list(zslices)
            if args.verbose:
                print(f"\tGenerating Z slices for {zslices}")
            make_depth_slices(data_fid, zslices, preset, save_path=save_path)
            reset()

        # Create cross sections normal to X axis and generate screenshots
        if args.xslices:
            # Do not go for 0 or 100 % because you might end up off the model
            xslices = list(range(5, 100, 10))
            if args.verbose:
                print(f"\tGenerating X-normal slices for {xslices}")
            make_cross_sections(data_fid, xslices, normal="x", preset=preset,
                                depth_cutoff_km=args.cutoff_km,
                                save_path=save_path)
            reset()

        # Create cross sections normal to Y axis and generate screenshots
        if args.yslices:
            # Possible to manually set the y-slice list here
            yslices = list(range(5, 100, 10))
            if args.verbose:
                print(f"\tGenerating Y-normal slices for {yslices}")
            make_cross_sections(data_fid, yslices, normal="y", preset=preset,
                                depth_cutoff_km=args.cutoff_km,
                                save_path=save_path)
            reset()

        # Create trench normal cross sections and generate screenshots
        if args.trench:
            if args.verbose:
                print("\tGenerating trench normal cross sections")
            make_trench_cross_sections(data_fid, preset,
                                       depth_cutoff_km=args.cutoff_km,
                                       save_path=save_path)
            reset()
           

