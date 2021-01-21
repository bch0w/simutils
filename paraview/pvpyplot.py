"""
A script to control Paraview using Python, setting up a figure and taking 
screenshots at various depths or cross-sections through the model.
Speficially designed for making screenshots of the NZNorth model during the 
Forest inversion

Must be run using PvPython from the CLI, which should be packaged with Paraview

.. rubric::
    Running this from the command line looks something like this:

    pvpython VTK_FILE.vtk -p model_vs -c -xyt -Z surface 10-50,1 -b 1000,3000

    -p model_vs: use the model_vs preset to define the colorbar
    -c: plot the coastline ontop of the model
    -xyt: create default cross sections normal to x-axis, y-axis and trench
    -Z ... : create depth slices for map view (surface) and for depths of 10 to
        50 km by steps of 1
    -b 1000,3000: hard-set the colorbar bounds from 1000 to 3000 m/s

"""
import os
import math
import string
import argparse
from paraview.vtk.numpy_interface import dataset_adapter
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
                             CreateRenderView, SetActiveView, GetSources, Axes,
                             ColorBy, PointDatasetInterpolator)



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
:bounds (bool or tuple or list): Can be overwritten using the '-b' arg
    * True: keep the colorbar bounds constant for every screenshot using the 
    min/max values of the entire volume. 
    * False: use whatever default bounds are set by Paraview, 
    * tuple/list: Manually set the bounds; must be in the form (min, max)
:num_table_values (int): Number of segmentations in the colorbar
"""
PRESETS = {
        "custom_vs":
            {"cbar_title": "Vs [m/s]",
             "colormap": "RdYlBu",
             "invert_cmap": False,
             "center_cmap": False,
             "range_label_format": "%.1f",
             "round_base": 10,
             "bounds": False,
             "num_table_values": 21,
             },
        "custom_vpvs":
            {"cbar_title": "Vp/Vs Ratio",
             "colormap": "Yellow - Gray - Blue",
             "invert_cmap": False,
             "center_cmap": False,
             "range_label_format": "%.2f",
             "round_base": None,
             "bounds": (1.55, 2.1),
             "num_table_values": 28,
             },
        "trench_vs":
            {"cbar_title": "Vs [m/s]",
             "colormap": "Rainbow Desaturated",
             "invert_cmap": True,
             "center_cmap": False,
             "range_label_format": "%.1f",
             "round_base": 10,
             "bounds": [2000, 5500],
             "num_table_values": 64,
             },
        "model_rho":
            {"cbar_title": "Density [kg m^-3]",
             "colormap": "Rainbow Desaturated",
             "invert_cmap": True,
             "center_cmap": False,
             "range_label_format": "%.1f",
             "round_base": 10,
             "bounds": False,
             "num_table_values": 64,
             },
        "model_vp":
            {"cbar_title": "Vp [m/s]",
             "colormap": "Rainbow Desaturated",
             "invert_cmap": True,
             "center_cmap": False,
             "range_label_format": "%.1f",
             "round_base": 10,
             "bounds": False,
             "num_table_values": 64,
             },
        "model_vs":
            {"cbar_title": "Vs [m/s]",
             "colormap": "Rainbow Desaturated",
             "invert_cmap": True,
             "center_cmap": False,
             "range_label_format": "%.1f",
             "round_base": 10,
             "bounds": False,
             "num_table_values": 64,
             },
        "gradient_vp_kernel":
            {"cbar_title": "Vp Gradient [m^-2 s^2]",
             "colormap": "Cool to Warm (Extended)",
             "invert_cmap": True,
             "center_cmap": True,
             "range_label_format": "%.1E",
             "round_base": None,
             "bounds": True,
             "num_table_values": 64,
             },
        "gradient_vs_kernel":
            {"cbar_title": "Vs Gradient [m^-2 s^2]",
             "colormap": "Cool to Warm (Extended)",
             "invert_cmap": True,
             "center_cmap": True,
             "range_label_format": "%.1E",
             "round_base": None,
             "bounds": True,
             "num_table_values": 64,
             },
        "update_vp":
            {"cbar_title": "Vp Net Update [ln(m/m00)]",
             "colormap": "Blue Orange (divergent)",
             "invert_cmap": True,
             "center_cmap": True,
             "range_label_format": "%.02f",
             "round_base": None,
             "bounds": True,
             "num_table_values": 64,
             },
        "update_vs":
            {"cbar_title": "Vs Net Update [ln(m/m00)]",
             "colormap": "Blue Orange (divergent)",
             "invert_cmap": True,
             "center_cmap": True,
             "range_label_format": "%.02f",
             "round_base": None,
             "bounds": [-.15, .15],
             # "bounds": False,
             "num_table_values": 64,
             },
        "ratio_poissons":
            {"cbar_title": "Poisson's Ratio",
             "colormap": "Blue - Green - Orange",
             "invert_cmap": False,
             "center_cmap": False,
             "range_label_format": "%.1f",
             "round_base": None,
             "bounds": True,
             "num_table_values": 64,
             },
        "ratio_vpvs":
            {"cbar_title": "Vp/Vs Ratio",
             "colormap": "Cool to Warm (Extended)",
             "invert_cmap": False,
             "center_cmap": False,
             "range_label_format": "%.2f",
             "round_base": None,
             "bounds": (1.55, 1.9),
             "num_table_values": 28,
             },
        }


# Pre-defined trench parallel cross sections w/ origins based on landmark
# locations in UTM -60. Trench normal is defined as 40deg from the X-axis
# This normal vector defines 40deg to the x-axis (40 deg from Donna 2015)
TRENCH_NORMAL = [-0.64, -0.76, 0.]
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
        if value in ["surface", "s"]:
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
    Show(slice_vtk, renderView)
    RenameSource(f"z_{depth}km", slice_vtk)
    Hide3DWidgets(proxy=slice_vtk)

    return slice_vtk


def cross_section(vtk, normal, origin, name):
    """
    Cut a vertical cross section through a volume, with the plane parallel to
    the Z axis.

    :type vtk: paraview.servermanager.LegacyVTKReader
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
    Hide3DWidgets(proxy=slice_vtk)

    return slice_vtk


def create_ruler(point1, point2, label="", ticknum=5, axis_color=None,
                 font_color=None, reg_name="ruler"):
    """
    Generate a ruler to be used as a scalebar

    :type point1: list of float
    :param point1: [x0, y0, z0] starting point of the ruler
    :type point 2: list of float
    :param point2: [x, y, z] ending point of the ruler
    :type label: str
    :param label: text label that will be annotated at the midpoint of the ruler
    :type ticknum: int
    :param ticknum: number of individual ticks to place on the ruler
    :type axis_color: str
    :param axis_color: color to make the ruler, in python color values, e.g. 'k'
    :type font_color: str
    :param axis_color: color to make the label, in python color values, e.g. 'k'
    :type reg_name: str
    :param reg_name: Name to register the ruler, useful for deleting the
        ruler afterwards using e.g. FindSource
    :rtype: Ruler
    :return: the created ruler object
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

    return reg_name

def create_text(s, position, fontsize=None, color=None, bold=False,
                reg_name="text"):
    """
    Create and show a text object at a given position

    :type position: list of float
    :param position: [x, y] where x and y are in percentage of the viewing
        window, e.g. [0.5, 0.5] places text in the middle of the screen
    :type fontsize: float
    :param fontsize: fontsize of the text to create, defaults to the fontsize
        constant defined at the top of the script
    :type color: list of float
    :param color: RGB color of the text, defaults to the constant color defined
        at top of the file
    :type bold: bool
    :param bold: if True, boldfaces the text
    :type reg_name: str
    :param reg_name: Name to register the ruler, useful for deleting the
        ruler afterwards using e.g. FindSource
    :rtype: tuple
    :return: (text object created, rendered view of the text object)
    """
    renderView = GetActiveView()

    text = Text(registrationName=reg_name)
    text.Text = s
    textDisplay = Show(text, renderView)
    textDisplay.Position = position
    textDisplay.Bold = int(bold)
    textDisplay.FontSize = fontsize or FONTSIZE
    textDisplay.Color = color or COLOR

    return text, reg_name


def create_ruler_grid_axes_depth_slice(src, tick_spacing_km=50, top=True,
                                       bottom=True, left=True, right=True):
    """
    Grid axes are a bit finagly because they can be difficult to control.
    Here we use rulers to outline the source domain, which allows us to
    fine tune the tick spacing, location etc., as well as which axes we want to
    plot on

    Works by generating rulers from minimum to maximum coordinates

    :type source: paraview.servermanager.Slice etc...
    :param src: the object we want to make grids for, can be the VTK volume,
        or resulting slice
    :type tick_spacing_km: int
    :param tick_spacing_km: the spacing of the ruler ticks in units of km
    :type top/bottom/left/right: bool
    :param top/bottom/left/right: turn on/off each of the corresponding 'axes'
    :rtype: list of str
    :return: the reg names for each of the generated rulers so they are easier
        to delete
    """
    tick_spacing_m = tick_spacing_km * 1E3
    x, y, z = get_coordinates(src)

    # Figure out an acceptable end point that allows us to have even grid
    # spacings but is as close to the actual volume end point as possible
    x_end_point = myround(max(x) - min(x), tick_spacing_m) + min(x)
    while x_end_point > max(x):
        x_end_point -= tick_spacing_m
    ruler_x_bot = [x_end_point, min(y), max(z)]
    ruler_x_top = [x_end_point, max(y), max(z)]
    tick_num_x = int((x_end_point - min(x)) // tick_spacing_m)

    y_end_point = myround(max(y) - min(y), tick_spacing_m) + min(y)
    while y_end_point > max(y):
        y_end_point -= tick_spacing_m
    ruler_y_lft = [min(x), y_end_point, max(z)]
    ruler_y_rgt = [max(x), y_end_point, max(z)]
    tick_num_y = int((y_end_point - min(y)) // tick_spacing_m)

    reg_names = []
    if bottom:
        create_ruler(point1=ruler_x_bot, point2= [min(x), min(y), max(z)],
                     label=f"[X] (dx={tick_spacing_km}km)",
                     reg_name="ruler_bottom", ticknum=tick_num_x)
        reg_names.append("ruler_bottom")

    if top:
        create_ruler(point1=[min(x), max(y), max(z)], point2=ruler_x_top,
                     reg_name="ruler_top", ticknum=tick_num_x)
        reg_names.append("ruler_top")

    if left:
        create_ruler(point1=[min(x), min(y), max(z)], point2=ruler_y_lft,
                     label=f"[Y]\n(dy={tick_spacing_km}km)",
                     reg_name="ruler_left", ticknum=tick_num_y)
        reg_names.append("ruler_left")
    if right:
        create_ruler(point1=ruler_y_rgt, point2=[max(x), min(y), max(z)],
                     reg_name="ruler_right", ticknum=tick_num_y)
        reg_names.append("ruler_right")

    return reg_names


def create_minmax_glyphs(src, glyph_type="2D Glyph", glyph_type_2d="Cross",
                        scale_factor=15000):
    """
    Create a glyph that denotes the position of the minimum and maximum values
    for a given source (usually a slice).

    :type src: paraview.servermanager.LegacyVTKReader or whatever
    :param src: source object to query for min max values
    :type glyph_type: str
    :param glyph_type: the glyph that you want to use to mark the values
    :type glyph_type_2d: str
    :param glyph_type_2d: if you choose '2D Glyph' as your glyph type, further
        choice about which 2D glyph
    :type scale_factor: int
    :param scale_factor: size of the corresponding glyphs
    """
    dataPlane = servermanager.Fetch(src)
    dataPlane = dataset_adapter.WrapDataObject(dataPlane)
    data_points = list(dataPlane.PointData[0])
    min_val_idx = data_points.index(min(data_points))
    max_val_idx = data_points.index(max(data_points))

    reg_names = []
    for i, idx in enumerate([min_val_idx, max_val_idx]):
        x, y, z = dataPlane.GetPoint(idx)

        # GLPYH: Create a reference point based on the location of the min/max
        point = PointSource(registrationName=f"point_minmax_{i}")
        point.Center = [x, y, z]
        glyph = Glyph(Input=point, GlyphType=glyph_type,
                      registrationName=f"glyph_minmax_{i}")
        if glyph_type == "2D Glyph":
            glyph.GlyphType.GlyphType = glyph_type_2d
        glyph.ScaleFactor = scale_factor
        Show(glyph, renderView, "GeometryRepresentation")

        reg_names.append(f"point_minmax_{i}")
        reg_names.append(f"glyph_minmax_{i}")

    return reg_names


def set_xsection_data_axis_grid(source, min_depth_km=0, max_depth_km=100,
                                dz_km=25, camera_parallel=False):
    """
    Create a standard looking axis grid for the trench cross sections which
    puts horizontal lines at pre-determined depths across the entire slice.

    :type source: paraview.servermanager.Slice
    :param source: the given slice that we want to make an axis grid for
    :type min_depth_km: int
    :param min_depth_km: the top (+Z) of the slice, used to define where to
        start generating labels
    :type max_depth_km: int
    :param max_depth_km: the bottom (-Z) of the slice, used to define where to
        stop generating labels, not inclusive
    :type dz_km: int
    :param dz_km: spacing for the labels
    :type camera_parallel: bool
    :param camera_parallel: turn on the parallel camera projection option,
        useful for when your camera is not facing an orthogonal direction (x, y
        or z), to keep the axis grid flat against the slice.
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


def set_colormap_colorbar(vtk, position, orientation, colormap=None,
                          title=None, fontsize=None, color=None):
    """
    Set the color transfer function based on preset values.
    Then, create a colorbar to match the colormap

    :type vtk: paraview.servermanager.LegacyVTKReader
    :param vtk: opened data file
    :type position: list of float
    :param position: [x, y] where x and y are in percentage of the viewing
        window, e.g. [0.5, 0.5] places text in the middle of the screen
    :type orientation: str
    :param orientation: colorbar orientation, 'Horizontal' or 'Vertical'
    :type colormap: str
    :param colormap: the name of the colormap LUT defined in Paraview. If None,
        defaults to the preset value defined by the user or filename
    :type title: str
    :param title: corresponding label for the colorbar
    :type fontsize: int
    :param fontsize: fontsize of the title, defaults to FONTSIZE
    :type color: str
    :param color: color of the title, defaults to COLOR
    """
    # Change the colormap to the desired preset value
    quantity = vtk.PointData.GetArray(0).Name  # e.g. model_init_vp
    vsLUT = GetColorTransferFunction(quantity)
    vsLUT.ApplyPreset(colormap or preset["colormap"], True)
    vsLUT.UseAboveRangeColor = 0
    vsLUT.UseBelowRangeColor = 0
    vsLUT.NumberOfTableValues = preset["num_table_values"]
    if preset["invert_cmap"]:
        vsLUT.InvertTransferFunction()

    # Create the colorbar and set a common look
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
    on 0, setting global bounds, or rounding bounds to some base value.

    :type vsLUT: ColorTransferFunction
    :param vsLUT: look up table that controls the colormap
    :type src: vtk object
    :param src: paraview.servermanager.'OBJECT', slice or vtk reader etc...
    :type vtk: paraview.servermanager.LegacyVTKReader
    :param vtk: The open volume used to find global min max values
    :type preset: dict
    :param preset: the preset choices for how to deal with the colorscale/ map
    """
    # Set global bounds based on the min and max of the entire model
    # or set based on the given source file
    if preset["bounds"]:
        if isinstance(preset["bounds"], (tuple, list)):
            vmin, vmax = preset["bounds"]
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
   
    # Manually set the values for the colorbar to include bounds and midpoint
    cbar = GetScalarBar(vsLUT)
    cbar.UseCustomLabels = 1
    cbar.CustomLabels = [vmin, (vmax + vmin) / 2, vmax]

    return vsLUT


def show_colorbar(src=None):
    """
    Sometimes colorbar is not set to visible, this function will make it visible
    """
    renderView = GetActiveView()
    if src is None:
        src = GetActiveSource()
    display = GetDisplayProperties(src, view=renderView)
    display.SetScalarBarVisibility(renderView, True)


def reset():
    """
    Delete any active sources (e.g. text, ruler) and reset the session,
    returning the pipeline to a blank slate. Not garuanteed to get rid of
    all objects

    :type reg_names: list
    :param reg_names: names of registered objects that need to be deleted
    """
    for src in [_[0] for _ in GetSources().keys()]:
        src = FindSource(src)
        Delete(src)
    # for x in GetSources().values():
    #     Delete(x[0])
    ResetSession()


def delete_temp_objects(ruler=False, point=False, glyph=False, text=False):
    """
    Delete objects that are made specifically for each screenshot, e.g. points,
    glyphs, rulers, text

    :type reg_names: list
    :param reg_names: names of registered objects that need to be deleted
    """
    for src in [_[0] for _ in GetSources().keys()]:
        if (ruler and "ruler" in src) or (point and "point" in src) or \
            (glyph and "glyph" in src) or (text and "text" in src):
            src = FindSource(src)
            Delete(src)


def get_coordinates(src):
    """
    Return the X, Y, Z coordinates of a given object

    :type: paraview.servermanager.LegacyVTKReader or whatever
    :param src: the object that you want to query
    :rtype: tuple of lists
    :return: the X, Y, Z coordinates of the src
    """
    dataPlane = servermanager.Fetch(src)
    x, y, z = [], [], []
    for j in range(dataPlane.GetNumberOfPoints()):
        x_, y_, z_ = dataPlane.GetPoint(j)
        x.append(x_)
        y.append(y_)
        z.append(z_)

    return x, y, z


def plot_coastline(point_size=2., color=None):
    """
    Plot a coastline file from an existing .VTK file

    :type point_size: float
    :param float_size: size of the points representing the coastline
    :type color: str
    :param color: color of the coastline, defaults to COLOR
    """
    if color is None:
        color = COLOR_TABLE["k"]
    renderView = GetActiveView()
    coast = OpenDataFile("/Users/Chow/Documents/academic/vuw/forest/utils/"
                         "vtk_files/coast.vtk")
    RenameSource("coast", coast)
    Show(coast, renderView)
    coastDisplay = GetDisplayProperties(coast, view=renderView)
    coastDisplay.SetScalarBarVisibility(renderView, False)
    coastDisplay.PointSize = point_size
    coastDisplay.AmbientColor = color or COLOR
    coastDisplay.DiffuseColor = color or COLOR
    ColorBy(coastDisplay, None)


def plot_events():
    """
    Plot FOREST INVERSION event locations (VTK file) as 2D circles normal to
    the Z axis
    """
    renderView = GetActiveView()
    srcs = OpenDataFile("/Users/Chow/Documents/academic/vuw/forest/utils/"
                         "vtk_files/srcs_d.vtk")
    RenameSource("srcs", srcs)
    glyph = Glyph(registrationName="srcs_g", Input=srcs,
                  GlyphType="2D Glyph")
    glyph.ScaleArray = ["POINTS", "No scale array"]
    glyph.ScaleFactor = 10000.
    glyph.GlyphMode = "All Points"
    glyph.GlyphType.GlyphType = "Circle"
    glyph.GlyphType.Filled = 1
    glyphDisplay = Show(glyph, renderView, "GeometryRepresentation")
    glyphDisplay.PointSize = 3.
    glyphDisplay.SetRepresentationType("Surface With Edges")
    glyphDisplay.LineWidth = 2.0
    glyphDisplay.Opacity = 1.
    vsLUT = GetColorTransferFunction("Z_Value")
    vsLUT.ApplyPreset("Inferno (matplotlib)")


def plot_stations():
    """
    Plot FOREST inversion station locations (VTK file) as 2D diamonds normal
    to the Z axis
    """
    renderView = GetActiveView()
    rcvs = OpenDataFile("/Users/Chow/Documents/academic/vuw/forest/utils/"
                         "vtk_files/rcvs.vtk")
    RenameSource("rcvs", rcvs)
    glyph = Glyph(registrationName="rcvs_g", Input=rcvs,
                  GlyphType="2D Glyph")
    glyph.ScaleArray = ["POINTS", "No scale array"]
    glyph.ScaleFactor = 10000.
    glyph.GlyphMode = "All Points"
    glyph.GlyphType.GlyphType = "Diamond"
    glyph.GlyphType.Filled = 1
    glyphDisplay = Show(glyph, renderView, "GeometryRepresentation")
    glyphDisplay.PointSize = 3.
    glyphDisplay.SetRepresentationType("Surface With Edges")
    glyphDisplay.LineWidth = 2.0
    glyphDisplay.Opacity = 1.


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

    # Annotate the depth at the bottom-right corner
    text, _ = create_text(s="", position=[0.65, 0.05], reg_name="text1",
                          fontsize=FONTSIZE * 2)

    # Annotate the file id at the top-left corner
    create_text(s=os.path.splitext(os.path.basename(data_fid))[0],
                position=[0.249, 0.9], reg_name="text2", fontsize=FONTSIZE * 2)

    # Create bounding axes with pre-defined tick marks using rulers
    create_ruler_grid_axes_depth_slice(vtk)

    # CAMERA: Hard set the camera as a top-down view over model, manually set
    ResetCamera()
    renderView = GetActiveView()
    renderView.InteractionMode = "2D"
    renderView.CameraPosition = [403648., 5595153., 1721019.]
    renderView.CameraFocalPoint = [403648., 5595153., 2300.0]
    renderView.CameraParallelScale = 326545.
    renderView.CameraViewUp = [0, 1, 0]
    Render()

    # This is good for a 3D projection but that's trash for viewing 2D planes
    # renderView.CameraPosition = [401399., 5567959., 1276057.]
    # renderView.CameraFocalPoint = [401399., 5567959., -197500.0]


    # Sit the colorbar partway down from the top left corner
    vsLUT, cbar = set_colormap_colorbar(vtk, position=[0.25, 0.7],
                                        orientation="Vertical",)

    # Generate slices through the volume at desired depth levels
    # Make special precautions if we're looking at surface projections
    for slice_ in slices:
        if slice_ == "surface":
            tag = "z_00map"
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
            create_minmax_glyphs(slice_vtk)
            tag = f"z_{slice_:0>2}km"
            rescale_colorscale(vsLUT, src=slice_vtk, vtk=vtk, preset=preset)

        text.Text = tag
        show_colorbar(slice_vtk)
        create_minmax_glyphs(slice_vtk)

        # Save screenshot, hide slice and move on, job done
        SaveScreenshot(os.path.join(save_path, f"{tag}.png"), renderView,
                       ImageResolution=VIEW_SIZE, TransparentBackground=1)
        Hide(slice_vtk, renderView)
        
        # Delete the min max value points because they'll change w/ each slice
        delete_temp_objects(glyph=True)

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

    vtk = OpenDataFile(fid)
    RenameSource("surface", vtk)
    Show(vtk)
    ResetCamera()
    renderView = GetActiveView()
    renderView.InteractionMode = "2D"
    Hide(vtk, renderView)

    # In order to set the camera, annotations, etc. in a general fashion,
    # we need to determine the absolute extent of the volume
    x, y, z = get_coordinates(vtk)
    if normal == "x":
        axis = y
        naxis = x
        parallel = "y"
        nvector = [1., 0., 0.]
    elif normal == "y":
        axis = x
        naxis = y
        parallel = "x"
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

        set_xsection_data_axis_grid(clip_vtk, camera_parallel=False)

        # Reset the colorbounds to data range
        active = GetActiveSource()
        display = GetDisplayProperties(active, view=renderView)
        display.RescaleTransferFunctionToDataRange(False, True)

        # Create rulers for use as scalebars
        tick_spacing_m = 50E3
        x, y, z = get_coordinates(clip_vtk)
        ruler_origin = [min(x), max(y), min(z)]  # bottom left
        dist_m_h = myround(max(axis) - min(axis), tick_spacing_m)
        while dist_m_h > (max(axis) - min(axis)):
            dist_m_h -= tick_spacing_m
        num_ticks_h = int(dist_m_h // tick_spacing_m)

        # Determine how long the ruler needs to be for each axis
        if normal == "x":
            ruler_h = [min(x), max(y) - dist_m_h, min(z)]
        elif normal == "y":
            ruler_h = [min(x) + dist_m_h, max(y), min(z)]

        # Vertical axis is based off the cutoff depth defined by the user
        dist_m_v = depth_cutoff_km * 1E3
        ruler_v = [min(x), max(y), min(z) + dist_m_v]

        create_ruler(point1=ruler_origin, point2=ruler_h,
                     label=f"[{parallel.upper()}] "
                           f"(d{parallel}={int(tick_spacing_m*1E-3)}km)",
                     reg_name="ruler1", ticknum=num_ticks_h)
        create_ruler(point1=ruler_v, point2=ruler_origin,
                     label="[Z]\n(dz=25km)", reg_name="ruler2",)

        # Annotate distance
        create_text(f"{normal}={(axis_dist-min(naxis))*1E-3:.2f}km",
                    position=[0.2, 0.25], reg_name="text1",
                    fontsize=int(FONTSIZE * 1.5))

        # Text showing the file name for easy id of data
        create_text(s=os.path.splitext(os.path.basename(data_fid))[0],
                    position=[0.6, 0.25], reg_name="text2",
                    fontsize=int(FONTSIZE * 1.5))

        # Generate and rescale the colorbar/ colormap
        vsLUT, cbar = set_colormap_colorbar(vtk, position=[0.4, 0.275],
                                            orientation="Horizontal", )
        cbar.TextPosition = "Ticks left/bottom, annotations right/top"
        rescale_colorscale(vsLUT, src=clip_vtk, vtk=vtk, preset=preset)
        display = GetDisplayProperties(clip_vtk, view=renderView)
        display.SetScalarBarVisibility(renderView, True)

        # Reset camera view to be normal to the plane. Specific to this plane
        renderView.CameraViewUp = [0.0, 0.0, 1.0]
        if normal == "x":
            renderView.CameraPosition = [-1100117., 5629267., -39303.]
            renderView.CameraFocalPoint = [402390.0, 5629267., -39303.]
            renderView.CameraParallelScale = 274961.
        elif normal == "y":
            renderView.CameraPosition = [378617., 4093007., -47332.]
            renderView.CameraFocalPoint = [378617., 5595515.0, -47332.]
            renderView.CameraParallelScale = 205602.
        Render()

        SaveScreenshot(os.path.join(save_path, f"{tag}.png"), renderView,
                       ImageResolution=VIEW_SIZE, TransparentBackground=1)

        # Clean up for next plot
        Hide(clip_vtk, renderView)
        delete_temp_objects(ruler=True, text=True)



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
    vtk = OpenDataFile(fid)
    RenameSource("surface", vtk)
    Show(vtk)
    ResetCamera()
    renderView = GetActiveView()
    renderView.InteractionMode = "2D"
    Hide(vtk, renderView)

    for i, (name, origin) in enumerate(TRENCH_XSECTIONS.items()):
        # Alphabetize the cross sections for easier identification
        ab = string.ascii_uppercase[i]
        tag = f"t_{ab.lower()}_{name.lower()}"

        # Slice across the given cross section plane
        slice_vtk = cross_section(vtk=vtk, normal=TRENCH_NORMAL,
                                  origin=origin, name=name)
        Hide(slice_vtk, renderView)

        # Cut off depths below a certain range as we are not interested in deep
        clip_vtk = Clip(slice_vtk)
        clip_vtk.ClipType = "Plane"
        clip_vtk.ClipType.Origin = [0., 0., abs(depth_cutoff_km) * -1E3]
        clip_vtk.ClipType.Normal = [0., 0., -1.]
        Show(clip_vtk, renderView, "UnstructuredGridRepresentation")

        # Put markers for depth values
        set_xsection_data_axis_grid(clip_vtk, camera_parallel=True)

        # Reset the colorbounds to data range
        active = GetActiveSource()
        display = GetDisplayProperties(active, view=renderView)
        display.RescaleTransferFunctionToDataRange(False, True)

        # In order to set the camera, annotations, etc. in a general fashion,
        # we need to determine the corner grid locations of the slice
        x, y, z = get_coordinates(clip_vtk)
        ruler_origin = [min(x), max(y), min(z)]  # bottom left

        # Find the optimal length of the slice to fit an even spacing of ticks
        tick_spacing_m = 50E3
        slice_length_m = ((max(y) - min(y)) ** 2 + (max(x) - min(x)) ** 2) ** .5
        dist_m_h = myround(slice_length_m, tick_spacing_m)
        # Sometimes the round overestimates so we just go back until its not
        while dist_m_h > slice_length_m:
            dist_m_h -= tick_spacing_m
        num_ticks_h = int(dist_m_h // tick_spacing_m)

        # Horizontal ruler dimensions must be found with the power of trig.
        # (it took me way too long to figure out how to do this properly, oh
        #  highschool trig...)
        angle_rad = math.atan(TRENCH_NORMAL[0] / TRENCH_NORMAL[1])
        ruler_h = [min(x) + dist_m_h * math.cos(angle_rad),
                   max(y) - dist_m_h * math.sin(angle_rad),
                   min(z)]

        # Keep the vertical scaling constant because it won't change
        dist_m_v = depth_cutoff_km * 1E3
        ruler_v = [min(x), max(y), min(z) + dist_m_v]

        # Generate appropriate rulers that act as the X and Y axes in this plane
        create_ruler(point1=ruler_origin, point2=ruler_h,
                     label=f"{dist_m_h/1E3:.0f}km " 
                           f"(dh={int(tick_spacing_m*1E-3)}km)",
                     reg_name="ruler1", ticknum=num_ticks_h)

        create_ruler(point1=ruler_v, point2=ruler_origin,
                     label="[Z]\n(dz=25km)", reg_name="ruler2")

        # Create a reference point based on the landmark location
        point = PointSource(registrationName="point1")
        point.Center = origin
        glyph = Glyph(Input=point, GlyphType="Cone", registrationName="glyph1")
        glyph.GlyphType.Direction = [0., 0., -1.]
        glyph.ScaleFactor = 10000
        Show(glyph, renderView, "GeometryRepresentation")

        # Annotate landmark location text to match glyph position
        create_text(f"{ab}. {name}", [0.2, 0.30], reg_name="text1",
                    fontsize=int(FONTSIZE * 1.5))

        # Text showing the file name for easy id of data
        create_text(s=os.path.splitext(os.path.basename(data_fid))[0],
                    position=[0.6, 0.30], reg_name="text2",
                    fontsize=int(FONTSIZE * 1.5))

        # Generate and rescale the colorbar/ colormap
        vsLUT, cbar = set_colormap_colorbar(vtk, position=[0.4, 0.325],
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
        Render()

        SaveScreenshot(os.path.join(save_path, f"{tag}.png"), renderView,
                       ImageResolution=VIEW_SIZE, TransparentBackground=1)

        # Clean up for next plot
        Hide(clip_vtk, renderView)
        delete_temp_objects(ruler=True, glyph=True, point=True, text=True)


def make_interface(fid, preset, save_path=os.getcwd()):
    """
    Project the volume  onto an arbitrarily defined 3D surface. In this case the
    surface defines the plate interface model from Charles Williams (2013).

    .. note::
        A neareset neighbor approach is used here (LinearKernel). The
        number of nearest neighbor points was explicitely chosen and
        seemed to provide the most realistic model values while simultaenously
        reducing the level of visual distortion/artefacts caused by
        interpolating a regular grid onto a smooth, dipping plane.
        The goal was to try to implement the least amount of
        smoothing/interpolation possible. Increasing the number of points will
        produce a smoother model but at the risk of damping out high/low
        amplitudes. Less number of points will provide a more accurate
        representation at the risk of looking patchy due to inconsistent grid
        spacing.
    """
    vtk = OpenDataFile(fid)
    RenameSource("surface", vtk)
    renderView = GetActiveView()
    renderView.InteractionMode = "3D"

    # Here we open a file I generated which is simply a VTK file that defines
    # an XYZ surface with the PointData being the Z values.
    interface = OpenDataFile("/Users/Chow/Documents/academic/vuw/data/"
                             "carto/interface/interface_utm60.vtk")

    interpolator = PointDatasetInterpolator(Input=vtk, Source=interface)

    # These specific values were chosen to best represent the surface, can be
    # changed to change the look of the final interpolated product
    interpolator.Kernel = "LinearKernel"
    interpolator.Kernel.KernelFootprint = "N Closest"
    interpolator.Kernel.NumberOfPoints = 8

    Show(interpolator, renderView, "UnstructuredGridRepresentation")

    # Hard set the camera as a top-down view over model, same as depth slices
    ResetCamera()
    renderView = GetActiveView()
    renderView.CameraPosition = [401399., 5567959., 1276057.]
    renderView.CameraFocalPoint = [401399., 5567959., -197500.0]
    renderView.CameraViewUp = [0, 1, 0]
    Render()

    # Standard accoutrements, colorbar, min/max values
    vsLUT, cbar = set_colormap_colorbar(interpolator, position=[0.249, 0.75],
                                        orientation="Vertical",)
    rescale_colorscale(vsLUT, src=interpolator, vtk=vtk, preset=preset)
    show_colorbar(interpolator)
    create_minmax_glyphs(interpolator, glyph_type="Sphere", scale_factor=20000)

    # Annotate the bottom-right corner to explain this is the interface
    text, _ = create_text(s="Interface", position=[0.625, 0.1],
                          reg_name="text1", fontsize=FONTSIZE * 2)

    # Annotate the file ID so we know what we're plotting
    create_text(s=os.path.splitext(os.path.basename(data_fid))[0],
                position=[0.25, 0.95], reg_name="text2", fontsize=FONTSIZE * 2)


    # Save screenshot, hide the surface and move on
    SaveScreenshot(os.path.join(save_path, f"interface.png"), renderView,
                   ImageResolution=VIEW_SIZE, TransparentBackground=1)
    Hide(interpolator, renderView)


def make_preplot(args):
    """
    Convenience function to plot extras such as coastline, srcs, rcvs
    """
    if args.coastline or args.extras:
        if args.verbose:
            print(f"\t\tPlotting coastline")
        plot_coastline()
    if args.events or args.extras:
        if args.verbose:
            print(f"\t\tPlotting event glyphs")
        plot_events()
    if args.stations or args.extras:
        if args.verbose:
            print(f"\t\tPlotting station glyphs")
        plot_stations()


if __name__ == "__main__":
    # Command line arguments to define behavior of script
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="+", help="file to plot with paraview")
    parser.add_argument("-o", "--output", type=str, help="output path",
                        default=os.getcwd())
    parser.add_argument("-p", "--preset", type=str, 
                        help="preset colormap and labels given file type")
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
    parser.add_argument("-i", "--interface", action="store_true",
                        help="plot projection of model onto plate interface "
                             "of Charles Williams",
                        default=False)
    parser.add_argument("-c", "--coastline", action="store_true",
                        help="plot the coastline", default=False)
    parser.add_argument("-e", "--events", action="store_true",
                        help="plot events as glyphs", default=False)
    parser.add_argument("-s", "--stations", action="store_true",
                        help="plot stations as glyphs", default=False)
    parser.add_argument("-E", "--extras", action="store_true",
                        help="shorthand to plot coastline, src glyphs and "
                             "rcv glyphs, same as '-ces'",
                        default=False,)
    parser.add_argument("-A", "--all", action="store_true", default=False,
                        help="shorthand to make all default slices, "
                             "same as '-xyzit'")
    parser.add_argument("-d", "--depth_cutoff_km", type=float, default=100,
                        help="For any vertical cross sections (Y-axis figure "
                             "normal to Z axis of volume), define the depth"
                             "cutoff of the screenshot as usually were not "
                             "interested in looking at the entire volume. "
                             "Units of km, positive values only.")
    parser.add_argument("-b", "--bounds", type=str,
                        help="Manually set the bounds of the colorbar, "
                             "overriding the default or preset bound values",
                        default=None)
    parser.add_argument("-v", "--verbose", action="store_true", default=False,
                        help="print output messages during the plotting")
    args = parser.parse_args()

    # Set up the active render view
    renderView = GetActiveViewOrCreate("RenderView")
    renderView.ViewSize = VIEW_SIZE
    renderView.InteractionMode = "2D"
    SetActiveView(renderView)

    for data_fid in args.files:
        if args.verbose:
            print(f"Plotting file {data_fid}")

        # ======================================================================
        # PATH CONFIGURATION
        # ======================================================================
        assert(os.path.exists(data_fid)), f"{data_fid} not found"

        save_fid = os.path.splitext(os.path.basename(data_fid))[0]
        save_path = os.path.join(args.output, save_fid)

        if not os.path.exists(save_path):
            os.makedirs(save_path)

        # ======================================================================
        # ASSIGN PRESET VALUES
        # ======================================================================
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
        if args.bounds:
            preset["bounds"] = [float(_) for _ in args.bounds.split(",")]

        # ======================================================================
        # DEPTH SLICES
        # ======================================================================
        zslices = []
        if args.default_zslices or args.all:
            zslices += ["surface", "2-20,2", "25-50,5"]
        if args.zslices:
            zslices += args.zslices

        if zslices:
            make_preplot(args)
            zslices = parse_slice_list(zslices)
            if args.verbose:
                print(f"\tGenerating Z slices for {zslices}")
            make_depth_slices(data_fid, zslices, preset, save_path=save_path)
            reset()

        # ======================================================================
        # X-NORMAL CROSS SECTIONS
        # ======================================================================
        if args.xslices or args.all:
            # Do not go for 0 or 100 % because you might end up off the model
            xslices = list(range(5, 100, 10))
            if args.verbose:
                print(f"\tGenerating X-normal slices for {xslices}")
            make_cross_sections(data_fid, xslices, normal="x", preset=preset,
                                depth_cutoff_km=args.depth_cutoff_km,
                                save_path=save_path)
            reset()

        # ======================================================================
        # Y-NORMAL CROSS SECTIONS
        # ======================================================================
        if args.yslices or args.all:
            # Possible to manually set the y-slice list here
            yslices = list(range(5, 100, 10))
            if args.verbose:
                print(f"\tGenerating Y-normal slices for {yslices}")
            make_cross_sections(data_fid, yslices, normal="y", preset=preset,
                                depth_cutoff_km=args.depth_cutoff_km,
                                save_path=save_path)
            reset()

        # ======================================================================
        # TRENCH NORMAL CROSS SECTIONS
        # ======================================================================
        if args.trench or args.all:
            if args.verbose:
                print("\tGenerating trench normal cross sections")
            make_trench_cross_sections(data_fid, preset,
                                       depth_cutoff_km=args.depth_cutoff_km,
                                       save_path=save_path)
            reset()

        # ======================================================================
        # PROJECTION ONTO INTERFACE
        # ======================================================================
        if args.interface or args.all:
            if args.verbose:
                print("\tGenerating interface projection")
            make_preplot(args)
            make_interface(data_fid, preset, save_path=save_path)
            reset()
           

