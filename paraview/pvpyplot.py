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
import sys
import math
import argparse
from subprocess import call
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
                             ColorBy, PointDatasetInterpolator,
                             GetOpacityTransferFunction, ImportPresets, Contour,
                             ProgrammableFilter)


Z = 2000.
# Strike parallel, 40 deg from X-axis
NORM = [-0.64, -0.76, 0.]
# NORM = [-0.45, -0.9, 0.]

# Pre-defined diagonals to make trench cross sections
DIAGONALS = {
    "Kaikoura":   ([226749., 5300535., Z], NORM),
    "Wellington": ([314007., 5426403., Z], NORM),
    "Flatpoint":  ([413092., 5433559., Z], NORM),
    "Castlepoint":([434724., 5471821., Z], NORM),
    "Akitio":     ([436980., 5511241., Z], NORM),
    # "Porangahau": ([467051., 5538717., Z], NORM),  # ORIGINAL
    "Porangahau": ([467051., 5535717., Z], NORM),
    "Intraplate": ([509140., 5515069., Z], [-0.23, -0.42, 0]),
    "Mahia_Psf":  ([577779., 5667301., Z], [-0.23, -0.42, 0]), 
    "Elsthorpe":  ([484394., 5581561., Z], NORM),
    "Napier":     ([489374., 5626518., Z], NORM),
    "Mohaka":     ([507922., 5670909., Z], NORM),
    "Mahia":      ([577779., 5667301., Z], NORM),
    "Gisborne":   ([588984., 5720001., Z], NORM),
    "Cvr":        ([417044., 5726459., Z], NORM),
    "Cook":       ([307699., 5384284., Z], NORM),
    "Psf_Cook":   ([307699., 5384284., Z], [-.45, -.45, 0]),
    "Psf_Cook_R": ([307699., 5384284., Z], [-0.45, 0.45, 0]),
    "Okataina":   ([463185., 5780787.0, Z], NORM),
    "Strait":     ([270065., 5402644., Z], [-0.4, -0.49, 0]),
    "Reyners":    ([314007., 5426403., Z], [-0.7, 0.7, 0]),
    "Henrys":     ([314007., 5426403., Z],[-0.9, 0.65, 0]),
    "Mahia_R":    ([577779., 5667301., Z], [-0.9, 0.65, 0]),
    "Cook_R":     ([307699., 5384284., Z], [-0.28, 0.27, 0]),
    # "Strait_R":   ([273419., 5377033., Z], [-0.49, 0.40, 0]),
    "Strait_R":   ([307699., 5384284., Z], [-0.25, 0.186, 0]),
    "Csfaults":   ([326337., 5409422., 0], [-0.324, 0.527, 0.0]),  # GOOD
    "Cspert_R":   ([307699., 5384284., Z], [-0.45, -0.45, 0.]),
    # "Csfaults":   ([316337., 5359422., 0], [-0.325, 0.325, 0.0]),  # DEEP
    "Csfaults_R": ([307699., 5384284., 0], [-0.527, -0.324, 0.0]),
    "Seamounts": ([475000, 5538861, Z], [-0.13, 0.11, 0]),
    "Seamounts_P":([466855, 5538861, Z], [-0.13, 0.11, 0]),
    "Seamounts_M":([577779, 5667301, Z], [-0.13, 0.11, 0]),
    "Tvz":        ([376534.0, 5650985.0, Z], [-0.2, 0.14, 0]),
    "Tvzpsf":    ([441230.0, 5776443.0, Z], [-0.125, 0.065, 0]),
    "Tvz_Rift1":   ([376534.3, 5650985.4, Z], [-0.13, 0.08, 0.0]),  # JGR Paper
    "Tvz_Rift2":   ([350000, 5650985.4, Z], [-0.13, 0.08, 0.0]),  # JGR Paper
    "Tvz_Rift3":   ([400000, 5650985.4, Z], [-0.13, 0.08, 0.0]),  # JGR Paper
    "Ruapehu_Psf": ([378927., 5651838., Z], NORM),
    "Rift":       ([463185., 5780787., Z], [-0.1298, 0.08665, 0.]),
    # "Rift":       ([246758., 5646174., Z], [-0.13, .216, 0.]),
    "Okataina_R": ([463185.0, 5780787.0, Z], [-0.63, 0.44, 0]),
    "Srcrcv":     ([255819.8, 5375745.8, Z], [-0.44, 0.11, 0]),
    "Centhawkbay": ([472312., 5683950.0, Z], [-.95, -.3, 0]),
}


def rgb_colors(c):
    """
    Paraview requires RGB color codes to set colors of items. This function
    converts python color characters to RGB to give to Paraview

    :type c: str
    :param c: python color code
    :return: list of float
    """
    return {"k": [0., 0., 0.], "w": [1., 1., 1.], "r": [1., 0., 0.],
            "b": [0., 0., 1.], "g": [0., 1., 0.], "y": [1., 1., 0.],
            "o": [1., .5, 0.], "c": [0., 1., 1.], "gray": [.5, .5, .5],
            "pink": [255., 0., 255.], "darkgray": [64., 64., 64.],
            "darkpurple":[102., 0., 104.],
            }[c]

class Preset(dict):
    """
    A more accesible dict object for storing preset values for defining
    colormaps, and the associated colorbars for different model types
    """
    def __init__(self, title="", cmap="Jet", invert=False, center=False, 
                 fmt="%.2f", rnd=10, bounds=False, nlabel=3, nvalues=28,
                 above_range=False, below_range=False, isosurfaces=None,
                 cdx=None, scale_units=1):
        """
        :title (str): Label for the colorbar
        :cmap (str): Colormap to define LUT
        :invert (bool): Invert the original bounds of the colormap
        :center (bool): Ensure that the middle value of the colormap is 0
        :fmt (str): String formatter for the range bounds on colorbar
        :rnd (int): Round range labels to a base value, if None, no rounding
        :bounds (bool or tuple or list): Can be overwritten using the '-b' arg
            * True: keep the colorbar bounds constant for every screenshot using 
                the min/max values of the entire volume.
            * False: use whatever default bounds are set by Paraview,
            * tuple/list: Manually set the bounds; must be in the form (min, max)
        :nlabel (int): Number of labels on the colorbar. Minimum is 2 to include
            min and max values
        :nvalues (int): Number of segmentations in the colorbar
        :above_range (list or None): Use above range color
        :below_range (list or None): Use below range color
        :isosurfaces (list of float): If contour lines are turned on, define the
            values which are contoured
        :cdx (int): If no isosurfaces are provided, cdx (contour dx) determines
            the spacing between isosurfaces. e.g. for a Vs model, a useful value
            would be 500, meaning isolines are separated by 500 m/s
        """
        self.title = title
        self.cmap = cmap
        self.invert = invert
        self.center = center
        self.fmt = fmt
        self.rnd = rnd
        self.bounds = bounds
        self.nlabel = nlabel
        self.nvalues = nvalues
        self.above_range = above_range
        self.below_range = below_range
        self.isosurfaces = isosurfaces
        self.cdx = cdx
        self.scale_units = scale_units

    def __setattr__(self, key, value):
        self[key] = value

    def __getattr__(self, key):
        return self[key]


PRESETS = {
    "tvz": Preset(
        title="Vs [m/s]", cmap="erdc_iceFire_H", invert=True, center=False,
        fmt="%.1f", rnd=10, bounds=False, nlabel=3, nvalues=33, cdx=500,
        above_range=rgb_colors("k"), below_range=rgb_colors("k")
    ),
    "henrys_ratio_vpvs": Preset(
        title="Vp/Vs Ratio", cmap="Black-Body Radiation",
        invert=False, center=False, fmt="%.2f", rnd=None,
        bounds=[1.5, 1.9], nlabel=4, nvalues=28,
    ),
    "perturbation": Preset(
        title="Perturbation [m/s]", cmap="Cool to Warm (Extended)",
        invert=True, center=True, fmt="%.1f", rnd=None, bounds=True, nlabel=3,
        nvalues=21, cdx=50
    ),
    "psf_custom": Preset(
        title="PSF [1E-7 $s^{3}m^{-1}$]", cmap="Black-Body Radiation", invert=True,
        center=False, fmt="%.2f", rnd=None, bounds=[0, 2.1], nlabel=2,
        nvalues=11, scale_units=1E7,
    ),
    "psfn_custom": Preset(
        title="PSF [1E-7 $s^{3}m^{-1}$]", cmap="Black-Body Radiation", invert=False,
        center=False, fmt="%.2f", rnd=None, bounds=[-6, 0], nlabel=2,
        nvalues=11, scale_units=1E7,
    ),
    # Zeroth Moment
    "zm": Preset(
        title="PSF [1E-6 $s^{3}m^{-1}$]", cmap="Black, Blue and White", invert=True,
        center=False, fmt="%.1f", rnd=None, bounds=[0, 5], nlabel=2,
        nvalues=21, isosurfaces=[.8],
        scale_units=1E7,
    ),
    # Incorrect finite-difference sign means some PSFVs have to be flipped
    "psfv_neg": Preset(
        title="PSFV [1E-6 m^3 s^2]", cmap="Black, Blue and White", invert=True,
        center=False, fmt="%.1f", rnd=None, bounds=[0, 10], nlabel=2,
        nvalues=21, # isosurfaces=[1., 2., 3.],
        scale_units=-1E6,
    ),
    "psf": Preset(
        title="PSF [1E-7 s^3 m^-1]", cmap="Cool to Warm (Extended)",
        invert=True, center=True, fmt="%.2f", rnd=None, # bounds=False,
        nlabel=3, nvalues=21, bounds=[-2, 2], scale_units=1E7,
    ),
    # Incorrect finite-difference sign means some PSFVs have to be flipped
    "psf_neg": Preset(
        title="PSF [1E-7 m^3 s^-2]", cmap="Cool to Warm (Extended)",
        invert=True, center=True, fmt="%.2f", rnd=None,  # bounds=False,
        nlabel=3, nvalues=33, bounds=[-4, 4], scale_units=-1E8
    ),
    "donna_vpvs": Preset(
        title="Vp/Vs Ratio", cmap="Cool to Warm (Extended)", invert=False,
        center=False, fmt="%.2f", rnd=None, bounds=[1.55, 1.9], nlabel=4,
        nvalues=21, isosurfaces=[1.5, 1.75, 2., 2.25]
    ),
    "trench_vs": Preset(
        title="Vs [m/s]", cmap="Rainbow Desaturated", invert=True, center=False,
        fmt="%.1f", rnd=10, bounds=[1750, 5500], nlabel=4, nvalues=33,
        cdx=500, above_range=False #rgb_colors("gray"),
    ),
    "trench_vp": Preset(
        title="Vs [m/s]", cmap="Rainbow Desaturated", invert=True, center=False,
        fmt="%.1f", rnd=10, bounds=[3500, 9250], nlabel=4, nvalues=33,
        cdx=500, above_range=False # rgb_colors("gray"),
    ),
    "model_rho": Preset(
        title="Density [kg m^-3]", cmap="Rainbow Desaturated", invert=True,
        center=False, fmt="%.1f", rnd=10, bounds=False, nlabel=3,
        nvalues=33
    ),
    "model_vp": Preset(
        title="Vp [m/s]", cmap="Rainbow Desaturated", invert=True, center=False,
        fmt="%.1f", rnd=10, bounds=False, nlabel=3, nvalues=33, cdx=500,
    ),
    "model_vp_kms": Preset(
        title="Vp [km/s]", cmap="Rainbow Desaturated", invert=True, center=False,
        fmt="%.2f", rnd=None, bounds=False, nlabel=3, nvalues=33, cdx=500,
        scale_units=1E-3
    ),
    "model_vs": Preset(
        title="Vs [m/s]", cmap="Rainbow Desaturated", invert=True, center=False,
        fmt="%.1f", rnd=10, bounds=False, nlabel=3, nvalues=33, cdx=500,
    ),
    "model_vs_kms": Preset(
        title="Vs [km/s]", cmap="Rainbow Desaturated", invert=True, 
        center=False, fmt="%.2f", rnd=None, bounds=False, nlabel=3, nvalues=33, 
        cdx=500, scale_units=1E-3
    ),
    "gradient_vp_kernel": Preset(
        title="Vp Gradient [m^-2 s^2]", cmap="Cool to Warm (Extended)",
        invert=True, center=True, fmt="%.1E", rnd=None, bounds=True,
        nlabel=3, nvalues=33
    ),
    "gradient_vs_kernel": Preset(
        title="Vs Grad. [1E-7 m^-2 s^2]", cmap="Cool to Warm (Extended)",
        invert=True, center=True, fmt="%.0f", rnd=None, bounds=[-5,5],
        nlabel=3, nvalues=33, scale_units=1E7, 
    ),
    "update_vp": Preset(
        title="Vp Update [ln(m28/m00)]", cmap="Blue Orange (divergent)",
        invert=True, center=True, fmt="%.02f", rnd=None, bounds=[-.1,.1],
        nlabel=3, nvalues=31,
    ),
    "update_vs": Preset(
        title="Vs Update [ln(m28/m00)]", cmap="Blue Orange (divergent)",
        invert=True, center=True, fmt="%.02f", rnd=None, bounds=[-.15,.15],
        nlabel=3, nvalues=31,
    ),
    "update_vpvs": Preset(
        title="Vp/Vs Update [ln(m28/m00)]", cmap="Blue Orange (divergent)",
        invert=True, center=True, fmt="%.02f", rnd=None, scale_units=1,
        bounds=True, nlabel=3, nvalues=33,
    ),
    "update_poissons": Preset(
        title="Poisson's Update [ln(m28/m00)]", cmap="Blue Orange (divergent)",
        invert=True, center=True, fmt="%.02f", rnd=None,
        bounds=[-.5, .5], nlabel=3, nvalues=33,
    ),
    "update_shear": Preset(
        title="Shear Modulus Update [ln(m28/m00)]",
        cmap="Blue Orange (divergent)",
        invert=True, center=True, fmt="%.02f", rnd=None,
        bounds=[-.5, .5], nlabel=3, nvalues=33,
    ),
    "ratio_poissons": Preset(
        title="Poisson's Ratio", cmap="Blue - Green - Orange", invert=False,
        center=False, fmt="%.1f", rnd=None, bounds=[0.1, 0.4], nlabel=3,
        nvalues=28, isosurfaces=[0.1, 0.2, 0.3, 0.4]
    ),
    "ratio_vpvs": Preset(
        title="Vp/Vs Ratio",
        # cmap="erdc_iceFire_H",
        cmap="Yellow - Gray - Blue",
        invert=False,
        center=False, fmt="%.2f", rnd=None,
        isosurfaces=[1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.1, 2.2]
    ),
    "modulus_shear": Preset(
        title="Shear Modulus [GPa]", cmap="Jet", #cmap="Inferno (matplotlib)",
        invert=True, center=False, fmt="%.0f", rnd=None,
        bounds=[5, 65], nlabel=3, nvalues=12, above_range=[.78, .78, .78],
        cdx=10,
    ),
}


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


def normal_to_angle(x, y):
    """
    Take two normal vectors and return the angle that they give.

    :type x: float
    :param x: x normal
    :type y: float
    :param y: y normal
    :rtype: float
    :return: angle created by the two normals
    """
    return math.atan2(y, x) * 180 / math.pi


def convert_coords(data, x0=171312, y0=5286950, scale=1e-3):
    """
    Convert landmark locatiosn from UTM60S to a domain with origin at mesh LLC
    """
    return [(data[0] - x0) * scale, (data[1] - y0) * scale, data[2] * scale]


def remove_whitespace(fid):
    """
    Use Imagemagick to remove any whitespace from a .png file
    """
    call([f"convert {fid} -fuzz 1% -trim +repage {fid}"], shell=True)


def project_point_trench_normal(point, origin=None, normal=None):
    """
    Project a point onto the trench normal which is defined by:
    normal = [-.7, .7., 0]
    origin = [455763., 5547040., 0.]
    or in math terms: -.7x + .7y + 0z = 3563893

    TO DO:
    make this generalized for a given normal, origin

    .. note::
        Assuming the plane normal is parallel to the Z-axis so that the Z point
        remains unchanged.
    """
    # if normal is None:
    #     normal == [-.7, .7, 0.]
    # if origin is None:
    #     origin = [455763, 5547040., 0.]
    #
    # NX, NY, NZ = normal
    # OX, OY, OZ = origin

    D = 3563893
    x, y, z = point
    t = (D + .7 * (x - y)) / (2 * .7 ** 2)
    x0 = x - .7 * t
    y0 = y + .7 * t

    return [x0, y0, z]


def parse_bounds(bounds):
    """
    Parse the 'bounds' argument, which is allowed to either be:
    1) two float values separated by a comma, negatives okay but there must not
        be a space after the command line input otherwise argparse will complain
        e.g. pvpython .... -b-5000,5000 ....
    2) a boolean where True means respect global bounds for any slice made, and
        False means use the local bounds of the slice to set colorbar
    3) None, meaning let Paraview decide the bounds automatically

    :type bounds: str
    :param bounds: string to parse
    """
    if bounds.lower() in "none":
        return None
    elif bounds.lower() in ["t", "true"]:
        return True
    elif bounds.lower() in ["f", "false"]:
        return False
    else:
        return [float(_) for _ in bounds.split(",")]


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


def scale_input(data, scale_units=1, scale_coords=1E-3):
    """
    Use the ProgrammableFilter to scale the input units, e.g. to get from m/s
    to km/s, and the coordinates from e.g. m to km. Also allows setting the 
    origin to 0 rather than an arbitrary value in UTM 60S.
    Should be run immediately after opening the input file.

    :type input
    :type scale_by: float
    :param scale_by: value to scale input units by.
    :return:
    """
    renderView = GetActiveViewOrCreate('RenderView')

    Show(data)
    progFilt = ProgrammableFilter(registrationName='ProgrammableFilter',
                                  Input=data)
    progFilt.Script = f"""from paraview import vtk
pdi = self.GetInput()
pdo = self.GetOutput()

numPoints = pdi.GetNumberOfPoints()
newPoints=vtk.vtkPoints()
min_x, _, min_y, _, _, _ = pdi.GetBounds()
for i in range(0, numPoints):
    coords = pdi.GetPoint(i)
    x, y, z = coords[:3]
    x -= min_x
    y -= min_y

    # Convert to units of km
    x *= {scale_coords}
    y *= {scale_coords}
    z *= {scale_coords}
    newPoints.InsertPoint(i, x, y, z)

# Set the new coordinate system
pdo.SetPoints(newPoints)

# Create a new array that will hold the converted values
ivals = pdi.GetPointData().GetScalars()

ca = vtk.vtkFloatArray()
ca.SetName(ivals.GetName())
ca.SetNumberOfComponents(1)
ca.SetNumberOfTuples(numPoints)

# Copy the values over element by element and convert
for i in range(0, numPoints):
  ca.SetValue(i, ivals.GetValue(i) * {scale_units})

# Set the new values
pdo.GetPointData().AddArray(ca)"""

    quantity = progFilt.PointData.GetArray(0).Name  # e.g. model_init_vp
    pfDisplay = Show(progFilt, renderView, "UnstructuredGridRepresentation")
    ColorBy(pfDisplay, quantity)
    Hide(data)
    Hide(progFilt)

    return progFilt


def clip_domain(vtk, xmin=None, xmax=None, ymin=None, ymax=None, zmin=None,
                zmax=None, color_by=True):
    """
    Clip away unwanted parts of the domain to effectively zoom in on a section
    of the volume. e.g. clipping the bottom to 30km by setting zmin=30

    :type vtk:
    :param vtk: opened data file
    :type depth: float
    :param depth: depth to clip in the units of the volume
    :param xmin: clip everything bigger than xmin,  +x direction
    :param xmax: clip everything smaller than xmax,  -x direction
    :param ymin: clip everything bigger than ymin, +y direction
    :param ymax: clip everything smaller than ymax, -y direction
    :param zmin: clip everything above zmin, +z direction (+Z UP)
    :param zmax: clip everything below zmax, -z direction (+Z UP)
    """
    renderView = GetActiveView()

    def clip_fx(vtk, origin, normal):
        clip_ = Clip(vtk)
        clip_.ClipType = "Plane"
        clip_.ClipType.Origin = origin
        clip_.ClipType.Normal = normal
        return clip_

    clip = vtk
    if xmin:
        clip = clip_fx(clip, origin=[xmin, 0., 0.], normal=[-1, 0., 0.])
    if xmax:
        clip = clip_fx(clip, origin=[xmax, 0., 0.], normal=[1, 0., 0.])
    if ymin:
        clip = clip_fx(clip, origin=[0., ymin, 0.], normal=[0., -1., 0.])
    if ymax:
        clip = clip_fx(clip, origin=[0., ymax, 0.], normal=[0., 1., 0.])
    if zmin:
        clip = clip_fx(clip, origin=[0., 0., zmin], normal=[0., 0., -1])
    if zmax:
        clip = clip_fx(clip, origin=[0., 0., zmax], normal=[0., 0., 1])

    # Make sure to se the color of the new clipped domain for lookup table
    try:
        quantity = clip.PointData.GetArray(0).Name  # e.g. model_init_vp
        clipDisplay = Show(clip, renderView, "UnstructuredGridRepresentation")
        # Tuple required otherwise in ColorBy() otherwise Paraview throws a
        # > RuntimeError: invalid association string 'NONE'
        if color_by:
            ColorBy(clipDisplay, ("POINTS", quantity))
        else:
            ColorBy(clipDisplay, None)
    except AttributeError:
        pass

    Hide(vtk)

    return clip


def depth_slice(vtk, depth):
    """
    Slice a horizontal plane through a volume normal to the Z axis

    :type vtk:
    :param vtk: opened data file
    :type depth: float
    :param depth: depth at which to slice through the model
    """
    args = parse_args()

    # Generate the Z-normal slice and set at a specific depth
    slice_vtk = Slice(vtk)
    slice_vtk.SliceType = "Plane"
    slice_vtk.SliceType.Normal = [0., 0., 1.]
    origin = slice_vtk.SliceType.Origin
    origin = [origin[0], origin[1], abs(depth) * -1E3 * args.coord_scale]
    slice_vtk.SliceType.Origin = origin

    # Apply the slice and render the new view
    renderView = GetActiveView()
    Show(slice_vtk, renderView)
    RenameSource(f"z_{depth}km", slice_vtk)
    Hide3DWidgets(proxy=slice_vtk)

    return slice_vtk


def clip_bottom(vtk, depth):
    """
    Clip away the bottom of the volume, e.g. if the mesh bottom is at 400km
    depth but you only want to plot top 30km, remove the remaining 370km

    :type vtk:
    :param vtk: opened data file
    :type depth: float
    :param depth: depth to clip in the units of the volume
    :return:
    """
    clip = Clip(vtk)
    clip.ClipType = "Plane"
    clip.ClipType.Origin = [0., 0., depth]
    clip.ClipType.Normal = [0., 0., -1.]

    return clip


def contour_lines(vtk, preset, scale=1., color=None, reg_name="contour",
                  line_width=1.):
    """
    Set countour lines for the given data file

    :type vtk:
    :param vtk: opened data file
    :type preset: dict
    :param preset: the preset choices for how to deal with the colorscale/ map
    :type dx: int
    :param dx: separation for contour lines
    :type color: list of float
    :param color: RGB color for contour lines
    :type reg_name: str
    :param reg_name: name to assign to the source object
    """
    # Contour lines require one of these presets to function
    if (preset.isosurfaces is None) and (preset.cdx is None):
        return

    renderView = GetActiveView()
    contour = Contour(registrationName=reg_name, Input=vtk)

    # Use server manager to determine min and max values of the source
    if preset.isosurfaces is None:
        dataPlane = servermanager.Fetch(vtk)
        dataPlane = dataset_adapter.WrapDataObject(dataPlane)
        data_points = list(dataPlane.PointData[0])
        min_point = myround(min(data_points), preset.cdx)
        while min_point > min(data_points):
            min_point -= preset.cdx
        max_point = myround(max(data_points), preset.cdx)
        while max_point < max(data_points):
            max_point += preset.cdx
        isosurfaces = list(range(min_point, max_point, preset.cdx))
        contour.Isosurfaces = isosurfaces
    else:
        contour.Isosurfaces = preset.isosurfaces

    # Turn off colorbar and make the contour lines a solid color
    contourDisplay = Show(contour, renderView, "GeometryRepresentation")
    contourDisplay.LineWidth = line_width
    # Work-around to avoid RuntimeError thrown by ColorBy() with arg 'None'
    contourDisplay.ColorArrayName = ["POINTS", ""]
    contourDisplay.Scale = [1., 1., scale]
    ColorBy(contourDisplay, None)

    if color is None:
        color = rgb_colors("w")
    contourDisplay.AmbientColor = color
    contourDisplay.DiffuseColor = color


def contour_overlay(fid, origin=None, normal=None, depth=None, color="pink",
                    scale=1, xmin=None, ymin=None, ymax=None, xmax=None,
                    zmin=None):
    """

    :param fid:
    :param origin:
    :param normal:
    :param scale:
    :return:
    """
    vtk = OpenDataFile(fid)
    vtk = scale_input(vtk, scale_coords=args.coord_scale)
    if xmin:
        vtk = clip_domain(vtk, xmin=xmin, ymin=ymin, ymax=ymax, xmax=xmax,
                          zmin=zmin)
    Hide(vtk)

    slice = Slice(vtk)
    slice.SliceType = "Plane"
    if normal:
        slice.SliceType.Normal = normal
    else:
        slice.SliceType.Normal = [0., 0., 1.] # normal

    if origin:
        slice.SliceType.Origin = origin
    elif depth:
        origin = slice.SliceType.Origin
        origin = [origin[0], origin[1], abs(depth) * -1E3 * args.coord_scale]
        slice.SliceType.Origin = origin
    Show(slice, GetActiveView())

    contour = Contour(Input=slice)
    contour.Isosurfaces = [2E-7]

    contourDisplay = Show(contour, GetActiveView(), "GeometryRepresentation")
    contourDisplay.LineWidth = 6.
    contourDisplay.Scale = [1., 1., scale]

    contourDisplay.ColorArrayName = ["POINTS", ""]
    ColorBy(contourDisplay, None)
    contourDisplay.AmbientColor = rgb_colors(color)
    contourDisplay.DiffuseColor = rgb_colors(color)

    Hide(slice, GetActiveView())

    # import pdb;pdb.set_trace()
def cross_section(vtk, normal, origin, name):
    """
    Cut a vertical cross section through a volume, with the plane parallel to
    the Z axis.

    :type vtk: paraview.servermanager.LegacyVTKReader
    :param vtk: opened data file
    :type normal: list of float
    :param normal: normal vector of the plane
    :type origin: list of float
    :param origin: origin points for the xsection
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


def create_ruler(point1, point2, label="", graduation=100E3, ticknum=None,
                 axis_color=None, font_color=None, justification="Left",
                 reg_name="ruler", line_width=3., fontsize=None):
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
    if ticknum:
        rulerDisplay.NumberOfTicks = ticknum
    else:
        rulerDisplay.RulerMode = 1
        rulerDisplay.Graduation = graduation
    rulerDisplay.AxisColor = axis_color or LINECOLOR
    rulerDisplay.Color = font_color or FONTCOLOR
    rulerDisplay.FontSize = fontsize or FONTSIZE
    rulerDisplay.AxisLineWidth = line_width
    rulerDisplay.Justification = justification

    Hide3DWidgets(proxy=ruler)

    return reg_name

def create_text(s, position, fontsize=None, color=None, bold=False, 
                justification="center", reg_name="text"):
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
    :type justification: str
    :param justification: text justification, 'left', 'right', 'center'
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
    textDisplay.Color = color or LINECOLOR
    textDisplay.Justification = justification.title()

    return text, reg_name


def create_cone_glyph(origin, reg_name_point="point", reg_name_glyph="glyph"):
    """
    Create a 3D cone glyph used to denote points at the surface
    :param origin:
    :param reg_name:
    :return:
    """
    point = PointSource(registrationName=reg_name_point)
    point.Center = origin
    glyph = Glyph(Input=point, GlyphType="Cone",
                  registrationName=reg_name_glyph)
    glyph.GlyphType.Direction = [0., 0., -1.]
    glyph.ScaleFactor = 10000
    Show(glyph, renderView, "GeometryRepresentation")


def mark_surface_point(origin, normal, color=None, z_value=None, size=10.,
                       reg_name_point="point", reg_name_glyph="glyph"):
    """
    Create a 3D cone glyph used to denote points at the surface

    :type origin: list of float
    :param origin: center point for mark
    :type normal: list of float
    :param normal: normal of the given slice, to rotate the glyph to orientation
    :type color: list of float
    :param color: rgb color for the mark, default black
    :type z_value: float
    :param z_value: height of the mark, defaults 9000. to push the glyph above
        any topography, might need to play with this value.
    """
    renderView = GetActiveView()

    if color is None:
        color = rgb_colors("k")
    if z_value is None:
        z_value = 9000.

    # Generate point at origin
    point = PointSource(registrationName=reg_name_point)
    point.Center = [origin[0], origin[1], z_value]

    # Create a 2D inverted triangle glyph, rotated to the correct orientation
    glyph = Glyph(Input=point, GlyphType="2D Glyph",
                  registrationName=reg_name_glyph)
    glyph.GlyphType.GlyphType = "Triangle"
    rotate = normal_to_angle(normal[0], normal[1])
    # glyph.GlyphTransform.Rotate = [-90, 90., rotate]
    glyph.GlyphTransform.Rotate = [90, 90., rotate]
    glyph.ScaleFactor = size

    # Set the look of the displayed glyphs
    glyphDisplay = Show(glyph, renderView, "GeometryRepresentation")
    glyphDisplay.AmbientColor = color
    glyphDisplay.DiffuseColor = color
    glyphDisplay.SetRepresentationType("Surface With Edges")
    glyphDisplay.LineWidth = 2.


def create_ruler_grid_axes_depth_slice(src, tick_spacing_km=50, top=True,
                                       bottom=True, left=True, right=True,
                                       fontsize=None, lw=6):
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
    if fontsize is None:
        fontsize = FONTSIZE
    tick_spacing_m = tick_spacing_km * 1E3
    x, y, z = get_coordinates(src)

    # Figure out an acceptable end point that allows us to have even grid
    # spacings but is as close to the actual volume end point as possible
    ruler_x_bot = [max(x), min(y), max(z)]
    ruler_x_top = [max(x), max(y), max(z)]

    ruler_y_lft = [min(x), max(y), max(z)]
    ruler_y_rgt = [max(x), max(y), max(z)]

    if bottom:
        create_ruler(point1=[min(x), min(y), max(z)], point2=ruler_x_bot,
                     graduation=tick_spacing_km, line_width=lw,
                     # label=f"[X] (dx={tick_spacing_km}km)",
                     reg_name="ruler_bottom", fontsize=fontsize)

    if top:
        create_ruler(point1=[min(x), max(y), max(z)], point2=ruler_x_top,
                     graduation=tick_spacing_km, line_width=lw,
                     reg_name="ruler_top", fontsize=fontsize)

    if left:
        create_ruler(point1=[min(x), min(y), max(z)], point2=ruler_y_lft,
                     graduation=tick_spacing_km, line_width=lw,
                     # label=f"[Y]\n(dy={tick_spacing_km}km)",
                     reg_name="ruler_left", fontsize=fontsize)

    if right:
        create_ruler(point1=[max(x), min(y), max(z)], point2=ruler_y_rgt,
                     graduation=tick_spacing_km * 1E3, line_width=lw,
                     reg_name="ruler_right", fontsize=fontsize)


# def create_ruler_grid_axes_cross_section(src, normal, tick_spacing_m=50E3,
#                                          depth_cutoff_km=100., scale=5.,
#                                          dz_km=10., xlabel=None,
#                                          flip_view=False):
#     """
#     Generate axis grids for cross section plots using rulers, used to define
#     heights and distances
#
#     :type mode: str
#     :param mode: depending on what plot is being made, different rulers need to
#         be defined. The available modes are:
#         x: normal to the X-axis, parallel to the Y-axis, slices are made for
#             fixed values of X,
#         y:normal to the Y-axis, parallel to the X-axis, slices are made for
#             fixed values of Y,
#         trench_normal: normal to the Hikurangi trench (trench defined as
#             40deg from X-axis)
#         trench_parallel: parallel to the Hikurangi trench
#     """
#     # In order to set the camera, annotations, etc. in a general fashion,
#     # we need to determine the corner grid locations of the slice
#     x, y, z = get_coordinates(src)
#     ruler_origin = [min(x), max(y), min(z) * scale]  # bottom left
#
#     origin
#
#     # Vertical scaling is constant for any cross sections
#     dist_m_v = depth_cutoff_km * 1E3 * args.coord_scale
#     ruler_v = [min(x), max(y), min(z) + dist_m_v]
#
#     # Determine the orientation of the horizontal axis
#     if normal in ["x", "y"]:
#         if normal == "x":
#             axis = y
#         elif normal == "y":
#             axis = x
#
#         # Create a ruler for axis grid scale
#         dist_m_h = max(axis) - min(axis)
#
#         if normal == "x":
#             ruler_h = [min(x), max(y) - dist_m_h, min(z) * scale]
#         elif normal == "y":
#             ruler_h = [min(x) + dist_m_h, max(y), min(z) * scale]
#
#     else:
#         if flip_view:
#             # Rulers for trench parallel differ because our viewing angle
#             # changes which swaps the definition of the positive direction on
#             # the horizontal axis
#             ruler_v = [min(x), min(y), min(z) + dist_m_v]
#             ruler_origin = [min(x), min(y), min(z) * scale]
#
#         # Use trig to figure out the length of the angled cut
#         slice_length_m = ((max(y) - min(y)) ** 2 + (max(x) - min(x)) ** 2) ** .5
#
#
#         # Calculate the end points of the hypotenuse based on the side lengths
#         angle_rad = math.atan(normal[0] / normal[1])
#         ruler_h = [ruler_origin[0] + slice_length_m * math.cos(angle_rad),
#                    ruler_origin[1] - slice_length_m * math.sin(angle_rad),
#                    ruler_origin[2]]
#
#     # Create the 'X' axis ruler along the bottom
#     create_ruler(point1=ruler_origin, point2=ruler_h,
#                  graduation=tick_spacing_m,
#                  label=f"{xlabel}={slice_length_m/1E3:.0f}km "
#                        f"(dh={int(tick_spacing_m * 1E-3)}km)",
#                  reg_name="ruler1")
#     # For the vertical axis, make the ruler custom so it matches up with the
#     # data axis grid
#     create_ruler(point1=ruler_v, point2=ruler_origin,
#                  ticknum=int(depth_cutoff_km / dz_km) + 1,
#                  label=f"Z={depth_cutoff_km}km\n"
#                        f"{int(scale)}x scale\n"
#                        f"(dz={int(dz_km)}km)",
#                  reg_name="ruler2", )

def ruler_axes_cross_section(src, normal, dz=50, dh=50, scale=5.,
                             flip_view=False):
    """
    Generate axis grids for cross section plots using rulers, used to define
    heights and distances

    :type mode: str
    :param mode: depending on what plot is being made, different rulers need to
        be defined. The available modes are:
        x: normal to the X-axis, parallel to the Y-axis, slices are made for
            fixed values of X,
        y:normal to the Y-axis, parallel to the X-axis, slices are made for
            fixed values of Y,
        trench_normal: normal to the Hikurangi trench (trench defined as
            40deg from X-axis)
        trench_parallel: parallel to the Hikurangi trench
    """
    # In order to set the camera, annotations, etc. in a general fashion,
    # we need to determine the corner grid locations of the slice
    x, y, z = get_coordinates(src)
    ruler_origin = [min(x), max(y), min(z) * scale]  # bottom left

    # Always stop the Z-axis ruler at 0 depth, easiest and intuitive
    ruler_v = [min(x), max(y), 0]

    if flip_view:
        # Rulers for trench parallel differ because our viewing angle
        # changes which swaps the definition of the positive direction on
        # the horizontal axis
        ruler_v = [min(x), min(y), 0]
        ruler_origin = [min(x), min(y), min(z) * scale]

    # Use trig to figure out the length of the angled cut
    slice_length_m = ((max(y) - min(y)) ** 2 + (max(x) - min(x)) ** 2) ** .5

    # Calculate the end points of the hypotenuse based on the side lengths
    angle_rad = math.atan(normal[0] / normal[1])
    ruler_h = [ruler_origin[0] + slice_length_m * math.cos(angle_rad),
               ruler_origin[1] - slice_length_m * math.sin(angle_rad),
               ruler_origin[2]]

    create_ruler(point1=ruler_origin, point2=ruler_h, reg_name="ruler1",
                 graduation=dh, label="")
    create_ruler(point1=ruler_v, point2=ruler_origin, label="",
                 reg_name="ruler2", graduation=dz * scale)


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


def set_depth_slice_data_axis_grid(src, xtitle="", ytitle="", round_to=50,
                                   spacing_km=100):
    """
    Create the data axis grid that labels the tick marks and axis
    """
    renderView = GetActiveView()
    display = GetDisplayProperties(src, view=renderView)
    display.DataAxesGrid.GridAxesVisibility = 1
    display.DataAxesGrid.XTitle = xtitle
    display.DataAxesGrid.YTitle = ytitle
    display.DataAxesGrid.ShowGrid = 0
    display.DataAxesGrid.ShowTicks = 1
    display.DataAxesGrid.ShowEdges = 1
    display.DataAxesGrid.XTitleFontSize = FONTSIZE
    display.DataAxesGrid.YTitleFontSize = FONTSIZE
    display.DataAxesGrid.XLabelFontSize = FONTSIZE
    display.DataAxesGrid.YLabelFontSize = FONTSIZE
    
    # Figure out how to label the tickmarks based on the edges of domain
    x, y, z = get_coordinates(src)
    # Check if we need to round the scaling values
    min_x, min_y, max_x, max_y = \
        int(min(x)), int(min(y)), int(max(x)), int(max(y))
    if min_x % round_to:
        min_x = int(myround(min_x, round_to, "up"))
    if max_x % round_to:
        max_x = int(myround(max_x, round_to, "down"))
    if min_y % round_to:
        min_y = int(myround(min_y, round_to, "up"))
    if max_y % round_to:
        max_y = int(myround(max_y, round_to, "down"))

    while max_x < max(x):
        max_x += spacing_km
    while max_y < max(y):
        max_y += spacing_km

    xlabels = list(range(min_x, max_x, spacing_km))
    # ylabels = list(range(int(min(y)) + spacing_km, int(max(y)), spacing_km))
    ylabels = list(range(min_y, max_y, spacing_km))

    display.DataAxesGrid.XAxisLabels = xlabels 
    display.DataAxesGrid.YAxisLabels = ylabels 

    # display.DataAxesGrid.XLabelFontFamily = "Arial"
    # display.DataAxesGrid.XLabelBold = 1



def scale_vertical_grid(source, scale=1, zgrid_km=None, camera_parallel=False):
    """
    For cross-sections. Scale the vertical grid (e.g. 2x exagerration), and set
    grid lines showing corresponding depth values.

    :type source: paraview.servermanager.Slice
    :param source: the given slice that we want to make an axis grid for
    :type z_grid_km: list of int
    :param z_grid_km: [min_z_km, max_z_km, dz_km]
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
    display.DataAxesGrid.ShowEdges = 0

    display.DataAxesGrid.XAxisUseCustomLabels = 1
    display.DataAxesGrid.XAxisLabels = []

    display.DataAxesGrid.YAxisUseCustomLabels = 1
    display.DataAxesGrid.YAxisLabels = []

    display.DataAxesGrid.ZAxisUseCustomLabels = 1

    display.Scale = [1., 1., scale]
    display.DataAxesGrid.Scale = [1., 1., scale]

    # Ensuring that all the depth values are formatted properly before ranging
    # Create grid lines based on min, max depth values
    if zgrid_km:
        min_z, max_z, dz = zgrid_km
        zrange = range(min_z, max_z + dz, dz)

        zrange = [_ * scale for _ in zrange]   # need to scale range too
        display.DataAxesGrid.ZAxisLabels = list(zrange)
        display.DataAxesGrid.ZLabelFontSize = 1


def set_colormap_colorbar(vtk, position, orientation, colormap=None,
                          title=None, fontsize=None, color=None, thickness=None,
                          length=None, text_position="top",
                          justification="Left"):
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
    # # Import any external colormaps that may be defined by XML files
    # for fid in glob("/Users/Chow/Documents/subduction/packages/simutils/"
    #                    "paraview/lut_presets/*.xml"):
    #     ImportPresets(fid)

    # Change the colormap to the desired preset value
    quantity = vtk.PointData.GetArray(0).Name  # e.g. model_init_vp
    vsLUT = GetColorTransferFunction(quantity)
    vsLUT.ApplyPreset(colormap or preset.cmap, True)
    vsLUT.UseAboveRangeColor = 0
    vsLUT.UseBelowRangeColor = 0
    vsLUT.NumberOfTableValues = preset.nvalues
    if preset.invert:
        vsLUT.InvertTransferFunction()

    # Create the colorbar and set a common look
    cbar = GetScalarBar(vsLUT)
    cbar.Title = title or preset.title
    cbar.ComponentTitle = ""
    cbar.AutoOrient = 0
    cbar.Orientation = orientation
    cbar.Position = position
    cbar.RangeLabelFormat = preset.fmt
    cbar.AutomaticLabelFormat = 0
    cbar.LabelFormat = preset.fmt

    cbar.AddRangeLabels = 1
    cbar.ScalarBarThickness = thickness or 35
    cbar.ScalarBarLength = length or 0.15
    cbar.TitleFontSize = fontsize or FONTSIZE
    cbar.LabelFontSize = fontsize or FONTSIZE
    cbar.TitleColor = color or FONTCOLOR
    cbar.LabelColor = color or FONTCOLOR
    cbar.TitleJustification = justification

    # !!! This doesnt work as expected?
    if text_position in ["top", "right"]:
        cbar.TextPosition = "Ticks left/bottom, annotations right/top"
    elif text_position in ["bottom", "left"]:
        cbar.TextPosition = "Ticks right/top, annotations left/bottom"
    else:
        print("cbar text position not recognized")

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
    # Set either global (vtk) or local (src) bounds
    if preset.bounds == True:
        vmin, vmax = vtk.PointData.GetArray(0).GetRange(0)
    else:
        vmin, vmax = src.PointData.GetArray(0).GetRange(0)

    # Allow user to override global points with preset.bounds
    if preset.bounds:
        if isinstance(preset.bounds, (tuple, list)):
            # Allows user to only set one bound using e.g. [0, None]
            if preset.bounds[0] is not None:
                vmin = preset.bounds[0]
            if preset.bounds[1] is not None:
                vmax = preset.bounds[1]

    # Some fields should be centered on 0 despite the actual data bounds
    if preset.center:
        vabsmax = max((abs(vmin), abs(vmax)))
        vmin, vmax = -1 * vabsmax, vabsmax
    # If desired, round the colobar bounds to some base value
    if preset.rnd:
        vmin = myround(vmin, preset.rnd)
        vmax = myround(vmax, preset.rnd)

    # Apply depth specific values
    vsLUT.RescaleTransferFunction(vmin, vmax)

    # Set above and below color range on the colorbar
    if preset.above_range:
        vsLUT.UseAboveRangeColor = 1
        vsLUT.AboveRangeColor = preset.above_range
    if preset.below_range:
        vsLUT.UseBelowRangeColor = 1
        vsLUT.BelowRangeColor = preset.below_range

   
    # Manually set the values for the colorbar to always include min and max
    cbar = GetScalarBar(vsLUT)
    cbar.UseCustomLabels = 1
    custom_labels = [vmin, vmax]

    # Make sure that labels are set at equal division points between min and max
    divisions = preset.nlabel - 1
    dv = (vmax - vmin) / divisions
    for i in range(1, divisions):
        custom_labels.insert(-1, vmin + dv * i)

    cbar.CustomLabels = custom_labels

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


def delete_temp_objects(reg_names=None):
    """
    Delete objects that are made specifically for each screenshot, e.g. points,
    glyphs, rulers, text

    :type reg_names: list
    :param reg_names: names of registered objects that need to be deleted,
        can be parts of names such as 'ruler', 'point', 'glyph' etc.
    """
    if reg_names is None:
        return

    for src in [_[0] for _ in GetSources().keys()]:
        for name in reg_names:
            if name in src:
                src = FindSource(src)
                Delete(src)
                break


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


def set_camera(choice):
    """
    CCC
    Pre-defined camera locations for each type of plot, hard coded to the
    NZ North mesh.

    :type choice: str
    :param choice: choice of
        z: depth slice
        x: camera facing +y
        y: camera facing -x
        i: same as depth slice 'z' but 3D projection
        d: camera facing +x +y
        d_flip: camera facing -x +y
    :return:
    """
    ResetCamera()
    camera = GetActiveCamera()
    renderView = GetActiveView()
    renderView.InteractionMode = "2D"
    if choice == "z":
        # Hard set to the NZNorth model
        renderView.CameraPosition = [240, 280, 1600]
        renderView.CameraFocalPoint = [240, 280, 0]
        renderView.CameraViewUp = [0, 1, 0]
    elif choice == "x":
        renderView.CameraPosition = [-1100117., 5629267., -39303.]
        renderView.CameraFocalPoint = [402390.0, 5629267., -39303.]
        renderView.CameraParallelScale = 274961.
        renderView.CameraViewUp = [0.0, 0.0, 1.0]
    elif choice == "y":
        renderView.CameraPosition = [378617., 4093007., -47332.]
        renderView.CameraFocalPoint = [378617., 5595515.0, -47332.]
        renderView.CameraParallelScale = 300000.
        renderView.CameraViewUp = [0.0, 0.0, 1.0]
    elif choice == "i":
        renderView.InteractionMode = "3D"
        renderView.CameraPosition = [401399., 5567959., 1276057.]
        renderView.CameraFocalPoint = [401399., 5567959., -197500.0]
        renderView.CameraViewUp = [0, 1, 0]
    elif choice == "francis":
        renderView.CameraPosition = [84., 88., -27.]
        renderView.CameraFocalPoint = [516., 501., -27.]
        renderView.CameraViewUp = [0.0, 0.0, 1.0]
        renderView.CameraParallelScale = 250
    elif choice == "francis_chb":
        renderView.CameraPosition = [-220., 166., -50.]
        renderView.CameraFocalPoint = [1550., 778., -50.]
        renderView.CameraViewUp = [0.0, 0.0, 1.0]
        renderView.CameraParallelScale = 200
    elif choice == "tvz":
        renderView.CameraParallelScale = 300  # 500?
        renderView.CameraPosition = [1300., -415., -200.]
        renderView.CameraFocalPoint = [-300., 766., -200.]
        renderView.CameraViewUp = [0.0, 0.0, 1.0]
    elif choice == "cook":
        renderView.CameraPosition = [2000., -1600., -200.]
        renderView.CameraFocalPoint = [155., 122., -200.]
        renderView.CameraViewUp = [0.0, 0.0, 1.0]
        renderView.CameraParallelScale = 400
    elif choice == "cook_r":
        renderView.CameraPosition = [-1165., -623., -200.]
        renderView.CameraFocalPoint = [267., 252., -200.]
        renderView.CameraViewUp = [0.0, 0.0, 1.0]
        renderView.CameraParallelScale = 400
    elif choice == "along_strike":
        renderView.CameraPosition = [1683., -1170., -200.]
        renderView.CameraFocalPoint = [214., 298., -200.]
        renderView.CameraViewUp = [0.0, 0.0, 1.0]
        renderView.CameraParallelScale = 500
    elif choice in ["wellington_as", "flatpoint_as"]:
        renderView.CameraPosition = [-592., -766., -200.]
        renderView.CameraFocalPoint = [489., 518., -200.]
        renderView.CameraViewUp = [0.0, 0.0, 1.0]
        renderView.CameraParallelScale = 500
    elif choice == "Porangahau":
        renderView.CameraPosition = [-642., -931., -150.]
        renderView.CameraFocalPoint = [375., 276., -150.]
        renderView.CameraViewUp = [0.0, 0.0, 1.0]
        renderView.CameraParallelScale = 500 # 350
    elif choice == "Intraplate":
        renderView.CameraPosition = [3., -337., -150.]
        renderView.CameraFocalPoint = [318., 238., -150.]
        renderView.CameraViewUp = [0.0, 0.0, 1.0]
        renderView.CameraParallelScale = 250
    elif choice in ["Mahia", "Napier"]:
        renderView.CameraPosition = [-134., -102., -150.]
        renderView.CameraFocalPoint = [327., 446., -150.]
        renderView.CameraViewUp = [0.0, 0.0, 1.0]
        renderView.CameraParallelScale = 350  # 500?
    elif choice == "Seamounts":
        renderView.CameraPosition = [987., -296., -150.]
        renderView.CameraFocalPoint = [318., 269., -150.]
        renderView.CameraViewUp = [0.0, 0.0, 1.0]
        renderView.CameraParallelScale = 350
    elif choice == "tvz_strike":
        renderView.CameraPosition = [-689., -1531., -150.]
        renderView.CameraFocalPoint = [231., 308., -150.]
        renderView.CameraViewUp = [0.0, 0.0, 1.0]
    elif choice == "ruapehu_psf":
        renderView.CameraPosition = [-372.8, -380., -150.]
        renderView.CameraFocalPoint = [714., 909., -150.]
        renderView.CameraViewUp = [0.0, 0.0, 1.0]
    elif choice == "d_flip":
        # renderView.CameraPosition = [1422985., 4461623., -139631.]
        # renderView.CameraFocalPoint = [346611., 5537997., -139631.]
        renderView.CameraViewUp = [0.0, 0.0, 1.0]
        # renderView.CameraParallelScale = 325507.0
    else:
        sys.exit("Camera choice not found, check script")

    Render()


def plot_point_cloud(fid, reg_name="", point_size=2., color=None, opacity=1.,
                     xmin=None, xmax=None, ymin=None, ymax=None):
    """
    Plot a point cloud from an existing .VTK file, used for, e.g., coastlines
    shape outlines, rift zones etc.

    :type fid: str
    :param fid: file identifier of the .vtk file that should be used for
        point cloud plotting
    :type point_size: float
    :param float_size: size of the points representing the coastline
    :type color: str
    :param color: color of the coastline, defaults to COLOR
    """
    if color is None:
        color = rgb_colors("k")

    renderView = GetActiveView()
    point_cloud = OpenDataFile(fid)
    RenameSource(reg_name, point_cloud)
    Show(point_cloud, renderView)

    # Allow arbitrary clipping of the point cloud, e.g. to cut down coastline
    # to a specific region
    if any([xmin, xmax, ymin, ymax]):
        Hide(point_cloud)
        def clip_fx(vtk, origin, normal):
            clip_ = Clip(vtk)
            clip_.ClipType = "Plane"
            clip_.ClipType.Origin = origin
            clip_.ClipType.Normal = normal
            return clip_
        if xmin:
            point_cloud = clip_fx(point_cloud, origin=[xmin, 0., 0.], 
                                  normal=[-1, 0., 0.])
        if xmax:
            point_cloud = clip_fx(point_cloud, origin=[xmax, 0., 0.], 
                                  normal=[1, 0., 0.])
        if ymin:
            point_cloud = clip_fx(point_cloud, origin=[0., ymin, 0.], 
                                  normal=[0., -1., 0.])
        if ymax:
            point_cloud = clip_fx(point_cloud, origin=[0., ymax, 0.], 
                                  normal=[0., 1., 0.])


    RenameSource(reg_name, point_cloud)
    Show(point_cloud, renderView)
    pcDisplay = GetDisplayProperties(point_cloud, view=renderView)
    pcDisplay.SetScalarBarVisibility(renderView, False)
    pcDisplay.PointSize = point_size
    pcDisplay.Opacity = opacity
    pcDisplay.AmbientColor = color or LINECOLOR
    pcDisplay.DiffuseColor = color or LINECOLOR
    ColorBy(pcDisplay, None)

    return point_cloud


def plot_events(fid=None, color=None, type="Circle", size=3, position=None,
                linewidth=6., show="All Points", filled=1, opacity=1,
                representation_type="Surface With Edges", max_number=None):
    """
    Plot FOREST INVERSION event locations (VTK file) as 2D circles normal to
    the Z axis
    """
    renderView = GetActiveView()
    if fid is None:
        fid = ("/Users/Chow/Documents/academic/vuw/forest/utils/" 
               "vtk_files/scaled/srcs_d.vtk")
    srcs = OpenDataFile(fid)
    Hide(srcs)
    RenameSource("srcs", srcs)
    glyph = Glyph(registrationName="srcs_g", Input=srcs,
                  GlyphType="2D Glyph")
    if max_number is not None:
        glyph.MaximumNumberOfSamplePoints = max_number

    glyph.ScaleArray = ["POINTS", "No scale array"]
    glyph.ScaleFactor = size  # 10000 not scaled
    glyph.GlyphMode = show
    glyph.GlyphType.GlyphType = type
    glyph.GlyphType.Filled = filled
    glyphDisplay = Show(glyph, renderView, "GeometryRepresentation")
    glyphDisplay.PointSize = size
    glyphDisplay.SetRepresentationType(representation_type)
    glyphDisplay.LineWidth = linewidth
    glyphDisplay.Opacity = opacity
    glyph.GlyphTransform.Rotate = [0.0, 0.0, 45.0]
    if position is None:
        glyphDisplay.Position = [0.0, 0.0, 10.0]
    else:
        glyphDisplay.Position = position
 
    if color is None:
        vsLUT = GetColorTransferFunction("Z_Value")
        vsLUT.ApplyPreset("Inferno (matplotlib)")
    else:
        ColorBy(glyphDisplay, None)
        glyphDisplay.AmbientColor = color
        glyphDisplay.DiffuseColor = color


def plot_stations(fid=None, type="Triangle", size=3, position=None,
                  linewidth=6., show="All Points", filled=1, opacity=1,
                  representation_type="Surface With Edges"):
    """
    Plot FOREST inversion station locations (VTK file) as 2D diamonds normal
    to the Z axis
    """
    renderView = GetActiveView()
    if fid is None:
        fid = ("/Users/Chow/Documents/academic/vuw/forest/utils/"
               "vtk_files/scaled/rcvs.vtk")
    rcvs = OpenDataFile(fid)
    RenameSource("rcvs", rcvs)
    glyph = Glyph(registrationName="rcvs_g", Input=rcvs,
                  GlyphType="2D Glyph")
    glyph.ScaleArray = ["POINTS", "No scale array"]
    glyph.ScaleFactor = size  # 10000 not scaled
    glyph.GlyphMode = show
    glyph.GlyphType.GlyphType = type
    glyph.GlyphTransform.Rotate = [0.0, 0.0, 180.0]
    glyph.GlyphType.Filled = filled
    glyphDisplay = Show(glyph, renderView, "GeometryRepresentation")
    glyphDisplay.PointSize = size
    glyphDisplay.SetRepresentationType(representation_type)
    glyphDisplay.LineWidth = linewidth
    glyphDisplay.Opacity = opacity
    if position is None:
        glyphDisplay.Position = [0.0, 0.0, 10.0]
    else:
        glyphDisplay.Position = position


def plot_sse_slip_patches():
    """
    Plot the slip patches for slow slip events provided by Laura Wallace,
    converted into VTK files for convenience
    """
    renderView = GetActiveView()
    sse = OpenDataFile("/Users/Chow/Documents/academic/vuw/forest/utils/"
                      "vtk_files/sse.vtk")
    RenameSource("sse", sse)
    Show(sse, renderView)
    csLUT = GetColorTransferFunction("cumulative_slip")
    csLUT.ApplyPreset("BrOrYl", True)
    csLUT.RescaleTransferFunction(50., 450.)
    csLUT.EnableOpacityMapping = 1  # lower values have higher opacity
    csPWF = GetOpacityTransferFunction("cumulative_slip")
    csPWF.UseLogScale = 1
    vtkDisplay = GetDisplayProperties(sse, view=renderView)
    vtkDisplay.SetScalarBarVisibility(renderView, False)
    vtkDisplay.Opacity = .75


def plot_srvtk(fid, src_depth=3E3):
    """
    For plotting sensitivity kernels, where only one source and one receiver
    are necessary, it's assumed that these files are produced by xspecfem3d,
    that there are only two points in the file, and that the first point is
    the source, the second the receiver
    """
    renderView = GetActiveView()
    lines = open(fid, "r").readlines()
    src = [float(_) for _ in lines[5].strip().split()]
    if src_depth is not None:
        src[-1] = src_depth
    rcv = [float(_) for _ in lines[6].strip().split()]

    # Glyph the source as a green sphere
    pointSource = PointSource(registrationName="PointSource1")
    pointSource.Center = src
    glyph = Glyph(registrationName="Glyph1", Input=pointSource)
    glyph.ScaleFactor = 12000
    glyph.GlyphType = "Sphere"
    glyphDisplay = Show(glyph, renderView, "GeometryRepresentation")
    glyphDisplay.AmbientColor = rgb_colors("g")
    glyphDisplay.DiffuseColor = rgb_colors("g")

    # Glyph the receiver as a white box
    pointSource = PointSource(registrationName="PointSource2")
    pointSource.Center = rcv
    glyph = Glyph(registrationName="Glyph2", Input=pointSource)
    glyph.ScaleFactor = 12000
    glyph.GlyphType = "Box"
    glyphDisplay = Show(glyph, renderView, "GeometryRepresentation")
    glyphDisplay.AmbientColor = rgb_colors("w")
    glyphDisplay.DiffuseColor = rgb_colors("w")
    glyphDisplay.SetRepresentationType('Surface With Edges')


def plot_point(origin, size, rotate=None):
    """
    Plot a point source and 2D glyph
    """
    renderView = GetActiveView()

    ps = PointSource(registrationName="PointSource1")
    ps.Center = origin
    glyph = Glyph(registrationName="Glyph1", Input=ps, 
                  GlyphType="2D Glyph")
    glyph.ScaleFactor = size
    glyph.GlyphType.GlyphType = "Circle"
    if rotate:
        glyph.GlyphTransform.Rotate = rotate


    glyphDisplay = Show(glyph, renderView, "GeometryRepresentation")

    glyphDisplay.AmbientColor = rgb_colors("g")
    glyphDisplay.DiffuseColor = rgb_colors("g")
    glyphDisplay.LineWidth = 4.


def plot_landmarks(src, color=None, glyph_type="Cross", size=10):
    """
    Plot landmark locations as glyphs, src needs to be a dictionary
    with values as lists of floats [x, y, z] denoting the location of the point
    :type src: dict
    """
    if color is None:
        color = LINECOLOR

    renderView = GetActiveView()

    for i, (name, origin) in enumerate(src.items()):
        # Glyph the source as a green sphere
        pointSource = PointSource(registrationName=f"PointSource{i}")
        pointSource.Center = origin
        glyph = Glyph(registrationName=f"Glyph{i}", Input=pointSource,
                      GlyphType="2D Glyph")
        glyph.ScaleArray = ["POINTS", "No scale array"]
        glyph.ScaleFactor = size
        glyph.GlyphMode = "All Points"
        glyph.GlyphType.GlyphType = glyph_type
        glyph.GlyphType.Filled = 1
        glyphDisplay = Show(glyph, renderView, "GeometryRepresentation")
        glyphDisplay.PointSize = 3.
        glyph.GlyphTransform.Rotate = [0.0, 0.0, 0.0]
        glyphDisplay.SetRepresentationType("Surface")
        glyphDisplay.LineWidth = 3.0
        glyphDisplay.Opacity = 1.
        glyphDisplay.AmbientColor = color
        glyphDisplay.DiffuseColor = color


def plot_checkers(src, slicez=-3000.):
    """
    Plot checkerboard as an overlay to see how the point spread function
    performs. Plots only one slice through the checkers and uses contour and
    a lower opacity so that the underlying model is still visible through
    the overlay

    :return:
    """
    renderView = GetActiveViewOrCreate('RenderView')
    checkers = OpenDataFile(src)
    RenameSource("checkers", checkers)
    checkersDisplay = Show(checkers, renderView,
                           "UnstructuredGridRepresentation")
    checkersDisplay.SetScalarBarVisibility(renderView, False)
    Hide(checkers, renderView)

    # Make a slice as we only want an overlay
    slice = Slice(registrationName="SliceC", Input=checkers)
    slice.SliceType.Origin = [402390.0, 5595515.0, slicez]
    slice.SliceType.Normal = [0., 0., 1.]
    sliceDisplay = Show(slice, renderView, "GeometryRepresentation")
    sliceDisplay.SetScalarBarVisibility(renderView, False)
    Hide(slice, renderView)

    # Make a contour so we can see through the checkers to the underlying model
    contour = Contour(registrationName="ContourC", Input=slice)
    contourDisplay = Show(contour, renderView, "GeometryRepresentation")
    contourDisplay.Position = [0.0, 0.0, 5000.0]
    contourDisplay.LineWidth = 1.
    contourDisplay.DataAxesGrid.Position = [0.0, 0.0, 5000.0]
    # contour.Isosurfaces = [1499.0, 1499.2, 1499.4, 1499.6, 1499.8, 1500.2,
    #                        1500.4, 1500.6, 1500.8, 1501.0]
    contour.Isosurfaces =[-200, -180, -160, -140, -120, -100, -80, -60, -40,
                          40, 60, 80, 100, 120, 140, 160, 180, 200]


def plot_cross_section_surface_trace(src, normal, origin, flip_view, color="w",
                                     tick_spacing_km=50, line_width=10.):
    """
    For depth slices, it's useful to visualize the surface trace of a cross
    section for quick reference on how the two views relate. This function
    takes the normal and origin that define the cross section and creates
    a ruler that matches the horizontal axis ruler created by
    create_ruler_grid_axes_cross_section()

    :param normal:
    :param origin:
    :return:
    """
    # Assuming that the src is a plane in XY, make a slice to figure out
    # the coordinates of the plane
    slice_vtk = cross_section(src, normal=normal, origin=origin, name="tmpslc")
    Hide(slice_vtk)

    x, y, z = get_coordinates(slice_vtk)

    if flip_view:
        ruler_origin = [min(x), min(y), 1]
    else:
        ruler_origin = [min(x), max(y), 1]

    # Use trig to figure out the length of the angled cut
    slice_length_m = ((max(y) - min(y)) ** 2 + (max(x) - min(x)) ** 2) ** .5

    # Calculate the end points of the hypotenuse based on the side lengths
    angle_rad = math.atan(normal[0] / normal[1])
    ruler_h = [ruler_origin[0] + slice_length_m * math.cos(angle_rad),
               ruler_origin[1] - slice_length_m * math.sin(angle_rad),
               ruler_origin[2]]

    # Create the 'X' axis ruler along the bottom
    create_ruler(point1=ruler_origin, point2=ruler_h,
                 graduation=tick_spacing_km,
                 reg_name="ruler1",  axis_color=rgb_colors(color),
                 line_width=line_width)

    Delete(slice_vtk)


def mark_point(origin, depth=5000, c="g"):
    """
    Mark a point on depth slices for reference
    """
    # Glyph the source as a green sphere
    pointSource = PointSource(registrationName="PointSource")
    pointSource.Center = [origin[0], origin[1], 5000]
    glyph = Glyph(registrationName=f"Glyph", Input=pointSource,
                  GlyphType="2D Glyph")
    glyph.ScaleArray = ["POINTS", "No scale array"]
    glyph.ScaleFactor = 10000.
    glyph.GlyphMode = "All Points"
    glyph.GlyphType.GlyphType = "Circle"
    glyph.GlyphType.Filled = 1
    glyphDisplay = Show(glyph, renderView, "GeometryRepresentation")
    glyphDisplay.PointSize = 3.
    glyph.GlyphTransform.Rotate = [0.0, 0.0, 180.0]
    glyphDisplay.SetRepresentationType("Surface With Edges")
    glyphDisplay.LineWidth = 3.0
    glyphDisplay.Opacity = 1.
    glyphDisplay.AmbientColor = rgb_colors(c)
    glyphDisplay.DiffuseColor = rgb_colors(c)


def label(s):
    """
    Add a text label to depth slices for publication figures, e.g. 'A' or '1'
    Will be put in the top left corner
    """
    create_text(s=f"[{s}]", position=[.2, .85], reg_name="text1",
                fontsize=int(FONTSIZE*3), color=rgb_colors("w"),
                justification="center")


def features():
    """
    Hacked together, text labels on figures showing locations of features
    """
    features = {"A": [.67, .585],  # Mahia
                "B": [.54, .43],  # Pora
                "C": [.34, .235],  # Cook
                "D": [.51, .67],  # TVZ
                "E": [.61, .4]}  # PORA
    for label, position in features.items():
        create_text(s=label, position=position, reg_name="text1",
            fontsize=int(FONTSIZE*1.5), color=rgb_colors("w"),
            justification="center")


def plot_2D_surface_projection(fid, origin, normal, flip_view,
                               scale=1, offset=1E3, pointsize=1.5, color=None,
                               position=None, **kwargs):
    """
    A hacky method of projecting a 3D surface onto a 2D plane.

    For example, we have a 3D plate interface model from Williams et al. (2013),
    and when we make a cross section we want to plot the line that defines this
    3D model on our 2D plane.

    This works but clipping a small buffer on either side of the 3D model,
    oriented in the same direction as the cross sections plane, essentially
    leaving a small 3D strip. If the strip is small enough, it should look like
    a line, and be visible on either side of the cross section. Vaya!

    :type origin: list of float
    :param origin: origin points for the xsection
    :type normal: list of float
    :param normal: normal vector of the plane
    :type offset: float
    :param offset: amount to offset x and y to generate the strip. The larger
        the value, the larger the strip, and probably the more 3D variation,
        which will look weird on a 2D plot. You may need to play around with
        this value to get the interface to look like a line rather than a
        strip
    :type pointsize: float
    :param pointsize: control the size of the points that plot the original
        interface. This along with offset may help you get a more line-looking
        interface
    :type flip_view: bool
    :param flip_view: normal viwe axis in the +x+y quadrant, flipping the view
        shifts the view angle 90 deg CCW to the -x+y quadarant
    :type color: list of float
    :param color: color of the interface, in rgb
    """
    if color is None:
        color = rgb_colors("w")

    renderView = GetActiveViewOrCreate('RenderView')
    vtk = OpenDataFile(fid)

    # Clip unused parts of the volume
    if kwargs:
        vtk = clip_domain(vtk, **kwargs)
        Hide(vtk)

    # Now clip away one side of the cross section
    # Too lazy to do this with trig so just add some km to x and y and pray!
    clip_1 = Clip(vtk)
    clip_1.ClipType = "Plane"
    if flip_view:
        clip_1.ClipType.Origin = [origin[0] - offset,
                                  origin[1] + offset,
                                  origin[2]]
    else:
        clip_1.ClipType.Origin = [origin[0] - offset,
                                  origin[1] - offset,
                                  origin[2]]

    clip_1.ClipType.Normal = normal

    # Now clip away the other side of the cross section
    clip_2 = Clip(clip_1)
    clip_2.ClipType = "Plane"
    if flip_view:
        clip_2.ClipType.Origin = [origin[0] + offset,
                                  origin[1] - offset,
                                  origin[2]]
    else:
        clip_2.ClipType.Origin = [origin[0] + offset,
                                  origin[1] + offset,
                                  origin[2]]

    clip_2.ClipType.Normal = [_ * -1 for _ in normal]

    clipDisplay = Show(clip_2, renderView, "UnstructuredGridRepresentation")
    clipDisplay.Scale = [1., 1., scale]
    clipDisplay.PointSize = pointsize
    if position:
        clipDisplay.Position = position

    try:
        ColorBy(clipDisplay, None)
        clipDisplay.AmbientColor = color
        clipDisplay.DiffuseColor = color
        return clip_2
    except RuntimeError:
        print("2D projection has clipped away the visible bounds\n"
              f"will not plot: {fid}")
        return None


def depth_slices(fid, slices, preset, save_path=os.getcwd()):
    """
    Main function for creating and screenshotting depth slices (slice plane
    normal/perpendicular to Z axis). Creates a standardized look for the
    depth slices with scale bar, colorbar and annotations
    """
    args = parse_args()

    # Open the model volume
    vtk = OpenDataFile(fid)
    vtk = scale_input(vtk, scale_units=preset.scale_units, 
                      scale_coords=args.coord_scale)

    RenameSource("surface", vtk)
    set_camera("z")
    Show(vtk)
    if not "surface" in slices:
        Hide(vtk)

    # Annotate file name and depth in title
    text, _ = create_text(s=f"", position=[0.19, 0.94], 
                          reg_name="text1", fontsize=int(FONTSIZE * 1.25), 
                          color=FONTCOLOR, justification="center")

    # Annotate text for X and Y labels
    create_text(s=f"X [km]", position=[0.45, 0.075], reg_name="text1", 
                fontsize=int(FONTSIZE), color=FONTCOLOR, 
                justification="center")
    create_text(s=f"Y [km]", position=[0.075, 0.575], reg_name="text1",
            fontsize=int(FONTSIZE), color=FONTCOLOR,
            justification="center")

    # Create bounding axes with pre-defined tick marks using rulers
    create_ruler_grid_axes_depth_slice(vtk)


    # Sit the colorbar partway down from the top left corner
    if not args.cbar_orientation or \
            args.cbar_orientation.lower() in ["horizontal", "h"]:
        vsLUT, cbar = set_colormap_colorbar(vtk,
                                            # position=[0.83, 0.6],
                                            position=[0.6, 0.07],
                                            length=.15, 
                                            orientation="Horizontal", 
                                            thickness=70, 
                                            text_position="top")
    elif args.cbar_orientation.lower in ["vertical", "v"]:
        vsLUT, cbar = set_colormap_colorbar(vtk, position=[0.78, 0.6], 
                                            length=.2, orientation="Vertical", 
                                            thickness=80)

    if args.surface_traces:
        for name in args.surface_traces:
            try:
                origin, normal = DIAGONALS[name.title()]
            except KeyError:
                print(f"{name} not in DIAGONALS")
                continue
            plot_cross_section_surface_trace(vtk, normal, origin,
                                            flip_view=bool(normal[1] > 0))

    # Plot an origin point on the map
    if args.mark:
        mark_point(DIAGONALS[args.mark.title()][0])

    # Generate slices through the volume at desired depth levels
    # Make special precautions if we're looking at surface projections
    for slice_ in slices:
        if slice_ == "surface":
            tag = "z_00surf"
            slice_vtk = vtk
            text.Text = f"{TITLE} [{slice_}]"

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
            if args.contour:
                contour_lines(slice_vtk, preset)
            if args.minmax:
                create_minmax_glyphs(slice_vtk)
            tag = f"z_{slice_:0>2}km"
            text.Text = f"{TITLE} [{tag.replace('_', '=')}]"
            rescale_colorscale(vsLUT, src=slice_vtk, vtk=vtk, preset=preset)

        if not args.cbar_off:
            show_colorbar(slice_vtk)


        set_depth_slice_data_axis_grid(slice_vtk)

        # Save screenshot, hide slice and move on, job done
        fid_out = os.path.join(save_path, f"{tag}.png")
        SaveScreenshot(fid_out, renderView, ImageResolution=VIEW_SIZE,
                       TransparentBackground=1)
        remove_whitespace(fid_out)
        Hide(slice_vtk, renderView)
        
        # Delete the min max value points because they'll change w/ each slice
        delete_temp_objects(reg_names=["glyph", "point", "contour"])


def cross_sections(fid, percentages, normal, preset, save_path=os.getcwd()):
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
    args = parse_args()
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
        naxis = x
        parallel = "y"
        nvector = [1., 0., 0.]
    elif normal == "y":
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
        clip_vtk = clip_bottom(slice_vtk, abs(args.depth_km) * -1E3)
        Show(clip_vtk, renderView, "UnstructuredGridRepresentation")

        if args.contour:
            contour_lines(clip_vtk, preset)

        scale_vertical_grid(clip_vtk, camera_parallel=False,
                                    min_depth_km=0., max_depth_km=args.depth_km)

        # Reset the colorbounds to data range
        active = GetActiveSource()
        display = GetDisplayProperties(active, view=renderView)
        display.RescaleTransferFunctionToDataRange(False, True)

        # Use rulers to generate scale bars
        create_ruler_grid_axes_cross_section(src=clip_vtk, normal=normal,
                                             xlabel=parallel, scale=args.zscale,
                                             depth_cutoff_km=args.depth_km,
                                             dz_km=args.dz)

        # Annotate distance
        anno = f"{normal.upper()}={(axis_dist-min(naxis))*1E-3:.2f}km"
        create_text(f"{anno}\n{TITLE}", color=FONTCOLOR,
                    position=[0.3, args.anno], reg_name="text1",
                    fontsize=int(FONTSIZE))

        # Generate and rescale the colorbar/ colormap
        vsLUT, cbar = set_colormap_colorbar(vtk,
                                            position=[0.5, args.anno + 0.01],
                                            orientation="Horizontal",
                                            thickness=25, length=.1)

        cbar.TextPosition = "Ticks left/bottom, annotations right/top"
        rescale_colorscale(vsLUT, src=clip_vtk, vtk=vtk, preset=preset)
        display = GetDisplayProperties(clip_vtk, view=renderView)
        display.SetScalarBarVisibility(renderView, not args.cbar_off)

        # Reset camera view to be normal to the plane. Specific to this plane
        set_camera(normal)


        fid_out = os.path.join(save_path, f"{tag}.png")
        SaveScreenshot(fid_out, renderView, ImageResolution=VIEW_SIZE,
                       TransparentBackground=1)
        remove_whitespace(fid_out)

        # Clean up for next plot
        Hide(clip_vtk, renderView)
        delete_temp_objects(reg_names=["ruler", "text", "contour"])


def interface(fid, preset, save_path=os.getcwd()):
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
    set_camera("i")

    # Standard accoutrements, colorbar, min/max values
    vsLUT, cbar = set_colormap_colorbar(interpolator, position=[0.249, 0.75],
                                        orientation="Vertical",)
    rescale_colorscale(vsLUT, src=interpolator, vtk=vtk, preset=preset)
    show_colorbar(interpolator)


    # Annotate the bottom-right corner to explain this is the interface
    text, _ = create_text(s="Interface", position=[0.625, 0.1], color=FONTCOLOR,
                          reg_name="text1", fontsize=FONTSIZE * 2)

    # Annotate the file ID so we know what we're plotting
    create_text(s=TITLE, color=FONTCOLOR, position=[0.25, 0.95], 
                reg_name="text2", fontsize=FONTSIZE * 2)


    # Save screenshot, hide the surface and move on
    fid_out = os.path.join(save_path, f"interface.png")
    SaveScreenshot(fid_out, renderView, ImageResolution=VIEW_SIZE,
                   TransparentBackground=1)
    remove_whitespace(fid_out)

    Hide(interpolator, renderView)
    

def diagonal(fid, preset, name=None, save_path=os.getcwd()):
    """
    Create a diagonal cross section through the domain that cuts through an
    arbitrary plane defined by an origin point and normal vector. Allows
    choosing of the viewing direction.

    :type fid: str
    :param fid: name of the file to cut diagonally
    :type preset: dict
    :param preset: colormap controls
    :type name: str
    :param name: name of the origin point for annotation
    :type save_path: str
    :param save_path: path to save figure

    .. note::
        flip_view will change view angle from +x+y to -x+y

        -x+y quadrant, e.g.

         -x+y | +x+y
         _____|_____
              |
              |
    """
    args = parse_args()

    # Normal is defined as 40deg from X-axis following Eberhart-Phillips
    if name is None:
        origin = args.origin
        normal = args.normal
        assert(origin is not None and normal is not None), \
            "if pre-defined name is not set, origin and normal must be given"
    else:
        try:
            origin, normal = DIAGONALS[name.title()]
        except KeyError:
            sys.exit(f"pre-defined origin {name} does not exist")

    tag = f"diagonal_{name.lower()}"
    # Decide the view angle based on the normal vector
    flip_view = bool(normal[1] > 0)

    vtk = OpenDataFile(fid)
    # vtk = clip_domain(vtk, xmax=255, ymax=255)
    RenameSource("surface", vtk)
    Show(vtk)
    ResetCamera()
    renderView = GetActiveView()
    renderView.InteractionMode = "2D"
    Hide(vtk, renderView)

    # Slice across the given cross section plane
    slice_vtk = cross_section(vtk=vtk, normal=normal, origin=origin, name=name)
    Hide(slice_vtk, renderView)

    # Cut off depths below a certain range as we are not interested in deep
    vtks = []
    clip_vtk = clip_domain(vtk, zmin=-1 * abs(args.depth_km))
    Show(clip_vtk, renderView, "UnstructuredGridRepresentation")
    vtks.append(clip_vtk)

    if args.contour:
        contour_lines(clip_vtk, preset, scale=args.zscale)

    # Reset the colorbounds to data range
    active = GetActiveSource()
    display = GetDisplayProperties(active, view=renderView)
    display.RescaleTransferFunctionToDataRange(False, True)


    # Use Paraviews grid lines to show depth values only, as arbitrary cross
    # sections will show the 3D outline of the volume on a 2D plane which
    # is confusing/ not helpful
    # !!! Fix this
    # scale_vertical_grid(src, camera_parallel=False, dz_km=dz_km,
    #                             min_depth_km= min(z) + dist_m_v,
    #                             max_depth_km=depth_cutoff_km,
    #                             scale=scale)
    # Use rulers to define the axis grids w/ tickmarks etc.
    # create_ruler_grid_axes_cross_section(src=clip_vtk, normal=normal,
    #                                      xlabel="H", scale=args.zscale,
    #                                      depth_cutoff_km=args.depth_km,
    #                                      dz_km=args.dz, flip_view=flip_view)

    # Generate and rescale the colorbar/ colormap
    vsLUT, cbar = set_colormap_colorbar(vtk, position=[0.5, args.anno + 0.01],
                                        orientation="Horizontal",
                                        thickness=80, length=.15)

    cbar.TextPosition = "Ticks left/bottom, annotations right/top"
    rescale_colorscale(vsLUT, src=clip_vtk, vtk=vtk, preset=preset)
    display = GetDisplayProperties(clip_vtk, view=renderView)
    display.SetScalarBarVisibility(renderView, not args.cbar_off)

    # Create a reference point based on the landmark origin location
    # mark_surface_point(origin, normal)

    # Plot the interface model of Williams et al. (2013)
    if args.williams:
        fid = ("/Users/Chow/Documents/academic/vuw/"
               "forest/utils/vtk_files/interface.vtk")
        vtks.append(plot_2D_surface_projection(fid, origin, normal, flip_view))

    # Create a mark wherever the coastline intersects this plane
    if args.outline:
        fid = ("/Users/Chow/Documents/academic/vuw/"
               "forest/utils/vtk_files/coast.vtk")
        vtks.append(plot_2D_surface_projection(fid, origin, normal, flip_view,
                                               offset=1E3, pointsize=4.5,
                                               color=rgb_colors("k"))
                    )

    # Annotate landmark location text and filename
    if not args.anno_off:
        create_text(f"{name} {TITLE}", [0.2, args.anno], color=FONTCOLOR,
                    reg_name="text1", fontsize=int(FONTSIZE))

    # Reset camera view to be normal to the plane. Specific to this plane
    if flip_view:
        set_camera("d_flip")
    else:
        set_camera("d")


    fid_out = os.path.join(save_path, f"{tag}.png")
    SaveScreenshot(fid_out, renderView, ImageResolution=VIEW_SIZE,
                   TransparentBackground=1)
    remove_whitespace(fid_out)

    # Clean up for next plot
    for vtk_ in vtks:
        Hide(vtk_, renderView)
    delete_temp_objects(reg_names=["ruler", "glyph", "point", "text",
                                   "contour"])

def make_preplot(args):
    """
    Convenience function to plot extras such as coastline, srcs, rcvs
    """
    util_dir = "/Users/Chow/Documents/academic/vuw/forest/utils/vtk_files/"
    if args.coord_scale != 1:
        util_dir = os.path.join(util_dir, "scaled")
    if args.outline:
        if args.verbose:
            print(f"\t\tPlotting coastline")
        plot_point_cloud(fid=os.path.join(util_dir, "coast.vtk"),
                         reg_name="coast", color=LINECOLOR, point_size=5)
        if args.verbose:
            print(f"\t\tPlotting Lake Taupo outline")
        plot_point_cloud(fid=os.path.join(util_dir, "taupo.vtk"),
                         reg_name="taupo", color=LINECOLOR, point_size=5)
    if args.sources:
        if args.verbose:
            print(f"\t\tPlotting event glyphs")
        plot_events(color=rgb_colors("g"))
    if args.receivers:
        if args.verbose:
            print(f"\t\tPlotting receiver glyphs")
        plot_stations()
    if args.faults:
        if args.verbose:
            print(f"\t\tPlotting active fault traces")
        plot_point_cloud(fid=os.path.join(util_dir, "faults.vtk"),
                         reg_name="faults", point_size=2., opacity=1.,
                         color=LINECOLOR)
    if args.slowslip:
        if args.verbose:
            print(f"\t\tPlotting SSE slip patches")
        plot_sse_slip_patches()
    if args.srvtk:
        plot_srvtk(fid=args.srvtk)
    if args.landmarks:
        plot_landmarks(color=rgb_colors("k"))
    if args.checkers:
        if isinstance(args.checkers, str):
            checkers = args.checkers
        else:
            checkers = ("/Users/Chow/Documents/academic/vuw/tomo/" 
                        "point_spread_test/dvs_stagger/checkers/"
                        "psf_stagger_vs.vtk")
        # plot_checkers(src=os.path.join(util_dir, "checkers.vtk"))
        plot_checkers(src=checkers)
    if args.williams:
        plot_point_cloud(fid=os.path.join(util_dir, "interface_contours.vtk"),
                         reg_name="intcont", point_size=3, opacity=.6,
                         color=rgb_colors("w"))
        # Super hard-coded locations for depth slices only, will not line up
        # if view size or camera location is changed
        contour_lines = {"3": [.71, .31], "6": [.68, .31], "9": [.64, .31],
                         "12": [.558, .31], "15": [.49, .31], "20": [.45, .31],
                         "30": [.39, .31], "40": [.35, .31], "50": [.315, .31],
                         "75": [.27, .31], "100": [.24, .31]}
        for s, pos in contour_lines.items():
            create_text(s=s, position=pos, fontsize=int(FONTSIZE / 2), 
                        color=rgb_colors("w"))
    if args.features:
        features()

    if args.label is not None:
        label(args.label)


# HACKED IN SPECIFIC FIGURE MAKERS
def manual_depth_slice(fid, preset, save_path=os.getcwd()):
    """
    ZZZ
    Main function for creating and screenshotting depth slices (slice plane
    normal/perpendicular to Z axis). Creates a standardized look for the
    depth slices with scale bar, colorbar and annotations
    DOMAIN BOUNDS:
    x: [0, 462.156005859375]
    y: [0, 617.1300048828125]
    z: [-400., 2.2788798809051514]
    """
    args = parse_args()
    choice = args.choice
    pert_origin = None

    if choice in ["porangahau", "mahia", "intraplate", "seamounts"]:
        xmin=175
        xmax=462.156
        ymin=100
        ymax=560
        depth_km = int(args.depth_km)
        camera_position = [(xmax+xmin)/2, (ymax+ymin)/2-20, 600]
        focal_point =[(xmax+xmin)/2, (ymax+ymin)/2-20, 0]
        title = [0.28, 0.875]
        xlabel = [0.425, 0.125]
        ylabel = [0.19, 0.53]
        cbar = [0.55, 0.125]  # PSF .6 .1
        cbar_length = .135
        cbar_thickness = 60
        FONTSIZE = 50
        spacing_km=100
        round_to=100
        ruler_tick_space=25

        # cross_sections = []
        #cross_sections = [choice.title()]
        cross_sections = ["Mahia_Psf"]
        cross_sections = ["Mahia", "Porangahau", "Seamounts"]
        earthquakes = False
        faults=False
        cbar_justification="Right"

        # Perturbations
        if choice == "porangahau":
            pert_origin = [295.54, 251.911, -11]
            pert_size = 8.85
        elif choice == "mahia":
            pert_origin = [406.68, 381.1, -11]
            pert_size = 21
        elif choice == "intraplate":
            pert_origin = [337.8, 228.1, -16]
            pert_size = 21
    elif choice in ["volumetric", "updates"]:
        xmin=0
        xmax=462.156
        ymin=0
        ymax=616.130
        depth_km = int(args.depth_km)
        camera_position = [(xmax+xmin)/2, (ymax+ymin)/2-20, 600]
        focal_point =[(xmax+xmin)/2, (ymax+ymin)/2-20, 0]
        title = [0.23, 0.88]
        xlabel = [0.475, 0.11]
        ylabel = [0.125, 0.53]
        cbar = [0.6, 0.1]  # PSF .6 .1
        cbar_length = .135
        cbar_thickness = 60
        # cbar_length = .2
        # cbar_thickness = 80
        FONTSIZE = 50
        spacing_km=100
        round_to=100
        ruler_tick_space=25
        cross_sections = []
        # cross_sections = ["Csfaults", "Cspert_R"]
        cross_sections = ["Reyners", "Porangahau"]
        # cross_sections = []
        earthquakes = False
        faults=True
        cbar_justification="Right"
    elif choice == "along_strike":
        xmin=0
        xmax=462.156
        ymin=0
        ymax=616.130
        depth_km = 10
        camera_position = [(xmax+xmin)/2, (ymax+ymin)/2-20, 600]
        focal_point =[(xmax+xmin)/2, (ymax+ymin)/2-20, 0]
        title = [0.23, 0.88]
        xlabel = [0.475, 0.11]
        ylabel = [0.125, 0.53]
        cbar = [0.4, 0.1]  # PSF .6 .1
        # cbar_length = .135
        # cbar_thickness = 60
        cbar_length = .5
        cbar_thickness = 80
        FONTSIZE = 50
        spacing_km=100
        round_to=100
        ruler_tick_space=25
        cross_sections = ["Reyners", "Wellington"]#, "Flatpoint"]
        # cross_sections = ["Csfaults", "Cspert_R"]
        # cross_sections = []
        earthquakes = False
        faults=True
        cbar_justification="Right"
    elif choice == "cook":
        xmin=0
        xmax=255
        ymin=0
        ymax=205
        depth_km = int(args.depth_km)
        camera_position = [(xmax+xmin)/2, (ymax+ymin)/2-20, 600]
        focal_point =[(xmax+xmin)/2, (ymax+ymin)/2-20, 0]
        title = [0.25, 0.74]
        xlabel = [0.475, 0.26]
        ylabel = [0.125, 0.505]
        cbar = [0.6, 0.26]  # PSF
        cbar_length = .135
        cbar_thickness = 60
        FONTSIZE = 50
        spacing_km=100
        round_to=100
        ruler_tick_space=25
        cross_sections = []
        # cross_sections = ["psf_cook", "psf_cook_r", "Csfaults"]
        # cross_sections = ["Csfaults_R", "Csfaults"]
        earthquakes = True
        faults = True
        cbar_justification="Right"
        pert_origin = [145.8, 75.18, -3]
        pert_size = 28
    elif choice == "tvz":
        # xmin=170
        # xmax=340
        # ymin=250
        # ymax=600
        xmin=50 # 50
        xmax=400 # 350
        ymin=325
        ymax=575
        depth_km = int(args.depth_km)
        camera_position = [200, 400, 600]
        focal_point = [200, 400, -100]
        title = [0.235, 0.8]
        xlabel = [0.5, 0.275]
        ylabel = [0.125, 0.525]
        cbar = [0.65, 0.275]
        # cbar = [0.6, 0.22]  # PSF
        cbar_length = .135
        cbar_thickness = 60
        FONTSIZE = 50
        spacing_km=100
        round_to=100
        ruler_tick_space=25
        volcanos = {'Ruapehu': [205.222, 364.035, 0],
                    'Taranaki': [75., 359., 0.],
                    'Tarawera': [280.530, 484.901, 0],
                    'White Island': [344.691, 560.552, 0],
                    'Tongariro': [212, 380, 0.],
                    # 'Ngauruhoe': [205, 363, 0.],
                    # 'Taupo': [247.81900000000002, 430.741, 0],
                    # 'Rotorua': [262.539, 491.572, 0],
                    # 'Reparoa': [270., 460., 0.],
                    # 'Whakatane': [327.409, 511.833, 0],
                    # 'Putauaki': [311.152, 496.12, 0],
                     }
        # cross_sections = ["tvz_rift", "tvz_rift2", "tvz_rift3"]
        cross_sections = ["tvzpsf"]
        cbar_justification="Left"

        pert_origin = [269.99, 489.493, -3]
        pert_origin = [269, 489, -5]
        pert_size = 20
    elif choice == "francis":
        xmin=250
        xmax=455
        ymin=300
        ymax=505
        depth_km = 2
        cross_sections = []  # ["mahia", "centhawkbay"]
        camera_position = [400, 400, 300]
        focal_point = [400, 400, -100]
        title = [0.2, 0.71]
        xlabel = [0.35, 0.24]
        ylabel = [0.1, 0.525]
        cbar = [0.46, 0.22]
        cbar_length = .135
        cbar_thickness = 60
        FONTSIZE = 50
        taupo=False
        spacing_km=100
        round_to=50
        ruler_tick_space=50
        cbar_justification="Left"
    elif choice == "ruapehu_psf":
        xmin=50
        xmax=350
        ymin=200
        ymax=450
        depth_km = int(args.depth_km)
        camera_position = [(xmax+xmin)/2, (ymax+ymin)/2-20, 600]
        focal_point =[(xmax+xmin)/2, (ymax+ymin)/2-20, 0]
        title = [0.23, 0.76]
        xlabel = [0.475, 0.25]
        ylabel = [0.125, 0.53]
        cbar = [0.6, xlabel[1]]  # PSF .6 .1
        cbar_length = .135
        cbar_thickness = 60
        FONTSIZE = 50
        spacing_km=100
        round_to=100
        ruler_tick_space=50
        cross_sections = []
        earthquakes = False
        faults=False
        cbar_justification="Right"

        pert_origin = [207.6, 365.88, -12]
        pert_size = 10.5
        cross_sections = ["Ruapehu_Psf"]
    else:
        print("unrecognized choice, using default")
        # Default Z plots
        xmin=0
        xmax=462.156
        ymin=0
        ymax=616.130
        depth_km = int(args.depth_km)
        camera_position = [(xmax+xmin)/2, (ymax+ymin)/2-20, 600]
        focal_point =[(xmax+xmin)/2, (ymax+ymin)/2-20, 0]
        title = [0.23, 0.88]
        xlabel = [0.475, 0.11]
        ylabel = [0.125, 0.53]
        cbar = [0.6, 0.1]  # PSF .6 .1
        cbar_length = .135
        cbar_thickness = 60
        FONTSIZE = 50
        spacing_km=100
        round_to=100
        ruler_tick_space=50
        cross_sections = []
        earthquakes = False
        faults=False
        cbar_justification="Right"

    # Open the model volume
    vtk = OpenDataFile(fid)
    vtk = scale_input(vtk, scale_units=preset.scale_units,
                      scale_coords=args.coord_scale)
    # Cook Strait
    # vtk = clip_domain(vtk, xmax=255, ymax=255)
    vtk = clip_domain(vtk, xmin=xmin, ymin=ymin, ymax=ymax, xmax=xmax)
    RenameSource("surface", vtk)

    # Set Camera
    ResetCamera()
    renderView = GetActiveView()
    renderView.InteractionMode = "2D"
    renderView.CameraPosition = camera_position  # [100, 100, 300]
    renderView.CameraFocalPoint = focal_point     # [100, 100, 0]
    renderView.CameraViewUp = [0, 1, 0]
    Render()
    Hide(vtk)

    # Modified coastline
    util_dir = "/Users/Chow/Documents/academic/vuw/forest/utils/vtk_files/scaled/"
    plot_point_cloud(fid=os.path.join(util_dir, "coast.vtk"),
                     reg_name="coast", color=LINECOLOR, point_size=5,
                     xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
    plot_point_cloud(fid=os.path.join(util_dir, "taupo.vtk"),
                     reg_name="taupo", color=LINECOLOR, point_size=5)
    if choice=="tvz":
        plot_landmarks(src=volcanos, glyph_type="Triangle", size=12)
    elif choice == "cook":
        if faults:
            plot_point_cloud(fid=os.path.join(util_dir, "faults.vtk"),
                             reg_name="faults", point_size=2, opacity=1.,
                             xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                             color=LINECOLOR)
        if earthquakes:
            fid= ("/Users/Chow/Documents/academic/vuw/data/events/"
                  "chamberlain_kaikoura/chamberlain_filtered_utm60s_converted.vtk")
            # plot_events(fid, color=rgb_colors("k"), size=2, type="Cross",
            #             linewidth=1, filled=1, representation_type="Surface",
            #             opacity=1.)
            plot_events(fid, color=rgb_colors("w"), size=1.25, 
                        linewidth=2, filled=1,  opacity=1)
            fid= ("/Users/Chow/Documents/academic/vuw/data/events/"
                  "geonet_cook_strait/cook_strait_events_0z_converted.vtk")
            plot_events(fid, color=rgb_colors("pink"), size=2.5, filled=1,
                        position=[0,0,20], linewidth=2)
            # plot_point_cloud(fid, color=rgb_colors("g"), point_size=8,
            #                  xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
    elif choice in ["volumetric", "updates"]:
        plot_events(color=rgb_colors("g"), size=8.)
        plot_stations(size=12.)
        if choice == "updates":
            features()
    elif choice == "seamounts":
        plot_point_cloud(fid=os.path.join(util_dir, "taupo.vtk"),
                         reg_name="taupo", color=LINECOLOR, point_size=5)

    if args.features:
        features()

    # if args.label:
    #     create_text(s=f"[{args.label}]", position=[.24, .82], reg_name="text1",
    #                 fontsize=int(FONTSIZE * 2), color=rgb_colors("k"),
    #                 justification="center")

    # Place perturbation as a 2D glyph
    if args.perturbation and pert_origin is not None:
        plot_point(pert_origin, pert_size)
    
    # Annotate file name and depth in title
    text, _ = create_text(s=f"", position=title,
                          reg_name="text1", fontsize=int(FONTSIZE * 1.25),
                          color=FONTCOLOR, justification="center")

    # Annotate text for X and Y labels
    if not args.xlabel_off:
        create_text(s="X [km]", position=xlabel, reg_name="text1",
            fontsize=int(FONTSIZE), color=FONTCOLOR,
            justification="center")
    if not args.ylabel_off:
        create_text(s="Y [km]", position=ylabel, reg_name="text1",
            fontsize=int(FONTSIZE), color=FONTCOLOR,
            justification="center")


    # Create bounding axes with pre-defined tick marks using rulers
    create_ruler_grid_axes_depth_slice(vtk, tick_spacing_km=ruler_tick_space)
    vsLUT, cbar = set_colormap_colorbar(vtk, position=cbar, length=cbar_length,
                                        orientation="Horizontal",
                                        thickness=cbar_thickness,
                                        text_position="top",
                                        justification=cbar_justification)


    for name in cross_sections:
        origin, normal = DIAGONALS[name.title()]
        plot_cross_section_surface_trace(vtk, normal, origin, color="w",
                                         flip_view=bool(normal[1] > 0),
                                         line_width=5.)


    # Generate slices through the volume at desired depth levels
    slice_ = depth_km
    slice_vtk = depth_slice(vtk, float(slice_))
    tag = f"z_{slice_:0>2}km"
    text.Text = f"{TITLE} [{tag.replace('_', '=')}]"
    if args.label:
        text.Text = f"[{args.label}] {text.Text}"
    rescale_colorscale(vsLUT, src=slice_vtk, vtk=vtk, preset=preset)
    show_colorbar(slice_vtk)

    # Plot the ZM contour
    if args.overlay:
        fid = ("/Users/Chow/Documents/academic/vuw/forest/spread/hessians/"
               "vs_vol_dvs_50ms.vtk")
        contour_overlay(fid=fid, depth=depth_km, xmin=xmin, xmax=xmax, ymin=ymin,
                        ymax=ymax, color="y")

    if args.contour:
        contour_lines(slice_vtk, preset, line_width=4, color=rgb_colors("y"))

    display = GetDisplayProperties(slice_vtk, view=renderView)
    display.SetScalarBarVisibility(renderView, not args.cbar_off)

    set_depth_slice_data_axis_grid(slice_vtk, spacing_km=spacing_km, 
                                   round_to=round_to)

    # Save screenshot, hide slice and move on, job done
    fid_out = os.path.join(save_path, f"{tag}.png")
    SaveScreenshot(fid_out, renderView, ImageResolution=VIEW_SIZE,
                   TransparentBackground=1)
    remove_whitespace(fid_out)
    Hide(slice_vtk, renderView)

    # Delete the min max value points because they'll change w/ each slice
    delete_temp_objects(reg_names=["glyph", "point", "contour"])


def manual_diagonal(fid, preset, save_path=os.getcwd()):
    """
    DDD
    """
    pert_origin = None
    args = parse_args()
    if args.choice is not None:
        choice = args.choice
    else:
        choice = "Mahia"
    if choice is None:
        a=1/0
    elif choice in ["volumetric"]:
        # these are defined in the zero-origin coord system
        # origin, normal = DIAGONALS["Porangahau"]; flip_view=False; choice="Porangahau"  # camera to 500
        # origin, normal = DIAGONALS["Napier"]; flip_view=False; choice="Napier"
        origin, normal = DIAGONALS["Reyners"]; flip_view=True; choice="along_strike"
        tag = choice
        xmin=0
        xmax=462.156
        ymin=0
        ymax=616.130
        zmin = -60
        scale = 2
        dz = 10
        dh = 50
        title = [.8, .7]
        cbar_position = [.8, .51]
        cbar_thickness = 50
        cbar_length = .14
        FONTSIZE = 30
        xsection_label = "A"
        xsection_1 = (.15, .71)
        xsection_2 = (.76, .71)
        interface_pointsize = 2
    elif choice in ["Porangahau", "Mahia", "Seamounts", "Napier", 
                    "Intraplate", "Above", "Mahia_Psf"]:
        origin, normal = DIAGONALS[choice]
        xmin=175
        xmax=462.156
        ymin=100
        ymax=560
        zmin = -30
        scale = 3
        dz = 10
        dh = 50
        cbar_thickness = 50
        cbar_length = .08
        FONTSIZE = 30
        interface_pointsize = 2
        if choice == "Porangahau":
            title = [.76, .72]
            cbar_position = [title[0], .6]
            xsection_label = "B"
            xsection_1 = (.21, title[1])
            xsection_2 = (.735, title[1])

            # pert_origin = [295.54, 251.911, -11 * scale]
            pert_origin = [290, 245, -11 * scale]
            pert_size = 3.5 * scale
        elif choice in ["Mahia", "Mahia_Psf"]:
            title = [.76, .72]
            cbar_position = [title[0], .6]
            xsection_label = "A"
            xsection_1 = (.25, title[1])
            xsection_2 = (.735, title[1])
            # pert_origin = [406.68, 381.1, -11 * scale]
            pert_origin = [405, 380, -11 * scale]
            pert_size = 5 * scale
        elif choice == "Napier":
            title = [.76, .72]
            cbar_position = [title[0], .6]
            xsection_label = "A"
            xsection_1 = (.25, title[1])
            xsection_2 = (.735, title[1])
        elif choice == "Seamounts":
            title = [.83, .72]
            cbar_position = [title[0], .6]
            xsection_label = "C"
            xsection_1 = (.175, title[1])
            xsection_2 = (.8, title[1])
        elif choice in ["Intraplate", "Above"]:
            title = [.83, .72]
            cbar_position = [title[0], .6]
            xsection_label = None # "C"
            xsection_1 = (.175, title[1])
            xsection_2 = (.8, title[1])

            # Perturbations
            if choice == "Intraplate":
                pert_origin = [335, 225, -17.5 * scale]
                pert_size = 7 * scale
                # Above
                pert_origin = [335, 225, -5 * scale]
                pert_size = 5 * scale
            elif choice == "Above":
                pert_origin = [335, 225, -5 * scale]
                pert_size = 5 * scale
        tag = choice
        if choice == "Seamounts":
            flip_view = True
        else:
            flip_view = False  # True means trench parallel
    elif choice in ["flatpoint_as", "wellington_as"]:
        # these are defined in the zero-origin coord system
        if choice == "wellington_as":
            origin, normal = DIAGONALS["Wellington"]
            title = [.71, .7]
            cbar_position = [.71, .58]
            xsection_label = "C"
            xsection_1 = (.3, .71)
            xsection_2 = (.675, .71)
        elif choice == "flatpoint_as":
            origin, normal = DIAGONALS["Flatpoint"]
            title = [.79, .7]
            cbar_position = [title[0], .58]
            xsection_label = "B"
            xsection_1 = (.235, title[1])
            xsection_2 = (.75, title[1])
        tag = choice
        flip_view = False
        xmin=0
        xmax=462.156
        ymin=0
        ymax=616.130
        zmin = -60
        scale = 2
        dz = 10
        dh = 50
        cbar_thickness = 50
        cbar_length = .08
        FONTSIZE = 30
        interface_pointsize = 2
    elif choice == 'a': # "along_strike":
        # these are defined in the zero-origin coord system
        origin, normal = DIAGONALS["Reyners"]
        tag = "along_strike"
        flip_view = True
        xmin=0
        xmax=462.156
        ymin=0
        ymax=616.130
        zmin = -60
        scale = 2
        dz = 10
        dh = 50
        title = [.8, .7]
        cbar_position = [.8, .58]
        cbar_thickness = 50
        cbar_length = .08
        FONTSIZE = 30
        xsection_label = "A"
        xsection_1 = (.15, .71)
        xsection_2 = (.76, .71)
        interface_pointsize = 2
    elif choice in ["cook", "cook_r"]:
        # these are defined in the zero-origin coord system
        if choice == "cook":
            origin, normal = DIAGONALS["Csfaults"]
            xsection_1 = (.28, .76)
            xsection_2 = (.7, xsection_1[1])
            title = [.65, .76]
            cbar_position = [title[0], title[1] - .12]
            flip_view = True
            xsection_label = None # "B"

        elif choice == "cook_r":
            origin, normal = DIAGONALS["Csfaults_R"]
            xsection_1 = (.28, .76)
            xsection_2 = (.7, xsection_1[1])
            title = [.73, .76]
            cbar_position = [title[0], title[1] - .12]
            flip_view = False
            xsection_label = None # "A"

        tag = "cook_strait"
        xmin=0
        xmax=255
        ymin=0
        ymax=205
        zmin=-30
        scale = 3
        dz = 10
        dh = 50
        cbar_thickness =50
        cbar_length=.08
        FONTSIZE = 30
        interface_pointsize = 2

        pert_origin = [145.8, 75.18, -3*scale]
        pert_size = 7 * scale
    elif choice == "tvz":
        # TVZ Rift parallel cross section limited to TVZ only
        # origin = [291.873, 493.837, 0]  # original
        # normal = [-.129802, .08665, 0]
        origin, normal = DIAGONALS["Tvzpsf"]
        # origin, normal = DIAGONALS["Tvz_Rift3"]
        tag = "tvz"
        flip_view = True
        xmin=50 # 50
        xmax=340 # 350
        ymin=325
        ymax=525
        # xmin=185 # 170
        # xmax=340
        # ymin=300 # 250
        # ymax=525 # 600
        zmin=-30
        scale = 2
        dz = 10
        dh = 50
        title = [.765, .86]
        cbar_position=[title[0], .7325]
        cbar_thickness =50
        cbar_length=.1
        FONTSIZE = 30
        xsection_label = None # "A"
        xsection_1 = (.365, title[1])
        xsection_2 = (.82, title[1])
        volcanos = {'Ruapehu': [205.222, 364.035, 0],
                    'Tongariro': [212, 380, 0.],
                    # 'White Island': [344., 560.552, 0],
                    'Tarawera': [280.53000000000003, 484.901, 0],
                    # 'Putauaki': [311.152, 496.12, 0],
                    }
        interface_pointsize = 2
        # TVZ
        pert_origin = [269.99, 489.493, -3 * scale]
        pert_size = 5 * scale
        # TAUPO
        pert_origin = [232.8, 415.7, -4*scale]
        pert_size = 8 * scale
    elif choice == "ruapehu_psf":
        origin, normal = DIAGONALS["Ruapehu_Psf"]
        tag = "ruapehu_psf"
        flip_view = False
        xmin=50
        xmax=350
        ymin=200
        ymax=450
        zmin=-30
        scale = 2
        dz = 10
        dh = 50
        title = [.765, .86]
        cbar_position=[title[0], .7325]
        cbar_thickness =50
        cbar_length=.1
        FONTSIZE = 30
        xsection_label = None # "A"
        xsection_1 = (.365, title[1])
        xsection_2 = (.82, title[1])
        interface_pointsize = 2
        pert_origin = [207.6, 365.88, -12 * scale]
        pert_size = 5 * scale
    elif choice == "tvz_strike":
        # Strike parallel multiple cross sections
        cross_sections = ["Porangahau", "Elsthorpe", "Napier", "Mohaka", "Mahia",
                          "Gisborne", ]
        pick = cross_sections[5]

        origin = DIAGONALS[pick][0]

        normal = [-.45, -.9, 0]
        tag = f"tvz_{pick}"
        flip_view = False
        xmin=50
        xmax=350
        ymin=275
        ymax=605
        zmin=-30
        scale = 3
        dz = 10
        dh = 50
        title = [.84, .825]
        cbar_position=[title[0], .68]
        cbar_thickness =50
        cbar_length=.1
        FONTSIZE = 30
        xsection_label = None
        xsection_1 = None
        xsection_2 = None
        volcanos = {}
        interface_pointsize = 2
        title = None
    elif choice == "tvz_original":
        # these are defined in the zero-origin coord system
        origin = [291.873, 493.837, 0]
        normal = [-.129802, .08665, 0]
        tag = "tvz"
        flip_view = True
        xmin=25
        xmax=390
        ymin=275
        ymax=605
        zmin=-30
        scale = 3
        dz = 5
        dh = 50
        title = [.9, .8]
        cbar_position=[.88, .66]
        cbar_thickness =50
        cbar_length=.1
        FONTSIZE = 30
        xsection_label = "B"
        xsection_1 = (.32, .8)
        xsection_2 = (.87, .8)
        volcanos = {'Ruapehu': [205.222, 364.035, 0],
                    'Tongariro': [212, 380, 0.],
                    'White Island': [344., 560.552, 0],
                    'Putauaki': [311.152, 496.12, 0],
                     }
        interface_pointsize = 2
    elif choice == "francis":
        origin = [361., 417, 2]
        normal = [-.64, -.76, 0]
        tag = "francis_2004_comp"
        flip_view = False
        xmin=300
        xmax=425
        ymin=300
        ymax=445
        zmin=-16
        scale = 4
        dz = 4
        dh = 20
        title = [.59, .565]
        cbar_position=[title[0], .425]
        cbar_thickness =50
        cbar_length=.1
        FONTSIZE = 30
        xsection_label = "C"
        xsection_1 = (.3, title[1])
        xsection_2 = (.57, title[1])
        interface_pointsize =2
    elif choice == "francis_chb":
        origin = [296., 397, 2]
        normal = [-.94, -.33, 0]
        tag = "francis_2004_comp_chb"
        flip_view = False
        xmin = 250
        xmax = 455
        ymin = 300
        ymax = 385
        zmin = -16
        scale = 3
        dz = 4
        dh = 20
        title = [.645, .64]
        cbar_position = [title[0], .505]
        cbar_thickness = 50
        cbar_length = .1
        FONTSIZE = 30
        xsection_label = "D"
        xsection_1 = (.37, title[1])
        xsection_2 = (.625, title[1])
        interface_pointsize = 2

    print(f"{choice}\nORIGIN: {origin}\nNORMAL: {normal}")
    args = parse_args()
    vtk = OpenDataFile(fid)
    vtk = scale_input(vtk, scale_units=preset.scale_units,
                      scale_coords=args.coord_scale)
    vtk = clip_domain(vtk, xmin=xmin, ymin=ymin, ymax=ymax, xmax=xmax,
                      zmin=zmin)
    RenameSource("surface", vtk)

    ResetCamera()
    renderView = GetActiveView()
    Hide(vtk, renderView)

    # Slice across the given cross section plane
    slice_vtk = cross_section(vtk=vtk, normal=normal, origin=origin, name="cs")
    Show(slice_vtk, renderView, "UnstructuredGridRepresentation")

    if args.contour:
        contour_lines(slice_vtk, preset, line_width=4, color=rgb_colors("y"),
                      scale=scale)

    # Reset the colorbounds to data range
    active = GetActiveSource()
    display = GetDisplayProperties(active, view=renderView)
    display.RescaleTransferFunctionToDataRange(False, True)

    if title:
        text, _ = create_text(s=TITLE, position=title,
                              reg_name="text1", fontsize=int(FONTSIZE * 1.25),
                              color=FONTCOLOR, justification="center")

    # A--A'
    if xsection_label:
        text, _ = create_text(s=f"{xsection_label}", position=xsection_1,
                              reg_name="text1", fontsize=int(FONTSIZE * 1.25),
                              color=FONTCOLOR, justification="center")
        text, _ = create_text(s=f"{xsection_label}'", position=xsection_2,
                              reg_name="text1", fontsize=int(FONTSIZE * 1.25),
                              color=FONTCOLOR, justification="center")

    # Use rulers to define the axis grids w/ tickmarks etc.
    ruler_axes_cross_section(src=slice_vtk, normal=normal, scale=scale,
                             dz=dz, dh=dh, flip_view=flip_view)

    scale_vertical_grid(slice_vtk, scale=scale, zgrid_km=[zmin, 0, dz])

    # Generate and rescale the colorbar/ colormap
    vsLUT, cbar = set_colormap_colorbar(vtk, position=cbar_position,
                                        orientation="Vertical",
                                        text_position="left",
                                        thickness=cbar_thickness,
                                        length=cbar_length, fontsize=FONTSIZE)

    rescale_colorscale(vsLUT, src=slice_vtk, vtk=vtk, preset=preset)
    display = GetDisplayProperties(slice_vtk, view=renderView)
    display.SetScalarBarVisibility(renderView, not args.cbar_off)

    if choice == "tvz":
        for origin_ in volcanos.values():
            mark_surface_point(origin_, normal, size=10., z_value=10)
        fid = ("/Users/Chow/Documents/academic/vuw/"
               "forest/utils/vtk_files/scaled/taupo.vtk")
        plot_2D_surface_projection(fid, origin, normal, flip_view=flip_view,
                                   scale=scale, offset=1,
                                   pointsize=4.5, color=rgb_colors("k"),
                                   position=[0., 0., 0], xmin=xmin, xmax=xmax,
                                   ymin=ymin, ymax=ymax, zmin=zmin,
                                   color_by=False)

    # Plot perturbation marker for point spread functions
    if args.perturbation and pert_origin is not None:
       angle = normal_to_angle(origin[0] - 5, origin[1]-5) 
       if not flip_view:
           angle *= -1
       plot_point(pert_origin, pert_size, [90, 0, angle])

       # angle = normal_to_angle(origin[1], origin[0]) 
       # plot_point(pert_origin, pert_size, [90, 0, angle+20])

    # Plot the ZM contour
    if args.overlay:
        fid = ("/Users/Chow/Documents/academic/vuw/forest/spread/hessians/"
               "vs_vol_dvs_50ms.vtk")
        contour_overlay(fid=fid, origin=origin, normal=normal, scale=scale,
                        xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin)

    # Plot the interface model of Williams et al. (2013)
    fid = ("/Users/Chow/Documents/academic/vuw/"
           "forest/utils/vtk_files/scaled/interface.vtk")
    plot_2D_surface_projection(fid, origin, normal, flip_view=flip_view,
                               scale=scale, offset=1,
                               xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                               zmin=zmin, pointsize=interface_pointsize,
                               color_by=False)
    fid = ("/Users/Chow/Documents/academic/vuw/"
           "forest/utils/vtk_files/scaled/coast.vtk")
    plot_2D_surface_projection(fid, origin, normal, flip_view=flip_view,
                               scale=scale, offset=1,
                               pointsize=4.5, color=rgb_colors("k"),
                               position=[0., 0., 0], xmin=xmin, xmax=xmax,
                               ymin=ymin, ymax=ymax, zmin=zmin,
                               color_by=False)


    # Reset camera view to be normal to the plane. Specific to this plane
    set_camera(choice)


    fid_out = os.path.join(save_path, f"{tag}.png")
    SaveScreenshot(fid_out, renderView, ImageResolution=VIEW_SIZE,
                   TransparentBackground=1)
    remove_whitespace(fid_out)


def parse_args():
    """
    Allow user defined arguments to fine tune figures
    """
    # Command line arguments to define behavior of script
    parser = argparse.ArgumentParser()

    # MAIN ARGUMENTS
    parser.add_argument("files", nargs="+", help="file to plot with paraview")
    parser.add_argument("-v", "--verbose", action="store_true", default=False,
                        help="print output messages during the plotting")
    parser.add_argument("-o", "--output", type=str, help="output path",
                        default=os.path.join(os.getcwd(), "output"))

    # COLORMAP
    parser.add_argument("-p", "--preset", type=str,
                        help="preset colormap and labels given file type")
    parser.add_argument("--nvalues", type=int, help="number of cbar values")
    parser.add_argument("-b", "--bounds", type=str,
                        help="Manually set the bounds of the colorbar, "
                             "overriding the default or preset bound values",
                        default=None)
    parser.add_argument("--cbar_orientation", type=str, default=None,
                        help="Select orientation of colorbar, 'horizontal' or "
                             "'vertical'")

    # GLOBAL DEFAULTS
    parser.add_argument("--fontcolor", type=str, default="k",
                        help="The default color to use for any text on plot")
    parser.add_argument("--linecolor", type=str, default="k",
                        help="The default color to use for any plot attributes")
    parser.add_argument("-F", "--fontsize", type=int, default=50,
                        help="The default color fontsize to use for any text.")
    parser.add_argument("-V", "--viewsize", type=list, default=[2100, 2100],
                        help="The default viewing size for all screenshots "
                             "made; controls the resolution as well as the "
                             "relative sizing of objects in the screenshot")

    # SLICES, CROSS SECTIONS, PROJECTIONS
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
    parser.add_argument("-d", "--diagonal", action="store_true",
                        help="plot diagonal cross-section with normal and "
                             "origin arguments",
                        default=False)
    parser.add_argument("-D", "--diagonals", nargs="+",
                        help="choose which diagonal cross sections, must"
                             "match the names defined in DIAGONALS")
    parser.add_argument("-i", "--interface", action="store_true",
                        help="plot projection of model onto plate interface "
                             "of Charles Williams",
                        default=False)
    parser.add_argument("--taupo", action="store_true",
                        help="plot trench parallel cross section that makes a "
                             "line between Ruapehu and White Island",
                        default=False)

    parser.add_argument("-A", "--all", action="store_true", default=False,
                        help="shorthand to make all default slices, "
                             "same as '-xyzit'")

    # PLOT EXTRAS
    parser.add_argument("--coord_scale", default=1E-3, type=float, 
                        help="Scale units of coordinate system, e.g. m -> km")
    parser.add_argument("-t", "--title", default=None,
                        help="Title to annotate somewhere on the figure, "
                             "defaults to the name of the file")
    parser.add_argument("-c", "--contour", action="store_true",
                        help="generate contour lines ontop of the slices",
                        default=False)
    parser.add_argument("--overlay", action="store_true",
                        help="zeroth moment contour overlay on manual slices",
                        default=False)
    parser.add_argument("-l", "--outline", action="store_true",
                        help="plot the coastline and shoreline of lake Taupo",
                        default=False)
    parser.add_argument("-s", "--sources", action="store_true",
                        help="plot events as glyphs", default=False)
    parser.add_argument("-r", "--receivers", action="store_true",
                        help="plot stations as glyphs", default=False)
    parser.add_argument("-S", "--srvtk", type=str, default=None,
                        help="plot an sr.vtk file for sensitivity kernels",)
    parser.add_argument("-f", "--faults", action="store_true",
                        help="plot north island active faults on surface proj "
                        "only", default=False)
    parser.add_argument("-e", "--slowslip", action="store_true",
                        help="plot slow slip slip event slip patches",
                        default=False)
    parser.add_argument("-P", "--perturbation", action="store_true", 
                        help="plot perturbations",
                        default=False)
    parser.add_argument("-w", "--williams", action="store_true", default=False,
                        help="Plot the Hikurangi interface model of Williams"
                             "et al. (2013). If depth slice, will plot contour"
                             "lines pre defined by the .vtk file. If xsection"
                             "will plot a line tracing the interface.")
    parser.add_argument("--mark", default=False,
                        help="Pick an origin point to mark on depth slices")
    parser.add_argument("--landmarks", action="store_true", default=False,
                        help="Plot landmarks as 2D glyphs for easier reference "
                             "on depth slices")
    parser.add_argument("--surface_traces", nargs="+", default=None,
                         help="Plot surface traces on depth slice plots for "
                              "given named origin/normal combinations for"
                              "reference")
    parser.add_argument("--checkers", type=str, default=False,
                        help="Plot checker overlay for point spread test")
    parser.add_argument("--minmax", action="store_true", default=False,
                        help="Mark the  min. and max. values on the slice. "
                             "Wont work with surface projections")
    parser.add_argument("--cbar_off", action="store_true", default=False,
                        help="Colorbar off, it is on by default")
    parser.add_argument("--anno_off", action="store_true", default=False,
                        help="Cross section annotations off, on by default")
    parser.add_argument("--xlabel_off", action="store_true", default=False,
                        help="Turn off xlabel")
    parser.add_argument("--ylabel_off", action="store_true", default=False,
                        help="Turn off ylabel")

    # SPECIFIC PUBLICATION ADDITIONS
    parser.add_argument("--choice", type=str, default=None)
    parser.add_argument("--features", action="store_true", default=False,
                        help="text labels for features A-E in paper")
    parser.add_argument("--label", type=str, default=None,
                        help="labels for multi-panel figures")
    parser.add_argument("--xlabel", type=int, default=1,
                        help="xlabel for depth slices, 0 for off, 1 for on")
    parser.add_argument("--ylabel", type=int, default=1,
                        help="ylabels for depth slices, 0 for off, 1 for on")

    parser.add_argument("--manualz", action="store_true", default=False,
                        help="Make hardcoded figures")
    parser.add_argument("--manuald", action="store_true", default=False,
                        help="Make hardcoded figures")


    # CROSS SECTION FINE TUNE
    parser.add_argument("--depth_km", type=float, default=60,
                        help="For any vertical cross sections (Y-axis figure "
                             "normal to Z axis of volume), define the depth"
                             "cutoff of the screenshot as usually were not "
                             "interested in looking at the entire volume. "
                             "Units of km, positive values only.")
    parser.add_argument("--zscale", type=int, default=2,
                        help="Set the scale for the Z-axis on cross sections")
    parser.add_argument("--dz", type=int, default=10,
                        help="Vertical grid spacing on cross sections in km")
    parser.add_argument("--anno", type=float, default=0.35,
                        help="Height of annotation and colorbar for cross "
                             "sections. Will change depending on the scale and"
                             "cutoff depth used.")
    parser.add_argument("--origin", nargs="+", default=None,
                        help="Origin for diagonal cross sections")
    parser.add_argument("--normal", nargs="+", default=None,
                        help="Normal vector for diagonal cross sections")

    return parser.parse_args()



if __name__ == "__main__":
    args = parse_args()

    # Set some global variables
    FONTCOLOR = rgb_colors(args.fontcolor)
    LINECOLOR = rgb_colors(args.linecolor)
    FONTSIZE = args.fontsize
    VIEW_SIZE = args.viewsize

    # Set up the active render view
    renderView = GetActiveViewOrCreate("RenderView")
    renderView.ViewSize = VIEW_SIZE
    renderView.InteractionMode = "2D"
    SetActiveView(renderView)
    Render()

    # Convert coordinates
    if args.coord_scale != 1:
        for key, val in DIAGONALS.items():
            origin, norm = val
            DIAGONALS[key] = (convert_coords(origin), norm)

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

        if args.title is None:
            TITLE = save_fid
        else:
            TITLE = args.title

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
        if args.bounds is not None:
            preset.bounds = parse_bounds(args.bounds)
        if args.nvalues is not None:
            preset.nvalues = args.nvalues

        # ======================================================================
        # PUBLICATION FIGURES HACKED IN
        # ======================================================================
        if args.manuald:
            manual_diagonal(data_fid, preset=preset, save_path=save_path)
        if args.manualz:
            manual_depth_slice(data_fid, preset=preset, save_path=save_path)


        # ======================================================================
        # DEPTH SLICES
        # ======================================================================
        zslices = []
        if args.default_zslices or args.all:
            zslices += ["surface", "2-20,2", "25-50,5", "60-100,20"]
        if args.zslices:
            zslices += args.zslices

        if zslices:
            make_preplot(args)
            zslices = parse_slice_list(zslices)
            if args.verbose:
                print(f"\tGenerating Z slices for {zslices}")
            depth_slices(data_fid, zslices, preset=preset, save_path=save_path)
            reset()

        # ======================================================================
        # X-NORMAL CROSS SECTIONS
        # ======================================================================
        if args.xslices or args.all:
            # Do not go for 0 or 100 % because you might end up off the model
            xslices = list(range(5, 100, 10))
            if args.verbose:
                print(f"\tGenerating X-normal slices for {xslices}")
            cross_sections(data_fid, xslices, normal="x", preset=preset,
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
            cross_sections(data_fid, yslices, normal="y", preset=preset,
                                save_path=save_path)
            reset()

        # ======================================================================
        # DIAGONAL CROSS SECTIONS
        # ======================================================================
        if args.diagonal:
            diagonal(data_fid, preset, save_path=save_path)

        # Parse any preset diagonal cross sections
        diagonal_names = []
        if args.diagonals or args.all:
            diagonal_names += DIAGONALS.keys()

        if args.diagonals:
            if args.verbose:
                print("\tGenerating diagonal cross sections")
            for name in args.diagonals:
                diagonal(data_fid, preset, name=name, save_path=save_path)


        # ======================================================================
        # PROJECTION ONTO INTERFACE
        # ======================================================================
        if args.interface or args.all:
            if args.verbose:
                print("\tGenerating interface projection")
            make_preplot(args)
            interface(data_fid, preset, save_path=save_path)
            reset()
           

