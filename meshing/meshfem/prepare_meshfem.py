#!/usr/bin/env python3
"""
Generate required files to run Specfem3D Cartesians internal mesher, Meshfem3D
while keeping track of the internal requirements of Meshfem3D and generally 
simplifying the meshing procedure.

LIMITATIONS:
    * This script is geared towards creating a very simple cubed-element mesh 
    with simple mesh refinement/coarsening layers.
    * Currently hardcoded for two doubling layers, both in the templates and in 
    the functions below. If more or less layers are required, one would need to 
    edit the templates to add more layers and interfaces, and then change the 
    string formatting in the write functions below

REQUIRES:
    Python>=3.6
    numpy
    template_Mesh_Par_file 
    template_interfaces.dat

MESHING NOTES:
    1) Doubling layers control the horizontal doubling of element area. Only x2
    2) Interfaces control vertical size of elements. Can be arbitrary scaling.
    3) Doubling layers placed above interfaces causes high skewness, thats bad!
    4) Put interfaces 1-2 layers above doubling layers if your elements are 
       square at the top, to reduce skewness.

USAGE:
    $ python prepare_meshfem.py  # looks for file parmesh.json

    OR

    $ python prepare_meshfem.py my_custom_parmesh.json  # must be in same dir.

|==============================================================================|
|                       JSON PARAMETER FILE PARAMETERS                         |
|                                                                              |
| The following parameters must be set in 'parmesh.json' (or user defined file)|
| NOTE: To set values as None in JSON, set them as: null                       |
|==============================================================================|

:type tag: str
:param tag: the name of your mesh, this tag is used to save the log file
:type dir_name: str
:param tag: directory to save all the output files to, this will end up being
    your meshfem3D_files/ directory in SPECFEM3D
:type lat_min: float
:param lat_min: minimum latitude in degrees OR minimum Y value in km 
:type lat_max: float
:param lat_max: maximum latitude in degrees OR maximum Y value in km
:type lon_min: float
:param lon_min: minimum longitude in degrees OR minimum X value in km
:type lon_max: float
:param lon_max: maximum longitude in degrees OR maximum X value in km
:type utm_projection: int
:param utm_prjection: UTM projection to convert lat/lon. 
    OPTIONAL: Not needed if `suppress_utm_proj==True`
    NOTE: requires `pyproj` to convert from lat/lon to UTM
:type suppress_utm_proj: bool
:param suppress_utm_proj: if True, will not try to convert from lat/lon to UTM,
    and coordinates can be in lat/lon or in cartesian XYZ.
    if False, lat_min etc. expected to be coordinates in units of degrees
:type mesh_depth_km: float
:param mesh_depth_km: mesh depth in units of kilometers
:type mesh_top_km: float
:param mesh_top_km: top of topography in units of kilometers. 
    If no topography then set to 0.
:type interfaces: list of float
:param interfaces: the depth locations of interfaces where doubling occurs
    in units of kilometers. 
    NOTE: Atleast one interface must be included otherwise no interfaces.dat
          file will be written out.
:type interface_fids: list of str
:param interface_fids: name of the output interface files corresponding to
    the list of interfaces.
    NOTE: Must match len(interfaces)
:type interface_increase: list of float OR float
:param interface_increase: An integer multiple of the previous vertical element
    size starting from the top of the mesh. For example if your mesh has 1 
    interface with dz_top=1km, an interface increase of 2 would make dz_bot=2km.

    NOTE: must match len(interfaces)
    NOTE: if given as a single float, each interface will be multipled by
          the float. e.g., if 2, each interface will cause vertical element
          size to double

    EXAMPLE:
    ...
    interfaces = [some_value]
    interface_fids = [interface]
    interface_increase = [2]
    ...

                LAYERS  DZ
    TOP OF MESH  1      1
                 2      1
    INTERFACE    3      2
                 4      2
:type regions: list of float
:param regions: depths at which to create regions, or a collection of layers. 
    These regions are related to the materials in the mesh (if using the simple 
    internal mesher method for assigning material properties). 
    Or they are used to assign external .xyz tomography models to the mesh. 
:type materials: list of list of float OR empty list
:param materials: materials to assign to each region. If empty, some random
    default values will be assigned. This is only important if you are setting
    `model=default` in the SPECFEM3D Par_file.

    NOTE each list must be of length 5, corresponding to the following values:
        [rho, vp, vs, qkappa, qmu]
    NOTE: len(materials) = len(regions) + 1 
    EXAMPLE: if e.g., regions = [20, 30] then
        materials = [[a, b, c, d, e], [f, g, h, i, j], [k, l, m, n, o, p]]
:type nproc: int
:param nproc: number of processors to split the mesh X and Y into. Usually good
    idea to make this an integer multiple of the number of cores per node
    on your cluster. 
:type shortest_period_s: float
:param shortest_period_s: the shortest period that the mesh is expected to
    resolve. Not a strict requirement, but gives an idea of the necessary
    element sizes required to resolve given features. A minimum
    resolvable period test must be done to get the actual resolution of the mesh.
:type top_vertical_element_ratio: float or None
:param top_vertical_element_ratio: Scale the vertical element size by the
    horizontal element size. Set to 1 if you want roughly cube-like elements
    at the top of your mesh. i.e.,
    
    dz = (dx + dy) / 2 * top_vertical_element_ratio

    NOTE: if None, vertical element size will be set as the 'v' in 
          grid_space_top_hv OR whatever is calculated in OPTION 2 (see below).

ELEMENT SIZE OPTION 1 - Manually set element sizes at top of mesh:

:type grid_space_top_hv: list of float or None
:param grid_space_top_hv: desired horizontal and vertical grid spacing at the
    top of the mesh [h, v] in units of km. 
    NOTE: If None, will use option 2
    NOTE: This number is NOT enforced, it is simply a recommendation or starting
          point for the script. You may have to play around with this parameter,
          plus nproc, plus the mesh dimensions to get elements of the desired
          size. The script will do its best to get close to the values of h, v.

ELEMENT SIZE OPTION 2 - Dynamically determine element size at top of mesh based
    on shortest wavelength and desired points per wavelength:

:type vs_min_km_per_s: float
:param vs_min_km_per_s: smallest expected wavespeed in the model, in km/s
    as with `shortest_period_s`, gives an idea of the necessary element sizes
    but is not a hard requirement. An MRP test will need to be performed
    NOTE: if grid_space_top_hv given, will not be used
:type points_per_wl: int
:param points_per_wl: desired points per wavelength for numerical accuracy.
    bare minimum = 2, SPECFEM recommended = 5

|==============================================================================|
|                       JSON PARAMETER FILE PARAMETERS                         |
|==============================================================================|
"""
import os
import sys
import json
import shutil
import logging
import traceback
import numpy as np

# Set up the logger with very simple format
logging.basicConfig(level=logging.DEBUG,
                    format="%(message)s",
                    filename=f"log_mesh_temp.out",
                    filemode="w")

# Write logging to console
console = logging.StreamHandler()
console.setLevel(logging.INFO)
logging.getLogger('').addHandler(console)
logger = logging.getLogger("mesher")


def myround(x, base=5, choice='near', rtype=int):                                           
    """                                                                          
    Round value x to nearest base, round 'up','down' or to 'near'est base        
                                                                                 
    :type x: float                                                               
    :param x: value to be rounded                                                
    :type base: int                                                              
    :param base: nearest integer to be rounded to                                
    :type choice: str                                                            
    :param choice: method of rounding, 'up', 'down' or 'near'                    
    :rtype roundout: float
    :return: rounded value                                                       
    """                                                                          
    if choice == 'near':                                                         
        roundout = rtype(base * round(float(x) / base))                            
    elif choice == 'down':                                                       
        roundout = rtype(base * np.floor(float(x) / base))                         
    elif choice == 'up':                                                         
        roundout = rtype(base * np.ceil(float(x) / base))                          
                                                                                 
    return roundout  


def lonlat_utm(lon_or_x, lat_or_y, utm_zone=-60, inverse=False):                 
    """                                                                          
    convert latitude and longitude coordinates to UTM projection                 
    From Pyatoa                                                                  
                                                                                 
    :type lon_or_x: float or int                                                 
    :param lon_or_x: longitude value in WGS84 or X in UTM-'zone' projection      
    :type lat_or_y: float or int                                                 
    :param lat_or_y: latude value in WGS84 or Y in UTM-'zone' projection         
    :type utm_zone: int                                                          
    :param utm_zone: UTM zone for conversion from WGS84                          
    :type inverse: bool                                                          
    :param inverse: if inverse == False, latlon => UTM, vice versa.              
    :rtype x_or_lon: float                                                       
    :return x_or_lon: x coordinate in UTM or longitude in WGS84                  
    :rtype y_or_lat: float                                                       
    :return y_or_lat: y coordinate in UTM or latitude in WGS84                   
    """                                                                          
    try:
        from pyproj import Proj
    except ModuleNotFoundError:
        sys.exit("To use UTM projections, you must install 'pyproj'")
                                                                                 
    # Determine if the projection is north or south                              
    if utm_zone < 0:                                                             
        direction = "south"                                                      
    else:                                                                        
        direction = "north"                                                      
    # Proj doesn't accept negative zones                                         
    utm_zone = abs(utm_zone)                                                     
                                                                                 
    # !!! Outdated method of calling Proj
    # Proj requires a string to tell it how to convert the coordinates           
    # projstr = (f"+proj=utm +zone={utm_zone}, +{direction} +ellps=WGS84"          
    #            " +datum=WGS84 +units=m +no_defs")                                
    # my_proj = Proj(projstr)                           

    # Initiate a Proj object and convert the coordinates                         
    my_proj = Proj(proj="utm", zone=abs(utm_zone), south=True, ellps="WGS84",
                   datum="WGS84", units="m", no_defs=True)
    x_or_lon, y_or_lat = my_proj(lon_or_x, lat_or_y, inverse=inverse)            
                                                                                 
    return x_or_lon, y_or_lat            


def set_parameters(fid="./parmesh.json"):
    """
    Define the necessary parameters here, these will be accessed by the
    constraint fuctions throughout the script. Defaults for a New Zealand North
    coarse mesh are set here.

    :rtype: dict
    :return: dictionary of parameters
    """
    try:
        with open(fid, "r") as f:
            parameters = json.load(f)
    except Exception as e:
        print("\n\nERROR READING JSON FILE") 
        traceback.print_exc()
        print("ERROR READING JSON FILE\n\n")
        sys.exit(-1)

    # Print the parameters so the User knows what they've chosen
    logger.info(f"PREPARING MESH FOR MESHFEM3D")
    logger.info(f"\n{parameters}\n")
    logger.info(f"PARAMETERS FOR: {parameters['tag']}")
    for key, item in parameters.items():
        logger.info(f"\t{key}: {item}")

    # Set some internal parameters 
    parameters["ndoublings"] = len(parameters["doubling_layers"])

    # Make sure num materials matches num regions, if materials given
    if parameters["materials"]:
        assert(len(parameters["materials"]) == len(parameters["regions"]) + 1
               ), ("number of regions must match number of materials, or else "
                   "materials must be an empty list")
    else:
        # Set None if empty to use default values
        parameters["materials"] = None

    # Make sure there are enough fids for each interface, or that none specified
    assert((len(parameters["interface_fids"]) == 0) or
           (len(parameters["interfaces"]) == len(parameters["interface_fids"]))
           ), "interfaces number must match interface_fids number"

    # Convert latlon to UTM
    if not parameters["suppress_utm_proj"]:
        x_min, y_min = lonlat_utm(parameters["lon_min"], parameters["lat_min"],
                                  parameters["utm_projection"])
        x_max, y_max = lonlat_utm(parameters["lon_max"], parameters["lat_max"],
                                  parameters["utm_projection"])
    else:
        x_min = parameters["lon_min"]
        x_max = parameters["lon_max"]
        y_min = parameters["lat_min"]
        y_max = parameters["lat_max"]

    # Set the UTM coordinates in the parameters
    parameters["x_min"] = x_min
    parameters["y_min"] = y_min
    parameters["x_max"] = x_max
    parameters["y_max"] = y_max
    parameters["x_length_km"] = (x_max - x_min) * 1E-3
    parameters["y_length_km"] = (y_max - y_min) * 1E-3

    return parameters


def minimum_grid_spacing(slowest_wavespeed, shortest_period, points_per_wl=2,
                         top_vertical_element_ratio=1,
                         grid_space_top_hv=None):
    """
    Define minimum grid spacing based on slowest wavespeed and shortest period

    :type slowest_wavespeed: float
    :param slowest_wavespeed: slowest wavespeed expected in model
    :type shortest_period: float
    :param shortest_period: shortest period to be resolved in mesh
    :rtype: float
    :return: minimum element spacing of the mesh
    """
    logger.info("CALCULATING MINIMUM GRID SPACE AT TOP OF MESH")
    if grid_space_top_hv:
        grid_space_h, grid_space_v = grid_space_top_hv
        logger.info("\tGrid space XYZ given")
        logger.info(f"\tMIN GRID SPACE HORIZONTAL: {grid_space_h}")
        logger.info(f"\tMIN GRID SPACE VERTICAL: {grid_space_v}")
        return grid_space_h, grid_space_v

    else:
        logger.info("\tManually calculating grid space based on T and V")
        logger.info(f"\tdesired grid points per wavelength: {points_per_wl}")
        min_grid_space = (shortest_period * slowest_wavespeed) / points_per_wl
        logger.info(f"\tT_min * V_min / points_per_wl = {min_grid_space} km")

        # I don't think we need to round this value? We can have non-integer
        # values of grid space 
        # grid_space_h = myround(min_grid_space, 1, 'down') or min_grid_space

        # Scale horizontal grid space to get vertical grid space
        grid_space_h = min_grid_space
        if top_vertical_element_ratio is not None:
            grid_space_v = grid_space_h / top_vertical_element_ratio
        else:
            grid_space_v = grid_space_h

        # Same as above, we can have non-integer values of grid space v
        # if grid_space_h == myround(min_grid_space, 1, 'down'):
        #     logger.info(f"\trounded to nearest factor of 2")

        logger.info(f"\tMIN GRID SPACE HORIZONTAL = {grid_space_h} km")
        logger.info(f"\tMIN GRID SPACE VERTICAL = {grid_space_v} km")
        
        return grid_space_h, grid_space_v


def cfl_condition(dx, vmax, C=1):
    """
    Courant-Friedrichs-Lewy Condition cointrolling the minimum time step given
    smallest velocity and grid spacing present in the numerical approximation
    schema. C gives the dimensionless Courant number which is defaulted to 1
    for an explicit (time-marching solver)
    
    :type dx: float
    :param dx: smallest element spacing in mesh (can be dx, dy or dz) in km
    :type vmax: float
    :param vmax: largest velocity in mesh, in km/s
    :type C: float
    :param C: Courant number, default 1
    """
    logger.info("COURANT CRITERION")
    logger.info(f"\tdt <= {dx*C/vmax:.2E}, for C={C}")


def number_of_processors(nproc, x_length, y_length):
    """
    Set the number of processors base on the ratio of mesh length

    :type nproc: int
    :param nproc: total number of MPI processors to be used
    :type x_length: float
    :param x_length: length of the mesh in x-direction (e.g. km)
    :type y_length: float
    :param y_length: length of the mesh in the y-direction
    :rtype nproc_x: int
    :return nproc_x: number of processors in the x-direction
    :rtype nproc_y: int
    :return nproc_y: number of the processors in the y-direction
    """
    logger.info("CALCULATING NUMBER OF PROCESSORS")
    # Determine the ratio
    ratio = min(x_length, y_length) / max(x_length, y_length)

    # Define which direction is shorter
    if x_length < y_length:
        short_direction = "x"
    else:
        short_direction = "y"

    logger.info(f"\tx-axis length = {x_length}")
    logger.info(f"\ty-axis length = {y_length}")
    logger.info(f"\tshort direction is '{short_direction}'")

    # Start guessing processor ratios at the square root, until 2 integers found
    guess_a = round(np.sqrt(ratio * nproc))
    while True:
        guess_b = nproc / guess_a
        if guess_b.is_integer():
            # Assign the short direction correctly
            if short_direction == "x":
                nproc_x = min(guess_a, guess_b)
                nproc_y = max(guess_a, guess_b)
            else:
                nproc_x = max(guess_a, guess_b)
                nproc_y = min(guess_a, guess_b)
            logger.info(f"\tNPROC_X = {int(nproc_x)}\n"
                        f"\tNPROC_Y = {int(nproc_y)}"
                        )
            return int(nproc_x), int(nproc_y)
        else:
            guess_a += 1


def adjust_vertical_element_size(dx, dy, top_vertical_element_ratio):
    """
    DX and DY are dynamically determined by the desired number of elements,
    mesh domain size, etc. Meanwhile DZ can be set as your heart desires.
    This function simply scales DZ to DX and DY to get cube-like elements, or
    other ratios if desired

    Simply returns: (dx + dy) / 2 * top_vertical_element_ratio

    :type dx: float
    :param dx: element size in x direction at top of mesh
    :type dy: float
    :param dy: element size in x direction at top of mesh
    :type top_vertical_element_ratio: float
    :param top_vertical_element_ratio: amount to scale DX by to get DZ
    :rtype: float
    :return: DZ as a scaled function of DX and DY
    """
    logger.info("SCALING VERTICAL ELEMENT SIZE BASED ON DX and DY")
    dz = (dx + dy) / 2 * top_vertical_element_ratio
    dz = myround(dz, 0.25, "near", rtype=float)  # round off to get clean values
    logger.info(f"\tDZ = (dx + dy)/2 * {top_vertical_element_ratio} ~= {dz}")
    return dz


def number_of_elements(nproc_x, nproc_y, x_length, y_length, grid_space,
                       shortest_period_s, ndoublings):
    """
    Define the number of elements in each horizontal direction at the top of the
    mesh based on number of processors

    :type nproc_x: int
    :param nproc_x: number of processors in the x-direction
    :type nproc_y: int
    :param nproc_y: number of processors in the y-direction
    :type x_length: float
    :param x_length: length of the mesh in the x-direction
    :type y_length: float
    :param y_length: length of the mesh in the y-direction
    :type grid_space: float
    :param grid_space:
    :type shortest_period_s: float
    :param shortest_period_s: shortest period to be resolved in mesh
    :rtype nex_x: number of elements in the x-direction
    :return nex_x: int
    :rtype nex_y: number of elements in the y-direction
    :return nex_y: int
    """
    logger.info("CALCULATING NUMBER OF ELEMENTS AT TOP OF MESH")

    # Based upon benchmarks against semi-analytical discrete wavenumber 
    # synthetic seismograms [Komatitsch et al., 2004], determined that a 
    # NEX_XI = 288 run is accurate to a shortest period of roughly 2 s. 
    # Therefore, since accuracy is determined by the number of grid points per 
    # shortest wavelength, for any particular value of NEX_XI the simulation 
    # will be accurate to a shortest period determined by 
    c = 1
    nex_x = 0
    nex_y = 0
    while ((nex_x < 2 * 288 / shortest_period_s) or \
            (nex_y < 2 * 288 / shortest_period_s)):

        # Meshfem requires the number of grid points be an integer multiple of 
        # 8 times the number of processors in a given direction
        nex_x = myround(x_length / grid_space, c * nproc_x * 8, "up")
        nex_y = myround(y_length / grid_space, c * nproc_y * 8, "up")
        
        # Increment by integers until we satisfy this requirement
        c += 1
        if c >= 50:
            raise OverflowError("nex_x seems too large for given nproc_x")

    # Ensure that the short direction is maintained
    while (nproc_x < nproc_y) and (nex_x > nex_y):
        nex_y += nproc_y * 8

    logger.info(f"\tinteger multiple found as c={c}")
    logger.info(f"\tNEX_X = {nex_x}\n"
                f"\tNEX_Y = {nex_y}")

    # Re-check assertions mandated by Meshfem 
    assert(nex_x >= 2 * 288 / shortest_period_s)
    assert(nex_y >= 2 * 288 / shortest_period_s)
    
    # NEX_PER_PROC_XI must be a multiple of 2 * 2**NDOUBLINGS
    if not ((nex_x / nproc_x)  %  2 * 2 ** ndoublings):
        logger.info("\t!!! enforcing NEX_PER_PROC_XI == 2c * 2 ** ndoublings")
        nex_x_new = 0
        c = 1
        while nex_x_new < nex_x:
            nex_x_new = c * nproc_x * (2 * 2 ** ndoublings)
            c += 1
        logger.info(f"\t\tNEX_X: {nex_x} -> {nex_x_new}")
        nex_x = nex_x_new
    if not ((nex_y / nproc_y)  %  2 * 2 ** ndoublings):
        logger.info("\t!!! enforcing NEX_PER_PROC_ETA == 2c * 2 ** ndoublings")
        nex_y_new = 0
        c = 1
        while nex_y_new < nex_y:
            nex_y_new = c * nproc_y * (2 * 2 ** ndoublings)
            c += 1
        logger.info(f"\t\tNEX_Y: {nex_y} -> {nex_y_new}")
        nex_y = nex_y_new

    # Minimum element number set by Eq. 3.1 in Specfem manual
    dx = x_length / nex_x
    dy = y_length / nex_y

    logger.info(f"\tDX = {dx:.4f}km\n"
                f"\tDY = {dy:.4f}km")

    logger.info("\tNote: if values below too large, adjust 'shortest_period_s'")
    logger.info("\tNote: if values below too small, check  'suppress_utm_proj'")
    logger.info(f"\tGRID_SPACE_H - DX = {abs(dx-grid_space)}")
    logger.info(f"\tGRID_SPACE_H - DY = {abs(dy-grid_space)}")

    return int(nex_x), int(nex_y), dx, dy


def calculate_nelements(nex_x, nex_y, layers, doubling_layers):
    """
    Calculates the total number of elements in the mesh, taking horizontal
    doubling layers into consideration
    
    :type nex_x: int
    :param nex_x: number of elements in the x direction at the surface
    :type nex_y: int
    :param nex_y: number of elements in the y direction at the surface
    :type layers: list of int
    :param layers: the number of elements in each doubling layer
    :rtype int:
    :return: the total number of elements in the mesh
    """
    logger.info("APPROXIMATING TOTAL NUMBER OF ELEMENTS IN MESH")
    nelements = 0
    # Figure out the grid spacing at the bottom of the mesh
    nelem_x = nex_x / (2 ** len(doubling_layers))
    nelem_y = nex_y / (2 ** len(doubling_layers))
    # Start counting from the bottom of the mesh
    for i in range(1, sum(layers) + 1):
        if i in doubling_layers:
            nelem_x *= 2
            nelem_y *= 2
        nelements += nelem_x * nelem_y

    assert(nelem_x == nex_x)
    assert(nelem_y == nex_y)
    nelements = int(nelements)
    logger.info(f"\tAPPRX_NUM_ELEM = {nelements}")
    
    return nelements


def interface_layers(top_of_mesh, depth, interfaces, grid_space_z, 
                     interface_increase=2):
    """
    Determine the number of elements in the vertical direction (depth) based on
    the desired grid spacing at the top of the mesh, the desired depth of the
    mesh and the desired number of interfaces which control vertical element
    size.

    :type top_of_mesh: float
    :param top_of_mesh: top of the mesh in km (negative for topography)
    :type depth: float
    :param depth: depth of the mesh in km (positive)
    :type interfaces: list of float
    :param interfaces: interfaces for vertical element size manipulations
    :type grid_space_z: float
    :param grid_space_z: desired vertical grid spacing at the top of the mesh
    :type interface_increase: int
    :param interface_increase: the amount to multiply vertical element spacing 
        at each interface
    :rtype layers: list of float
    :return: number of elements in each layer corresponding to each desired
        interface. layers start counting from the bottom of the mesh
    """
    logger.info("CALCULATING NUMBER OF VERTICAL LAYERS (interfaces)")

    # Start interfaces from the top, include the bottom interface of depth
    # all_interfaces = [top_of_mesh] + interfaces + [depth]
    all_interfaces =  interfaces + [depth]
    all_interfaces.sort()

    layers = []
    # Allows for both list-like interface increases and single valued
    if isinstance(interface_increase, int):
        interface_increase = [interface_increase] * len(all_interfaces)

    # Set up how the vertical element will change at each interface    
    # First entry needs to be 1 because no double after topography (layer 0)
    int_inc_tmp = interface_increase
    grid_space_vert = grid_space_z

    # Stop at the bottom, no layers after bottom of mesh
    for i in range(len(all_interfaces) - 1):
        j = i + 1
        grid_space_vert *= int_inc_tmp[i]
        num_layers = (all_interfaces[j] - all_interfaces[i]) / grid_space_vert
        layers.append(myround(num_layers, 1, "near"))
        logger.info(f"\t{layers[i]} layers of {grid_space_vert}km between "
                    f"{all_interfaces[i]}km and {all_interfaces[j]}km "
                    f"({layers[i] * grid_space_vert} km)")

    # Get an even number of elements, place new layers on top
    while sum(layers) < myround(sum(layers), 2, "up"):
        layers[0] += 1

    # Start counting layers from bottom
    layers = layers[::-1]

    logger.info(f"\tNLAYERS = {len(layers)}")
    logger.info(f"\tNELEMENTS_Z = {sum(layers)}")
    logger.info(f"\tNUM ELEM PER LAYER FROM BOTTOM = {layers}")

    return layers


def approx_element_number(depth_km, layers_from_top, top, bottom, 
                          grid_space_top, interface_increase=2):
    """
    Return an approximate element number based on the depth, the grid space
    at the top of the mesh, and the known interfaces.

    NOTE: Element numbering starts from the bottom!

    :type depth_km: float
    :param depth_km: depth of the requested element, depth increase positive
    :type interface_increase: int
    :param interface_increase: the amount to multiply vertical element spacing
        at each interface
    
    :rtype: int
    :return: approximate element
    """
    num_elem = sum(layers_from_top)

    # Return mesh bottom if requested
    if depth_km == bottom:
        logger.info(f"\tdepth {depth_km}km is bottom of mesh {bottom}")
        return 1
    elif depth_km == top:
        logger.info(f"\tdepth {depth_km}km is top of mesh {top}km")
        return num_elem

    # Ensure we don't change the grid_space_z variable
    grid_space_vertical = grid_space_top

    # Reverse the Z axis so increasing depth is positive
    current_depth = top

    # Start counting down from the top layer
    element = num_elem
    for num_lay in layers_from_top:
        # Count through the number of layer in each grouping
        for _ in range(1, num_lay + 1):
            current_depth += grid_space_vertical
            if current_depth > depth_km:
                logger.info(f"\tdepth {depth_km}km is approx element "
                            f"{element} at depth {current_depth:.2f}km")
                return element
            else:
                element -= 1
        # The next layer means vertical increase (two-fold, three-fold etc.)
        # grid_space_vertical *= 3  # ORIGINAL
        grid_space_vertical *= interface_increase 

    logger.warning("You weren't supposed to see this... sorry, goodbye")
    sys.exit()


def nmaterials_nregions_ndoublings(doubling_layers, regions, layers, nex_xi,
                                   nex_eta, top, bottom, grid_space_z,
                                   materials, interface_increase=2):
    """
    Define the string that gives the doubling layers which control the
    horizontal doubling of element area, as well as the strings that define
    nmaterials and nregions which control different material properties of
    different sections of the mesh.

    There should atleast be one doubling layer, one material and one region.
    """
    # Ensure that we start from the top layer
    layers_from_top = layers[::-1]
    cumulative_layers = np.cumsum(layers)

    logger.info("FORMATTING HORIZONTAL DOUBLING LAYERS (ndoublings)")
    db_fmt = ""
    db_template = "NZ_DOUBLING_{j}                   = {value}\n"
    doubling_layers.sort(reverse=True)
    nz_per_doubling_layer = []
    for i, dl in enumerate(doubling_layers):
        elem_num = approx_element_number(depth_km=dl, 
                                         layers_from_top=layers_from_top, 
                                         top=top, bottom=bottom,
                                         grid_space_top=grid_space_z,
                                         interface_increase=interface_increase,)
        # This isnt necessary, matched interface and doubling layers is good
        # if elem_num in cumulative_layers:
        #     logger.info("\t!!! doubling layer matches interface layer, "
        #                 "placing above interface")
        #     elem_num -= 1
        nz_per_doubling_layer.append(elem_num)
        db_fmt += db_template.format(j=i+1, value=elem_num)

    # Doubling layers can't be next to one another    
    nz_per_dbl = np.array(nz_per_doubling_layer)
    layer_differences = abs(nz_per_dbl[1:] - nz_per_dbl[:-1])
    assert(1 not in layer_differences), \
                        "Doubling layers must be more than 1 element apart"

    logger.info(f"\tNDOUBLING LAYERS = {nz_per_doubling_layer}")

    # Format the reigons and materials files
    logger.info("FORMATTING NREGIONS AND NMATERIALS")
    # material_id, rho, vp, vs, Qk, Qm, anisotropy_flag, domain_id
    materials_template = "{j}  {rho}  {vp}  {vs}  {qk}  {qm}  0  2\n"
    # nex_xi_beg, nex_xi_end, nex_eta_beg, nex_eta_end, nz_beg, nz_end, mat_id
    regions_template = ("{nex_begin:<5} {nex_xi:<5} {nex_begin:<5} "
                        "{nex_eta:<5} {nz_beg:<5} {nz_end:<5} {j:<5}\n")

    # Regions define interfaces so nregions is 1 additional
    nregions = nmaterials = len(regions) + 1
    logger.info(f"\tnregions = nmaterials = {nregions}")

    # Ensure we start counting from the bottom of the mesh, element number 1
    nz_begin = 1

    # Include the top of the mesh in the regions
    regions.sort(reverse=True)
    regions += [top]

    # Determine the material parameters that will be placed in each region
    input_materials_default = [5500, 7500, 4500, 9999., 1000.]
    if materials is None:
        input_materials = []
        # Dynamically create a list of lists to match number of regions + 1
        for i in range(len(regions) + 1):
            input_materials.append(input_materials_default)
    else:
        input_materials = materials

    materials_out, regions_out, rgn_elem = "", "", []
    for i, reg in enumerate(regions):
        j = i + 1
        # Materials don't need to be specific as they will be overwritten
        # by external tomography files. Otherwise User will have to manual set

        rho, vp, vs, qk, qm = input_materials[i]
        materials_out += materials_template.format(j=j, rho=rho, vp=vp, vs=vs,
                                                   qk=qk, qm=qm)
        # Step down the values of the input materials so that they still remain
        # physical, otherwise meshfem will whine
        # !!! Commented out for Ristau1D model, manual input material 
        # input_materials = [_-200 for _ in input_materials]
        # for _ in input_materials:
        #     assert(_ > 0), \
        #         "Number of material layers has caused default values to be neg."

        # Regions are described by the interfaces of the external tomo files
        nz_end = approx_element_number(depth_km=reg, 
                                       layers_from_top=layers_from_top, 
                                       top=top, bottom=bottom,
                                       grid_space_top=grid_space_z,
                                       interface_increase=interface_increase)
        
        regions_out += regions_template.format(nex_begin="1",
                                               nex_xi=nex_xi,
                                               nex_eta=nex_eta,
                                               nz_beg=nz_begin,
                                               nz_end=nz_end,
                                               j=j)
        rgn_elem.append(nz_end)
    
        # Ensure that the next layer starts on the next element
        nz_begin = nz_end + 1

    rgn_elem.sort(reverse=True)
    logger.info(f"\tREGION STARTING ELEMENT = {rgn_elem}")

    return db_fmt, nz_per_doubling_layer, nregions, rgn_elem, nmaterials, \
           regions_out, materials_out


def write_interfaces(template, dir_name, layers, interfaces, lat_min, lon_min,
                     suppress_utm_proj, fids=[]):
    """
    Write the interfaces.dat file as well as the corresponding flat interface
    layers. Topo will need to be written manually

    :type template: str
    :param template: fid of the interfaces.dat template, preformatted
    :type dir_name: str
    :param dir_name: directory to store output files
    :type layers: list
    :param layers: number of elements in each layer
    :type interfaces: list
    :param interfaces: depths in km of each interface
    :param lat_min: minimum latitude of the mesh
    :param lon_min: minimum longitude of the mesh
    :type suppress_utm_proj: bool
    :param suppress_utm_proj: if True, will use Lat/Lon coordinates, else
        uses a UTM projection conversion
    :type fids: list
    :param fids: for custom interface names, starting from the top going down,
        excluding topography and the bottom of the mesh. If left blank, defaults
        to 1, 2 ... until the number of layers
    """
    # Ensure counting layers from bottom
    layers_from_bottom = layers

    # Template for setting a flat layer
    flat_layer = (f".{str(suppress_utm_proj).lower()}. 2 2 {lon_min:.1f}d0 "
                  f"{lat_min:.1f}d0 180.d0 180.d0")

    # Write to a new file
    with open(os.path.join(dir_name, "interfaces.dat"), "w") as f:
        f.write("# number of interfaces\n")
        f.write(" {ninterfaces}\n".format(ninterfaces=len(interfaces)))

        # Write flat layers up to topo
        for i, (interface, fid) in enumerate(zip(interfaces, reversed(fids))):
            f.write("# interface number {}\n".format(i+1))
            f.write(f" {flat_layer}\n")
            f.write(f" {fid}\n")

        # Write number of elements in each layer
        f.write("# number of spectral elements per layer (from bottom)\n")
        for j, layer in enumerate(layers_from_bottom):
            f.write(f" {layer}\n")

    # Write the individual interface files
    for fid, interface in zip(fids, interfaces):
        logger.info(f"WRITING INTERACE {fid}")
        with open(os.path.join(dir_name, fid), "w") as f:
            for i in range(4):
                f.write("{}\n".format(-1 * abs(int(interface * 1E3))))


def pprint_mesh_stats(interfaces, doublings, regions, dx, dy, dz, top,
                      interface_increase=2):
    """
    A printing function that shows element numbers, sizes and doubling layers
    and interfaces, to get an idea of what the mesh will look like without
    having to run the xmeshfem3D binary
    """
    logger.info("MESH HAS THE APPROXIMATE FORM")
        
    # Total number of vertical layers in the mesh
    layers_nz = sum(interfaces)
    # Figure out the element numbers corresponding to interfaces
    interfaces = np.cumsum(interfaces)
    c, r = 0, 1
    dx_, dy_, dz_ = dx, dy, dz
    depth = top
    int_idx = 0  # interface index 
    for i in range(layers_nz, 0, -1):
        j = layers_nz - i 
        msg = ""
        if i == layers_nz:
            msg += "\n\t\tTop of Mesh"
        if i in interfaces:
            msg += f"\n\t\tInterface (vertical x{interface_increase[int_idx]})"
            # dz *= 3  # ORIGINAL
            # Increment the interface increase number
            dz *= interface_increase[int_idx]
            int_idx += 1
        if i in doublings:
            msg += "\n\t\tDoubling (horizontal doubling)"
            dx *= 2
            dy *= 2
        if i in regions:
            msg += f"\n\t\tRegion {r} begins"
            r += 1
        logger.info(f"\t{i:0>2}/{j:0>2}; Z~={depth:6.2f}km")
        depth -= dz
        if msg:
            logger.info(msg[1:])
            logger.info(f"\t\tdx={dx:.2f}, dy={dy:.2f}, dz={dz}")
        

def prepare_meshfem(parameter_file, mesh_par_file_template,
                    interfaces_template):
    """
    Run all the above functions and output the necessary meshfem3D files
    """
    assert(os.path.exists(parameter_file))
    assert (os.path.exists(mesh_par_file_template))
    assert (os.path.exists(interfaces_template))

    pars = set_parameters(parameter_file)

    try:
        # Put the outputs into a directory for easy transfer
        if not os.path.exists(pars["dir_name"]):
            os.makedirs(pars["dir_name"])

        # Determine the grid space based on wavespeed/min period or by user set
        grid_space_h, grid_space_v = minimum_grid_spacing(
            grid_space_top_hv=pars["grid_space_top_hv"],
            slowest_wavespeed=pars["vs_min_km_per_s"],
            shortest_period=pars["shortest_period_s"],
            top_vertical_element_ratio=pars["top_vertical_element_ratio"],
            points_per_wl=pars["points_per_wl"]
                                          )

        # Number of processors controlled by given User set parameters
        nproc_x, nproc_y = number_of_processors(nproc=pars["nproc"],
                                                x_length=pars["x_length_km"],
                                                y_length=pars["y_length_km"]
                                                )

        # Number of elements controlled by Specfem requirements
        nex_x, nex_y, dx, dy = number_of_elements(
            nproc_x=nproc_x, nproc_y=nproc_y, x_length=pars["x_length_km"],
            y_length=pars["y_length_km"], grid_space=grid_space_h,
            shortest_period_s=pars["shortest_period_s"], 
            ndoublings=len(pars["doubling_layers"])
        )

        # Adjust vertical element spacing to new grid spacing, which is new
        if pars["top_vertical_element_ratio"] is not None:
            grid_space_v = adjust_vertical_element_size(
                    dx=dx, dy=dy, top_vertical_element_ratio=pars[
                        "top_vertical_element_ratio"]
                    )

        # Interfaces control vertical height of elements
        layers = interface_layers(
            top_of_mesh=pars["mesh_top_km"], depth=pars["mesh_depth_km"],
            interfaces=pars["interfaces"], grid_space_z=grid_space_v,
            interface_increase=pars["interface_increase"]
                                               )

        # Write out the specific format for doubling, materials and regions
        dbl_str, nz, nregions, regions, nmaterials, rgn_str, mat_str = \
            nmaterials_nregions_ndoublings(
                doubling_layers=pars["doubling_layers"], layers=layers,
                regions=pars["regions"], nex_xi=nex_x, nex_eta=nex_y,
                top=pars["mesh_top_km"], bottom=pars["mesh_depth_km"],
                grid_space_z=grid_space_v, materials=pars["materials"]
        )

        # Print the mesh stats for easier qualification of mesh
        pprint_mesh_stats(interfaces=layers, doublings=nz, regions=regions,
                          dx=dx, dy=dy, dz=grid_space_v, 
                          top=pars["mesh_top_km"],
                          interface_increase=pars["interface_increase"])

        # Format the template Mesh_Par_file
        logger.info("WRITING Mesh_Par_file")
        with open(mesh_par_file_template, "r") as f:
            lines = f.read()
        lines = lines.format(lat_min=pars["lat_min"], lat_max=pars["lat_max"],
                             lon_min=pars["lon_min"], lon_max=pars["lon_max"],
                             depth=pars["mesh_depth_km"],
                             utm_projection=pars["utm_projection"],
                             suppress_utm_proj=str(
                                    pars["suppress_utm_proj"]).lower(),
                             nex_xi=nex_x, nex_eta=nex_y, nproc_xi=nproc_x,
                             nproc_eta=nproc_y,
                             ndoublings=len(pars["doubling_layers"]),
                             doubling_layers=dbl_str,
                             nmaterials=nmaterials, nregions=nregions,
                             materials=mat_str, regions=rgn_str
                             )
        with open(os.path.join(pars["dir_name"], "Mesh_Par_file"), "w") as f:
            f.write(lines)

        # Format the interfaces file
        if pars["interfaces"]:
            logger.info("WRITING interfaces.dat")
            write_interfaces(template=interfaces_template,
                             dir_name=pars["dir_name"], layers=layers, 
                             interfaces=pars["interfaces"],
                             suppress_utm_proj=pars["suppress_utm_proj"],
                             lat_min=pars["lat_min"], lon_min=pars["lon_min"],
                             fids=pars["interface_fids"]
                             )

        calculate_nelements(nex_x, nex_y, layers, nz)
    except Exception as e:
        traceback.print_exc()
        pass

    # Hacky way to get the output log in the right name
    shutil.move("log_mesh_temp.out", f"{pars['dir_name']}/{pars['tag']}.out")
    shutil.copy(parameter_file, 
                f"{pars['dir_name']}/{os.path.basename(parameter_file)}")


if __name__ == "__main__":
    # Allow custom parameter file names
    try:
        parfile = sys.argv[1]
    except IndexError:
        parfile = "parmesh.json"
    prepare_meshfem(parameter_file=parfile,
                    mesh_par_file_template="./template_Mesh_Par_file",
                    interfaces_template="./template_interfaces.dat",
                    )
