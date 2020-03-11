"""
A script to generate the necessary files for Specfem3D Cartesians internal
mesher, Meshfem3D. So far hardcoded for two doubling layers, both in the
templates and in the functions below. If more or less layers are required,
one would need to edit the templates to add more layers and interfaces, and then
change the string formatting in the write functions below

NOTE:
1) Doubling layers control the horizontal doubling of element area
2) Interfaces control vertical doubling of element length

Parameters must be set in the meshpar.json file (or user defined file)
The JSON parameter file can contain the following parameters

:type tag: str
:param tag: tag to save outputs to
:type lat_min: float
:param lat_min: minimum latitude in degrees
:type lat_max: float
:param lat_max: maximum latitude in degrees
:type lon_min: float
:param lon_min: minimum longitude in degrees
:type lon_max: float
:param lon_max: maximum longitude in degrees
:type utm_projection: int
:param utm_prjection: UTM projection to convert lat/lon with
:type mesh_depth_km: float
:param mesh_depth_km: mesh depth in units of kilometers
:type interfaces: list of float
:param interfaces: the depth locations of interfaces where doubling occurs
  in units of kilometers
:type interface_fids: list of str
:param interface_fids: name of the output interface files corresponding to
  the list of interfaces
:type mantle_from_interface_idx: int
:param mantle_from_interface_idx: allows the mantle to have doubling layers
  while retaining all the other properties.
:type nproc: int
:param nproc: number of processors to split up the mesh into
:type shortest_period_s: float
:param shortest_period_s: the shortest period that the mesh is expected to
  resolve. Not a strict requirement, but gives an idea of the necessary
  element sizes required to resolve given features. A minimum
  resolvable period test must be done to get the actual resolution of the mesh
:type vs_min_km_per_s: float
:param vs_min_km_per_s: smallest expected wavespeed in the model, in km/s
  as with `shortest_period_s`, gives an idea of the necessary element sizes
  but is not a hard requirement. An MRP test will need to be performed
"""
import os
import sys
import json
import logging
import traceback
import numpy as np

# Import from the utilities above
import sys
sys.path.append("..")
from mesh_utils import myround, lonlat_utm


logging.basicConfig(level=logging.DEBUG,
                    format="%(message)s",
                    filename=f"log_mesh_temp.out",
                    filemode="w")

# Write logging to console
console = logging.StreamHandler()
console.setLevel(logging.INFO)
logging.getLogger('').addHandler(console)
logger = logging.getLogger("mesher")


def set_parameters(fid="./parmesh.json"):
    """
    Define the necessary parameters here, these will be accessed by the
    constraint fuctions throughout the script. Defaults for a New Zealand North
    coarse mesh are set here.

    :rtype: dict
    :return: dictionary of parameters
    """
    with open(fid, "r") as f:
        parameters = json.load(f)

    # Print the parameters so the User knows what they've chosen
    logger.info(f"PREPARING MESH FOR MESHFEM3D")
    logger.info(f"\n{parameters}\n")
    logger.info(f"PARAMETERS FOR: {parameters['tag']}")
    for key, item in parameters.items():
        logger.info(f"\t{key}: {item}")

    parameters["ndoublings"] = len(parameters["doubling_layers"])

    # Make sure there are enough fids for each interface, or that none specified
    assert((len(parameters["interface_fids"]) == 0) or
           (len(parameters["interfaces"]) == len(parameters["interface_fids"]))
           ), "interfaces number must match interface_fids number"

    # Convert latlon to UTM
    x_min, y_min = lonlat_utm(parameters["lon_min"], parameters["lat_min"],
                              parameters["utm_projection"])
    x_max, y_max = lonlat_utm(parameters["lon_max"], parameters["lat_max"],
                              parameters["utm_projection"])

    # Set the UTM coordinates in the parameters
    parameters["x_min"] = x_min
    parameters["y_min"] = y_min
    parameters["x_max"] = x_max
    parameters["y_max"] = y_max
    parameters["x_length_km"] = (x_max - x_min) * 1E-3
    parameters["y_length_km"] = (y_max - y_min) * 1E-3

    return parameters


def minimum_grid_spacing(slowest_wavespeed, shortest_period,
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

        # Value of 2 comes from two points per wavelength
        min_grid_space = (shortest_period * slowest_wavespeed) / 2
        logger.info(f"\tT_min/V_min = {min_grid_space}")
        grid_space_h = myround(min_grid_space, 2, 'down') or min_grid_space
        grid_space_v = grid_space_h / 2
        if grid_space_h == myround(min_grid_space, 2, 'down'):
            logger.info(f"\trounded to nearest factor of 2")
        logger.info(f"\tMIN GRID SPACE HORIZONTAL = {grid_space_h} km")
        logger.info(f"\tMIN GRID SPACE VERTICAL = {grid_space_v} km")
        
        return grid_space_h, grid_space_v


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

    logger.info(f"\tshort direction is {short_direction}")

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


def number_of_elements(nproc_x, nproc_y, x_length, y_length, grid_space,
                       shortest_period_s):
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
    # meshfem requires the number of grid points be an integer multiple of 8
    # times the number of processors in a given direction
    c = 1
    nex_x = 0
    while nex_x < 2 * 288 / shortest_period_s:
        nex_x = myround(x_length / grid_space, c * nproc_x * 8, "up")
        nex_y = myround(y_length / grid_space, c * nproc_y * 8, "up")
        c += 1
        if c >= 50:
            raise OverflowError("nex_x seems too large for given nproc_x")

    # ensure that the short direction is maintained
    while (nproc_x < nproc_y) and (nex_x > nex_y):
        nex_y += nproc_y * 8

    logger.info(f"\tinteger multiple found as c={c}")
    logger.info(f"\tNEX_X = {nex_x}\n"
                f"\tNEX_Y = {nex_y}")

    assert(nex_x >= 2 * 288 / shortest_period_s)
    assert(nex_y >= 2 * 288 / shortest_period_s)

    # Minimum element number set by Eq. 3.1 in Specfem manual
    dx = x_length / nex_x
    dy = y_length / nex_y

    logger.info(f"\tDX = {dx:.4f}km\n"
                f"\tDY = {dy:.4f}km")

    return int(nex_x), int(nex_y)


def calculate_nelements(nex_x, nex_y, layers, ):
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
    logger.info("CALCULATING NUMBER OF ELEMENTS")
    nelements = 0
    layers_from_top = layers
    layers_from_top.sort()
    for i, nz in enumerate(layers_from_top):
        j = i + 1
        elements_in_layer = (nex_x * nex_y) / (2 * j) * nz
        logger.info(f"\t{elements_in_layer} in layer {j}")
        nelements += elements_in_layer

    nelements = int(nelements)
    logger.info(f"\tNUM_ELEM = {nelements}")
    
    return nelements


def vertical_doubling_proportions(top_of_mesh, depth, interfaces, grid_space_z):
    """
    Set the number of vertical doubling layers
    :param depth:
    :param ndoublings:
    :param interfaces:
    :return:
    """
    logger.info("CALCULATING NUMBER OF VERTICAL DOUBLING LAYERS (interfaces)")
    # Start interfaces from the top, include the bottom interface of depth
    all_interfaces = [top_of_mesh] + interfaces + [depth]
    all_interfaces.sort()

    layers = []
    for i in range(len(interfaces) + 1):
        j = i + 1
        num_layers = ((all_interfaces[j] - all_interfaces[i]) /
                      ((2 ** i) * grid_space_z))
        layers.append(myround(num_layers, 1, "near"))

    # Get an even number of elements, place new layers on top
    while sum(layers) < myround(sum(layers), 2, "up"):
        layers[0] += 1

    # Start counting layers from bottom
    layers.sort(reverse=True)

    logger.info(f"\tNLAYERS = {len(layers)}")
    logger.info(f"\tNELEMENTS_Z = {sum(layers)}")
    logger.info(f"\tLAYER SIZE FROM BOTTOM = {layers}")

    return layers


def approx_element_number(depth_km, layers, top, bottom, grid_space_top):
    """
    Return an approximate element number based on the depth, the grid space
    at the top of the mesh, and the known vertical doubling layers.

    NOTE: Element numbering starts from the bottom!

    :type depth_km: float
    :param depth_km: depth of the requested element, depth increase positive
    :rtype: int
    :return: approximate element
    """
    num_elem = sum(layers)

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
    for num_lay in layers:
        # Count through the number of layer in each grouping
        for _ in range(1, num_lay + 1):
            current_depth += grid_space_vertical
            if current_depth > depth_km:
                logger.info(f"\tdepth {depth_km}km is roughly element "
                            f"{element} at depth {current_depth}km")
                return element
            else:
                element -= 1
        # The next layer means vertical doubling
        grid_space_vertical *= 2

    logger.warning("You weren't supposed to see this... sorry")
    sys.exit()


def nmaterials_nregions_ndoublings(doubling_layers, regions, layers, nex_xi,
                                   nex_eta, top, bottom, grid_space_z):
    """
    Define the string that gives the doubling layers which control the
    horizontal doubling of element area, as well as the strings that define
    nmaterials and nregions which control different material properties of
    different sections of the mesh.

    There should atleast be one doubling layer, one material and one region.
    """
    # Ensure that we start from the top layer
    layers.sort()

    logger.info("FORMATTING HORIZONTAL DOUBLING LAYERS (ndoublings)")
    db_fmt = ""
    db_template = "NZ_DOUBLING_{j}                   = {value}\n"
    doubling_layers.sort(reverse=True)
    nz_per_doubling_layer = []
    for i, dl in enumerate(doubling_layers):
        elem_num = approx_element_number(depth_km=dl, layers=layers, top=top,
                                         bottom=bottom,
                                         grid_space_top=grid_space_z)
        nz_per_doubling_layer.append(elem_num)
        db_fmt += db_template.format(j=i+1, value=elem_num)
    logger.info(f"\tNDOUBLING LAYERS = {nz_per_doubling_layer}")

    # Format the reigons and materials files
    logger.info("FORMATTING NREGIONS AND NMATERIALS")
    # material, id, rho, vp, vs, Qk, Qm, anisotropy_flag, domain_id
    materials_template = "{j}  9999.  9999.  9999.  9999.  9999.  0  2\n"
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

    materials_out = ""
    regions_out = ""
    # Use range because
    for i, reg in enumerate(regions):
        j = i + 1
        # Materials don't need to be specific as they will be overwritten
        # by external tomography files. Otherwise User will have to manual set
        materials_out += materials_template.format(j=j)

        # Regions are described by the interfaces of the external tomo files
        nz_end = approx_element_number(depth_km=reg, layers=layers, top=top,
                                       bottom=bottom,
                                       grid_space_top=grid_space_z)
        regions_out += regions_template.format(nex_begin="1",
                                               nex_xi=nex_xi,
                                               nex_eta=nex_eta,
                                               nz_beg=nz_begin,
                                               nz_end=nz_end,
                                               j=j)

        # Ensure that the next layer starts on the next element
        nz_begin = nz_end + 1

    return db_fmt, nregions, nmaterials, regions_out, materials_out


def write_interfaces(template, dir_name, layers, interfaces, lat_min, lon_min,
                     fids=[]):
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
    :type fids: list
    :param fids: for custom interface names, starting from the top going down,
        excluding topography and the bottom of the mesh. If left blank, defaults
        to 1, 2 ... until the number of layers
    """
    # Ensure counting layers from bottom
    layers_from_bottom = layers
    layers_from_bottom.sort(reverse=True)

    # Template for setting a flat layer
    flat_layer = f".false. 2 2 {lon_min:.1f}d0 {lat_min:.1f}d0 180.d0 180.d0"

    # Hardcoded topo layer, CHANGE THIS
    topo_fid = "interface_topo.dat"
    topo = ".false. 720 720 173.d0 -43.d0 0.00833d0 0.00833d0"

    # Write to a new file
    logger.info("WRITING interfaces.dat")
    with open(os.path.join(dir_name, "interfaces.dat"), "w") as f:
        f.write("# number of interfaces\n")
        f.write(" {ninterfaces}\n".format(ninterfaces=len(interfaces) + 1))

        # Write flat layers up to topo
        for i, (interface, fid) in enumerate(zip(interfaces, reversed(fids))):
            f.write("# interface number {}\n".format(i+1))
            f.write(f" {flat_layer}\n")
            f.write(f" {fid}\n")

        # Write topo interface
        f.write("# interface number {} (topo)\n".format(i+2))
        f.write(f" {topo}\n")
        f.write(f" {topo_fid}\n")

        # Write number of elements in each layer
        f.write("# number of spectral elements per layer (from bottom)\n")
        for j, layer in enumerate(layers_from_bottom):
            f.write(f" {layer}\n")

    # Write the individual interface files
    for fid, interface in zip(fids, interfaces):
        # Skip top interface (topography)
        logger.info(f"WRITING INTERACE {fid}")
        with open(os.path.join(dir_name, fid), "w") as f:
            for i in range(4):
                f.write("-{}\n".format(abs(int(interface * 1E3))))


def prepare_meshfem(parameter_file, mesh_par_file_template,
                    interfaces_template, dir_name):
    """
    Run all the above functions and output the necessary meshfem3D files
    """
    assert(os.path.exists(parameter_file))
    assert (os.path.exists(mesh_par_file_template))
    assert (os.path.exists(interfaces_template))
    # Put the outputs into a directory for easy transfer
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

    try:
        pars = set_parameters(parameter_file)

        # Determine the grid space based on wavespeed/min period or by user set
        grid_space_h, grid_space_v = minimum_grid_spacing(
            grid_space_top_hv=pars["grid_space_top_hv"],
            slowest_wavespeed=pars["vs_min_km_per_s"],
            shortest_period=pars["shortest_period_s"]
                                          )

        # Number of processors controlled by given User set parameters
        nproc_x, nproc_y = number_of_processors(nproc=pars["nproc"],
                                                x_length=pars["x_length_km"],
                                                y_length=pars["y_length_km"]
                                                )

        # Number of elements controlled by Specfem requirements
        nex_x, nex_y = number_of_elements(
            nproc_x=nproc_x, nproc_y=nproc_y, x_length=pars["x_length_km"],
            y_length=pars["y_length_km"], grid_space=grid_space_h,
            shortest_period_s=pars["shortest_period_s"]
        )

        # Doubling layers control doubling in height of elements
        layers = vertical_doubling_proportions(
            top_of_mesh=pars["mesh_top_km"], depth=pars["mesh_depth_km"],
            interfaces=pars["interfaces"], grid_space_z=grid_space_v
                                               )

        # Write out the specific format for doubling, materials and regions
        doubling_layers, nregions, nmaterials, regions, materials = \
            nmaterials_nregions_ndoublings(
                doubling_layers=pars["doubling_layers"], layers=layers,
                regions=pars["regions"], nex_xi=nex_x, nex_eta=nex_y,
                top=pars["mesh_top_km"], bottom=pars["mesh_depth_km"],
                grid_space_z=grid_space_v
        )

        # Format the template Mesh_Par_file
        logger.info("WRITING Mesh_Par_file")
        with open(mesh_par_file_template, "r") as f:
            lines = f.read()
        lines = lines.format(lat_min=pars["lat_min"], lat_max=pars["lat_max"],
                             lon_min=pars["lon_min"], lon_max=pars["lon_max"],
                             depth=pars["mesh_depth_km"],
                             utm_projection=pars["utm_projection"],
                             nex_xi=nex_x, nex_eta=nex_y, nproc_xi=nproc_x,
                             nproc_eta=nproc_y,
                             ndoublings=len(pars["doubling_layers"]),
                             doubling_layers=doubling_layers,
                             nmaterials=nmaterials, nregions=nregions,
                             materials=materials, regions=regions
                             )
        with open(os.path.join(dir_name, "Mesh_Par_file"), "w") as f:
            f.write(lines)

        # Format the interfaces file
        write_interfaces(template=interfaces_template, dir_name=dir_name,
                         layers=layers, interfaces=pars["interfaces"],
                         lat_min=pars["lat_min"], lon_min=pars["lon_min"],
                         fids=pars["interface_fids"]
                         )

        calculate_nelements(nex_x, nex_y, layers)
    except Exception as e:
        traceback.print_exc()
        pass

    # Hacky way to get the output log in the right name
    os.rename("log_mesh_temp.out", f"./{dir_name}/{pars['tag']}.out")


if __name__ == "__main__":
    prepare_meshfem(parameter_file="./parmesh.json",
                    mesh_par_file_template="./template_Mesh_Par_file",
                    interfaces_template="./template_interfaces.dat",
                    dir_name="./meshfem3D_files"
                    )
