"""
A script to generate the necessary file for Geocubit
"""
import numpy as np
import sys
sys.path.append('..')
from mesh_utils import myround, lonlat_utm


def set_parameters():
    """
    Define the necessary parameters here, these will be accessed by the
    constraint fuctions throughout the script. Defaults for a New Zealand North
    coarse mesh are set here.

    :rtype: dict
    :return: dictionary of parameters
    """
    parameters = {"x_or_lon_min": 125E3,
                  "x_or_lon_max": 725E3,
                  "y_or_lat_min": 5150E3,
                  "y_or_lat_max": 5950E3,
                  "utm_projection": -60,
                  "mesh_depth_km": 400.,
                  "interfaces": [8, 33, 100],
                  "nproc": 40,
                  "shortest_period_s": 10,
                  "vs_min_km_per_s": 1.,
                  "working_dir": "./working_dir",
                  "output_dir": "./output_dir",
                  "topo_fid": "./topo_utm60_x600_y800_4km.xyz",
                  "moho_fid": "./moho33_utm60_x600_y800_4km.xyz"
                  }
    
    # Set ntriplings
    parameters["ntriplings"] = len(parameters["interfaces"])
        
    # Convert latlon to UTM
    if parameters["x_or_lon_min"] < 1E3:
        x_min, y_min = lonlat_utm(parameters["lon_min"], parameters["lat_min"],
                                  parameters["utm_projection"])
        x_max, y_max = lonlat_utm(parameters["lon_max"], parameters["lat_max"],
                                  parameters["utm_projection"])

        # Set the UTM coordinates in the parameters
        parameters["x_min"] = x_min
        parameters["y_min"] = y_min
        parameters["x_max"] = x_max
        parameters["y_max"] = y_max

    # If mesh dimensions already given in UTM projection
    else:
        parameters["x_min"] = parameters["x_or_lon_min"]
        parameters["x_max"] = parameters["x_or_lon_max"]
        parameters["y_min"] = parameters["y_or_lat_min"]
        parameters["y_max"] = parameters["y_or_lat_max"]

    # Determine the length of each mesh dimension
    parameters["x_length_km"] = (
            (parameters["x_max"] - parameters["x_min"]) * 1E-3
    )
    parameters["y_length_km"] = (
            (parameters["y_max"] - parameters["y_min"]) * 1E-3
    )
    return parameters


def minimum_grid_spacing(slowest_wavespeed, shortest_period):
    """
    Define minimum grid spacing based on slowest wavespeed and shortest period

    :type slowest_wavespeed: float
    :param slowest_wavespeed: slowest wavespeed expected in model
    :type shortest_period: float
    :param shortest_period: shortest period to be resolved in mesh
    :rtype: float
    :return: minimum element spacing of the mesh
    """
    # Value of 2 comes from two points per wavelength
    min_grid_space = (shortest_period * slowest_wavespeed) / 2
    print(f"\tminimum grid spacing calculated as {min_grid_space}")
    grid_space = myround(min_grid_space, 2, 'down')
    print(f"\t\trounded down to {grid_space}")

    return grid_space


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
    # Determine the ratio
    ratio = min(x_length, y_length) / max(x_length, y_length)

    # Define which direction is shorter
    if x_length < y_length:
        short_direction = "x"
    else:
        short_direction = "y"

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
            print(f"\tNPROC_X = {nproc_x}\n\tNPROC_Y = {nproc_y}")
            return int(nproc_x), int(nproc_y)
        else:
            guess_a += 1


def number_of_elements(nproc_x, nproc_y, x_length, y_length, grid_space,
                       shortest_period_s):
    """
    Define the number of elements in each horizontal direction based on number
    of processors

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
    # meshfem requires the number of grid points be an integer multiple of 8
    # times the number of processors in a given direction
    nex_x = myround(x_length / grid_space, nproc_x * 8, "near")
    nex_y = myround(y_length / grid_space, nproc_y * 8, "near")

    # ensure that the short direction is maintained
    while (nproc_x < nproc_y) and (nex_x > nex_y):
        nex_y += nproc_y * 8

    # Minimum element number set by Eq. 3.1 in Specfem manual
    assert(nex_x >= 2 * 288 / shortest_period_s)

    dx = x_length / nex_x
    dy = y_length / nex_y

    print(f"\tNEX_X = {nex_x}\n\tNEX_Y = {nex_y}\n"
          f"\tdx = {dx:.2f}km\n\tdy = {dy:.2f}km"
          )

    return int(nex_x), int(nex_y)


def vertical_tripling_proportions(depth, ntriplings, interfaces, grid_space):
    """
    Set the number of vertical doubling layers

    :param depth:
    :param ntriplings:
    :param interfaces:
    :return:
    """
    # Start interfaces from the top, include the bottom interface of depth
    all_interfaces = [0] + interfaces + [depth]
    all_interfaces.sort()

    layers = []
    for i in range(ntriplings + 1):
        j = i + 1
        num_layers = ((all_interfaces[j] - all_interfaces[i]) /
                      ((3 ** i) * grid_space))
        layers.append(myround(num_layers, 1, "near"))

    # Get an even number of elements, place new layers on top
    while sum(layers) < myround(sum(layers), 2, "up"):
        layers[0] += 1

    # Start counting layers from bottom
    layers.sort(reverse=True)

    print("\t{nlay} sections; bottom to top {layers} elements per layer".format(
        nlay=len(layers), layers=layers))
    print("\t{} total layers".format(sum(layers)))

    return layers


def write_config(template, lat_min, lat_max, lon_min, lon_max,
                 depth, nex_x, nex_y, nproc_x, nproc_y,
                 ntriplings, layers, grid_space, working_dir,
                 output_dir, topo_fid, moho_fid):
    """
    Write the Mesh_Par_file with the given values and a template script

    TO DO:
        add the capability to dynamically set the number of doubling layers

    :type template: str
    :param template: fid of the Mesh_Par_file template, preformatted
    :param lat_min: minimum latitude of the mesh
    :param lat_max: maximum latitude of the mesh
    :param lon_min: maximum longitude of the mesh
    :param lon_max: maximum longitude of the mesh
    :param depth: depth of the mesh in km
    :param nex_x: number of elements in the x direction
    :param nex_y: number of elements in the y direction
    :param nproc_x: number of processors in the x direction
    :param nproc_y: number of processors in the y direction
    :param ntriplings: number of doubling layers
    :type layers: list
    :param layers: number of elements in each layer

    """
    with open(template, "r") as f:
        lines = f.read()

    # Format the list of tripling layers
    ntripl_layers = ""
    for layer in layers[:-1]:
        ntripl_layers += f"{layer},"
    ntripl_layers = ntripl_layers[:-1]

    # Get the approximate size of a bottom element
    bottom_element = grid_space * 3 ** ntriplings

    lines = lines.format(lat_min=lat_min, lat_max=lat_max, lon_min=lon_min,
                         lon_max=lon_max, depth_km=-1E3 * depth,
                         nex_xi=nex_x, nex_eta=nex_y, nproc_xi=nproc_x,
                         nproc_eta=nproc_y, ntriplings=ntriplings,
                         ntripl_layers=ntripl_layers,
                         bottom_element=bottom_element * 1E3,
                         working_dir=working_dir, output_dir=output_dir,
                         nlayer=ntriplings + 1, topo_fid=topo_fid,
                         moho_fid=moho_fid
                         )

    print("\twriting Config")
    with open("./geocubit_config.cfg", "w") as f:
        f.write(lines)


def prepare_meshfem():
    """
    Run all the above functions and output the necessary meshfem3D files
    :return:
    """
    config_template = "./template_config.cfg"

    print("Preparing Geocubit config")
    pars = set_parameters()
    grid_space = minimum_grid_spacing(slowest_wavespeed=pars["vs_min_km_per_s"],
                                      shortest_period=pars["shortest_period_s"]
                                      )
    nproc_x, nproc_y = number_of_processors(nproc=pars["nproc"],
                                            x_length=pars["x_length_km"],
                                            y_length=pars["y_length_km"]
                                            )
    nex_x, nex_y = number_of_elements(
        nproc_x=nproc_x, nproc_y=nproc_y, x_length=pars["x_length_km"],
        y_length=pars["y_length_km"], grid_space=grid_space,
        shortest_period_s=pars["shortest_period_s"]
    )
    layers = vertical_tripling_proportions(depth=pars["mesh_depth_km"],
                                           ntriplings=pars["ntriplings"],
                                           interfaces=pars["interfaces"],
                                           grid_space=grid_space
                                           )
    write_config(template=config_template,
                 lat_min=pars["y_min"], lat_max=pars["y_max"],
                 lon_min=pars["x_min"], lon_max=pars["x_max"],
                 depth=pars["mesh_depth_km"],
                 nex_x=nex_x, nex_y=nex_y, nproc_x=nproc_x,
                 nproc_y=nproc_y, ntriplings=pars["ntriplings"],
                 layers=layers, grid_space=grid_space,
                 working_dir=pars["working_dir"],
                 output_dir=pars["output_dir"],
                 topo_fid=pars["topo_fid"], moho_fid=pars["moho_fid"]
                 )


if __name__ == "__main__":
    prepare_meshfem()
