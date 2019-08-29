"""
A script to generate the necessary files for Specfem3D Cartesians internal
mesher, Meshfem3D. So far hardcoded for two doubling layers, both in the
templates and in the functions below. If more or less layers are required,
one would need to edit the templates to add more layers and interfaces, and then
change the string formatting in the write functions below
"""
import os
import numpy as np


def set_parameters():
    """
    Define the necessary parameters here, these will be accessed by the
    constraint fuctions throughout the script. Defaults for a New Zealand North
    coarse mesh are set here.

    :rtype: dict
    :return: dictionary of parameters
    """
    parameters = {"lat_min": -43.,
                  "lat_max": -37.,
                  "lon_min": 173.,
                  "lon_max": 179.,
                  "utm_projection": -60,
                  "mesh_depth_km": 400.,
                  "interfaces": [8, 36, 100, 200],
                  "interface_fids": ["interface_crustshallow.dat",
                                     "interface_moho.dat",
                                     "interface_100km.dat",
                                     "interface_200km.dat"],
                  "mantle_from_interface_idx": 2,
                  "nproc": 1,
                  "shortest_period_s": 10,
                  "vs_min_km_per_s": 1.,
                  }

    parameters["ndoublings"] = len(parameters["interfaces"])

    # Make sure there are enough fids for each interface, or that none specified
    assert((len(parameters["interface_fids"]) == 0) or
           (len(parameters["interfaces"]) == len(parameters["interface_fids"]))
           )

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


def vertical_doubling_proportions(depth, ndoublings, interfaces, grid_space):
    """
    Set the number of vertical doubling layers

    :param depth:
    :param ndoublings:
    :param interfaces:
    :return:
    """
    # Start interfaces from the top, include the bottom interface of depth
    all_interfaces = [0] + interfaces + [depth]
    all_interfaces.sort()

    layers = []
    for i in range(ndoublings + 1):
        j = i + 1
        num_layers = ((all_interfaces[j] - all_interfaces[i]) /
                      ((2 ** i) * grid_space))
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


def write_mesh_par_file(template, lat_min, lat_max, lon_min, lon_max,
                        depth, utm_projection, nex_x, nex_y, nproc_x, nproc_y,
                        ndoublings, layers, mantle_from=None):
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
    :param utm_projection: utm projection to convert lat lon to
    :param nex_x: number of elements in the x direction
    :param nex_y: number of elements in the y direction
    :param nproc_x: number of processors in the x direction
    :param nproc_y: number of processors in the y direction
    :param ndoublings: number of doubling layers
    :type layers: list
    :param layers: number of elements in each layer
    :type mantle_from: int
    :param mantle_from: if given, sets the index of the mantle layer, allowing
        for internal interfaces/doubling layers within the mantle, but not
        changing the number of regions or materials
    """
    with open(template, "r") as f:
        lines = f.read()

    # Allows the mantle to contain interfaces/doubling layers, but still retain
    # the same region/material id within. The nz_3b is a hack to pick the right
    # end layer, but this is hardcoded to only have 3 layers, should be changed
    if mantle_from:
        nmaterials = nregions = mantle_from + 1
        nz_3b = sum(layers)
    else:
        nmaterials = nregions = ndoublings + 1
        nz_3b = sum(layers[:3])

    lines = lines.format(lat_min=lat_min, lat_max=lat_max, lon_min=lon_min,
                         lon_max=lon_max, depth=depth,
                         utm_projection=utm_projection,
                         nex_xi=nex_x, nex_eta=nex_y, nproc_xi=nproc_x,
                         nproc_eta=nproc_y, ndoublings=ndoublings,
                         nz_doubling_1=sum(layers[:1]),
                         nz_doubling_2=sum(layers[:2]),
                         nmaterials=nmaterials, nregions=nregions,
                         nz_1a=1,
                         nz_1b=sum(layers[:1]),
                         nz_2a=sum(layers[:1]) + 1,
                         nz_2b=sum(layers[:2]),
                         nz_3a=sum(layers[:2]) + 1,
                         nz_3b=nz_3b,
                         )

    # Hacky way to incorporate different doubling layers
    print("\n!!! You must add the following to the Mesh_Par_file !!!")
    for i, layer in enumerate(layers[:-1]):
        j = i + 1
        print('NZ_DOUBLING_{}\t\t\t\t\t= {}'.format(j, sum(layers[:j])))
    print("!!! You must add the following to the Mesh_Par_file !!!\n")

    # Put the outputs into a directory for easy transfer
    base = './meshfem3D_files'
    if not os.path.exists(base):
        os.makedirs(base)

    print("\twriting Mesh_Par_file")
    with open(os.path.join(base, "Mesh_Par_file"), "w") as f:
        f.write(lines)


def write_interfaces(template, layers, interfaces, lat_min, lon_min, fids=[]):
    """
    Write the interfaces.dat file as well as the corresponding flat interface
    layers. Topo will need to be written manually

    :type template: str
    :param template: fid of the interfaces.dat template, preformatted
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
    # Put the outputs into a directory for easy transfer
    base = './meshfem3D_files'
    if not os.path.exists(base):
        os.makedirs(base)

    # Template for setting a flat layer
    flat_layer = f" .false. 2 2 {lon_min:.1f} {lat_min:.1f} 180.d0 180.d0\n"
    topo = ".false. 720 720 173.d0 -43.d0 0.00833d0 0.00833d0\n"

    # Write interfaces from bottom to top
    interfaces_reverse = interfaces.copy()
    interfaces_reverse.sort(reverse=True)

    # Write to a new file
    print("\twriting interfaces.dat")
    with open(os.path.join(base, "interfaces.dat"), "w") as f:
        f.write("# number of interfaces\n")
        f.write(" {ndoublings}\n".format(ndoublings=len(interfaces)))

        # Write flat layers up to topo
        for i, interface in enumerate(interfaces[:-1]):
            f.write("# interface number {}\n".format(i+1))
            f.write(flat_layer)

        # Write topo interface
        f.write("# interface number {} (topo)\n".format(i+2))
        f.write(topo)

        # Write number of elements in each layer
        f.write("# number of spectral elements per layer\n")
        for j, layer in enumerate(layers):
            f.write(f" {layer}\n")

    # Write the individual interface files
    for fid, interface in zip(fids, interfaces):
        # Skip top interface (topography)
        print(f"\twriting {fid}")
        with open(os.path.join(base, fid), "w") as f:
            for i in range(4):
                f.write("-{}\n".format(abs(int(interface * 1E3))))


def prepare_meshfem():
    """
    Run all the above functions and output the necessary meshfem3D files
    :return:
    """
    mesh_par_file_template = "./template_Mesh_Par_file"
    interfaces_template = "./template_interfaces.dat"

    print("Preparing Meshfem3D files")
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
    layers = vertical_doubling_proportions(depth=pars["mesh_depth_km"],
                                           ndoublings=pars["ndoublings"],
                                           interfaces=pars["interfaces"],
                                           grid_space=grid_space
                                           )
    write_mesh_par_file(template=mesh_par_file_template,
                        lat_min=pars["lat_min"], lat_max=pars["lat_max"],
                        lon_min=pars["lon_min"], lon_max=pars["lon_max"],
                        depth=pars["mesh_depth_km"],
                        utm_projection=pars["utm_projection"],
                        nex_x=nex_x, nex_y=nex_y, nproc_x=nproc_x,
                        nproc_y=nproc_y, ndoublings=pars["ndoublings"],
                        layers=layers,
                        mantle_from=pars["mantle_from_interface_idx"]
                        )
    write_interfaces(template=interfaces_template, layers=layers,
                     interfaces=pars["interfaces"], lat_min=pars["lat_min"],
                     lon_min=pars["lon_min"], fids=pars["interface_fids"]
                     )


def myround(x, base=5, choice='near'):
    """
    Round value x to nearest base, round 'up','down' or to 'near'est base

    :type x: float
    :param x: value to be rounded
    :type base: int
    :param base: nearest integer to be rounded to
    :type choice: str
    :param choice: method of rounding, 'up', 'down' or 'near'
    :rtype roundout: int
    :return: rounded value
    """
    if choice == 'near':
        roundout = int(base * round(float(x) / base))
    elif choice == 'down':
        roundout = int(base * np.floor(float(x) / base))
    elif choice == 'up':
        roundout = int(base * np.ceil(float(x) / base))

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
    from pyproj import Proj

    # Determine if the projection is north or south
    if utm_zone < 0:
        direction = "south"
    else:
        direction = "north"
    # Proj doesn't accept negative zones
    utm_zone = abs(utm_zone)

    # Proj requires a string to tell it how to convert the coordinates
    projstr = (f"+proj=utm +zone={utm_zone}, +{direction} +ellps=WGS84"
               " +datum=WGS84 +units=m +no_defs")

    # Initiate a Proj object and convert the coordinates
    my_proj = Proj(projstr)
    x_or_lon, y_or_lat = my_proj(lon_or_x, lat_or_y, inverse=inverse)

    return x_or_lon, y_or_lat


if __name__ == "__main__":
    prepare_meshfem()
