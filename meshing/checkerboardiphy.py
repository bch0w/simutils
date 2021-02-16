"""
Generate a checkerboard tomography model from an existing tomography .xyz file
For use with Specfem3D Cartesian.
"""
import os
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

plt.rcParams['image.cmap'] = 'seismic'


def xyz_reader(xyz_fid, save=True):
    """
    read .xyz file used as input for specfem and return information
    header parsing specifically related to specfem
    data is delineated as :
    x, y, z, vp[m/s], vs[m/s], rho[kg/m**3], Qp, Qs
    -origin and endpoints are the coordinates in meters
    -spacing given in units of meters for UTM projections
    -nx, ny, nz = [(end_x - orig_x)/spacing_x] + 1
    :type xyz_fid: str
    :param xyz_fid: file to read in
    :type save: bool
    :param save: to save the read files into numpy arrays
    :rtype header: dict
    :return header: a dictionary with relevant header informadtion
    :rtype data: np.ndarray
    :return data: all the data contained in the xyz file
    """
    if os.path.exists(xyz_fid + ".npy") and os.path.exists(xyz_fid + ".npz"):
        print("numpy data and header exist")
        header = np.load(xyz_fid + ".npz")
        data = np.load(xyz_fid + ".npy")
    elif os.path.exists(xyz_fid + ".npy"):
        print("numpy data exists, parsing header")
        data = np.load(xyz_fid + ".npy")
        header = parse_data_to_header(data)
        np.savez(xyz_fid, **header)  # save header in .npz
    else:
        print("no numpy files exist, reading xyz file")
        with open(xyz_fid) as f:
            lines = f.readlines()

            # parse the header and make sure the values are floats
            orig_x, orig_y, orig_z, end_x, end_y, end_z = \
                [float(_) for _ in lines[0].strip().split()]
            spacing_x, spacing_y, spacing_z = \
                [float(_) for _ in lines[1].strip().split()]
            nx, ny, nz = [float(_) for _ in lines[2].strip().split()]
            vp_min, vp_max, vs_min, vs_max, rho_min, rho_max = \
                [float(_) for _ in lines[3].strip().split()]

            header = {"orig_x": orig_x, "orig_y": orig_y, "orig_z": orig_z,
                      "end_x": end_x, "end_y": end_y, "end_z": end_z,
                      "spacing_x": spacing_x, "spacing_y": spacing_y,
                      "spacing_z": spacing_z, "nx": nx, "ny": ny, "nz": nz,
                      "vp_min": vp_min, "vp_max": vp_max, "vs_min": vs_min,
                      "vs_max": vs_max, "rho_min": rho_min, "rho_max": rho_max,
                      }

            data = np.array([_.strip().split() for _ in lines[4:]], dtype=float)
            if save:
                np.save(xyz_fid, data)  # save data in .npy
                np.savez(xyz_fid, **header)  # save header in .npz

    return header, data


def checkerboardiphy(xyz_fid, spacing_x, spacing_y=None, checker_z=None,
                     perturbation=0.02, mode="apply", apply_to=None,
                     taper_signal=None, no_incompletes=True, plot_fid=None,
                     **kwargs):
    """
    Read in the data from an XYZ tomography file and create a checkerboard
    overlay which has +/- {perturbation} checkers overlain on the tomography
    file. The checkers also include tapering so that there are no sharp
    transitions inside the checkerboard

    :type xyz_fid: str
    :param xyz_fid: path to file to read
    :type spacing_x: float
    :param spacing_x: spacing in meters for the x direction
    :type spacing_y: float
    :param spacing_x: spacing in meters for the y direction, optional,
        will default to spacing_x if not given
    :type checker_z: dict of floats
    :param checker_z: {"origin": starting depth for depth checkering in meters,
                       "spacing": spacing in meters of depth checkers}
    :type perturbation: float
    :param perturbation: perturbation to give to checkers
    :type mode: str
    :param mode: 'apply' - apply the given perturbation to the given velocity
        model to generate a perturbed 'target' velocity model
        'extract' - extract the given perturbations leaving the intermediate
            spaces as 0's to create a point spread function checkerboard
    :type apply_to: list of str
    :param apply_to: define which parameters to apply the checkerboard to.
        acceptable values: vp, vs, rho, qp, qs. By default, vp, vs.
    :type taper_signal: scipy.signal.function
    :param taper_signal: scipy signal to taper checkers by
    :type no_incompletes: bool  
    :param no_incompletes: if True, only create a checker if it doesnt cut off
    :type plot_fid: bool
    :param plot_fid: plot the overlay for confirmation
    :return:
    """
    # Pre-defined indices in the data array
    apply_dict = {"vp": 3, "vs": 4, "rho": 5, "qp": 6, "qs": 7}
    if not apply_to:
        apply_to = ["vp", "vs"]

    if not spacing_y:
        spacing_y = spacing_x

    # Read in the data
    header, data = xyz_reader(xyz_fid=xyz_fid, save=True)

    # Scipy signal kwargs
    if "std" in kwargs:
        std_x = kwargs["std"] / header["spacing_x"]
        std_y = kwargs["std"] / header["spacing_y"]

    # Initialize starting values
    checker_overlay = np.zeros(len(data))
    x = 1

    print(xyz_fid)
    # Loop through the x-axis, setting left and right boundaries
    for x_left in np.arange(header["orig_x"], header["end_x"], spacing_x):
        y = 1
        print(f"x_left: {x_left:.3E}/{header['end_x']:.3E}\t{x:+d}")
        x_right = x_left + spacing_x
        
        # Do not create incomplete checkers
        if no_incompletes and (x_right > header["end_x"]):
            continue

        # Create the tapered checker for a given checker defined by bounds
        x_checker = np.arange(x_left, x_right, header["spacing_x"])
        if "std" in kwargs:
            kwargs["std"] = std_x
        x_window = taper_signal(len(x_checker), **kwargs)

        # Loop through the y-axis, setting lower and upper boundaries
        for y_bot in np.arange(header["orig_y"], header["end_y"], spacing_y):
            print(f"\ty_bot: {y_bot:.3E}/{header['end_y']:.3E}\t{y:+d}")
            y_top = y_bot + spacing_y

            # Do not create incomplete checkers
            if no_incompletes and (y_top > header["end_y"]):
                continue

            # Create the tapered checker for the given checker defined by bounds
            y_checker = np.arange(y_bot, y_top, header["spacing_y"])
            if "std" in kwargs:
                kwargs["std"] = std_y
            y_window = taper_signal(len(y_checker), **kwargs)

            # Determine the sign of the overall checker, which is set by the
            # alternating sign of each inner checker axis
            checker = x * y

            # Determine where the data falls within this checker's bounds
            checker_indices = np.where(
                (data[:, 0] >= x_left) & (data[:, 0] < x_right) &
                (data[:, 1] >= y_bot) & (data[:, 1] < y_top)
            )
            # For each of the given indices, figure out the resulting overlay
            # Which is a multiplication of the overall sign, and the combination
            # of the x and y tapers within the checker
            for ind in checker_indices[0]:
                try:
                    checker_overlay[ind] = checker * (
                            x_window[np.where(x_checker == data[ind, 0])[0]] *
                            y_window[np.where(y_checker == data[ind, 1])[0]]
                    )
                except ValueError:


            # Flip the sign of the y-axis checker
            y *= -1
        # Flip the sign of the x-axis checker
        x *= -1

    # Checker in the vertical direction.
    # We pre-set the depth checkers based on "origin" and "spacing" to avoid any
    # conflicts with overlapping tomography files. Anything above "origin" will
    # be positive checker, this should be okay if topography << spacing
    if checker_z:
        print("\nDepth Layers")
        z = 1
        for z_top in np.arange(checker_z["origin"], header["orig_z"], 
                               -1 * checker_z["spacing"]):
            z_bottom = z_top - checker_z["spacing"]
            # Case for values above the origin, include w/ first checker
            if header["end_z"] > checker_z["origin"] and \
                                        z_top == checker_z["origin"]:
                z_top = header["end_z"]

            # Create the tapered checker for the vertical direction,
            # sample it fine so that we naturally interpolate to deal with the 
            # fact that we have chosen an arbitrary origin and so the sampling
            # points may not line up between our model and our checkers
            z_checker = np.arange(z_top, z_bottom, 
                                  max(-1000, -1 * header["spacing_z"])
                                  )
            if "std" in checker_z:
                kwargs["std"] = checker_z["std"] / header["spacing_z"]
            z_window = taper_signal(len(z_checker), **kwargs)

            print(f"\t{z_top:.3E} to {z_bottom:.3E}\t{z:+d}")
            
            # Determine where to apply the checker based on depth values
            checker_indices = np.where(
                (data[:, 2] <= z_top) & (data[:, 2] > z_bottom))
            for ind in checker_indices[0]:
                checker_overlay[ind] *= (z * z_window[np.where(
                                                z_checker == data[ind, 2])[0]])
            z *= -1

        # Case for values below the bottom, fill in the rest same as above
        if z_bottom > header["orig_z"]:
            z *= -1
            checker_indices = np.where(
                (data[:, 2] <= z_bottom) & (data[:, 2] > header["orig_z"]))
            print(f"z_top: {z_bottom:.3E} to {header['orig_z']:.3E}\t{z:+d}")
            for ind in checker_indices[0]:
                checker_overlay[ind] *= (z * z_window[np.where(
                                                z_checker == data[ind, 2])[0]])

    # Generate a quick plot to show a representation of the checkerboard
    if plot_fid is not None:
        print("Plot top down view")
        # Top down view only plots if surface is included
        z_ind = np.where(data[:, 2] ==  data[:, 2].max())[0]
        if z_ind.any():
            plt.scatter(x=data[z_ind, 0], y=data[z_ind, 1],
                        c=checker_overlay[z_ind], s=0.1)
            plt.xlabel("UTM-60 EASTING (m)")
            plt.ylabel("UTM-60 NORTHING (m)")
            plt.title(f"{xyz_fid}\n +/- {perturbation} perturbation\n"
                      f"{taper_signal.__name__} taper\n"
                      f"xy spacing; {spacing_x}m by {spacing_y}m")

            # plot coastline if possible
            coastline_fid = "./nz_resf_coast_mod_utm60H_xyz.npy"
            if os.path.exists(coastline_fid):
                coastline = np.load(coastline_fid)
                coastline = coastline[
                        np.where((coastline[:, 0] > header["orig_x"]) &
                                 (coastline[:, 0] < header["end_x"]) & 
                                 (coastline[:, 1] > header["orig_y"]) &
                                 (coastline[:, 1] < header["end_y"])
                                 )[0]]
                plt.scatter(coastline[:, 0], coastline[:, 1], c='k', marker='.',
                            s=1.)
            plt.gca().ticklabel_format(style='sci', axis='both')
            plt.gca().set_aspect(1)
            plt.savefig("{}.png".format(plot_fid))
            plt.close()
        
    # Side on view for checkers with depth
    if checker_z:
        print("Plot side-on view")
        # Determine where the maximum of the first set of checkers occurs
        checker_max = np.where(checker_overlay == checker_overlay.max())[0]
        y_checker_max = data[checker_max, 1].min()
        y_ind = np.where(data[:, 1] == y_checker_max)[0]

        # Plot the side view of a cut where the first maximum row occurs
        plt.scatter(x=data[y_ind, 0], y= data[y_ind, 2],
                    c=checker_overlay[y_ind], s=2)
        plt.xlabel("UTM-60 EASTING(m)")
        plt.ylabel("DEPTH (m)")
        z_space = checker_z["spacing"]
        plt.title(f"{xyz_fid}\n +/- {perturbation} perturbation\n"
                  f"{taper_signal.__name__} taper\n"
                  f"z spacing; {z_space}m")
        plt.gca().ticklabel_format(style='sci', axis='both')
        plt.savefig("{}_depth.png".format(plot_fid), dpi=100)
        plt.close()

    # Apply the peturbations, only to Vp and Vs values, not to rho or Q
    data_out = np.copy(data)
    checker_overlay *= perturbation
    for val in apply_to:
        i = apply_dict[val]
        if mode == "apply":
            offset = data_out[:, i]
        elif mode == "extract":
            offset = 0
        data_out[:, i] = offset + (checker_overlay * data_out[:, i])

    return checker_overlay, data_out


def parse_data_to_header(data):
    """
    Make a new header dictionary with adjusted data
    :param data:
    :return:
    """
    x_values = np.unique(data[:, 0])
    y_values = np.unique(data[:, 1])
    z_values = np.unique(data[:, 2])
    parsed_header = {
        "orig_x": x_values.min(), "orig_y": y_values.min(),
        "orig_z": z_values.min(), "end_x": x_values.max(),
        "end_y": y_values.max(), "end_z": z_values.max(),
        "spacing_x": abs(x_values[1] - x_values[0]),
        "spacing_y": abs(y_values[1] - y_values[0]),
        "spacing_z": abs(z_values[1] - z_values[0]),
        "nx": len(x_values), "ny": len(y_values), "nz": len(z_values),
        "vp_min": data[:, 4].min(), "vp_max": data[:, 3].max(),
        "vs_min": data[:, 4].min(), "vs_max": data[:, 4].max(),
        "rho_min": data[:, 5].min(), "rho_max": data[:, 5].max(),
    }
    return parsed_header


def write_xyz(header, data, fid_out):
    """
    Write out a new xyz file with proper header and data
    """
    with open(fid_out, "w") as f:
        # write header
        f.write("{:.1f} {:.1f} {:.1f} {:.1f} {:.1f} {:.1f}\n".format(
            header["orig_x"], header["orig_y"], header["orig_z"],
            header["end_x"], header["end_y"], header["end_z"])
        )
        f.write("{:.1f} {:.1f} {:.1f}\n".format(
            header["spacing_x"], header["spacing_y"], header["spacing_z"])
                )
        f.write("{:d} {:d} {:d}\n".format(
            header["nx"], header["ny"], header["nz"])
        )
        f.write("{:.1f} {:.1f} {:.1f} {:.1f} {:.1f} {:.1f}\n".format(
            header["vp_min"], header["vp_max"], header["vs_min"],
            header["vs_max"], header["rho_min"], header["rho_max"])
        )
        # write data by line
        for line in data:
            x, y, z, vp, vs, rho, qp, qs = line.tolist()
            f.write(
             "{:.1f} {:.1f} {:.1f} {:.1f} {:.1f} {:.1f} {:.1f} {:.1f}\n".format(
                   x, y, z, vp, vs, rho, qp, qs)
            )


if __name__ == "__main__":
    # =========================== Parameter set ================================
    fid_template="tomography_model_{}.xyz"
    taper_signal=signal.gaussian
    spacing_x=47E3
    spacing_y=48E3
    dict_z = {
        "shallow": {"origin": 3E3, "spacing": 5E3, "std": 500},
        "crust": {"origin": 7E3, "spacing": 5E3, "std": 500 },
        "mantle": {"origin": 47E3, "spacing": 5E3, "std": 500 }
    }
    perturbation=0.05
    std=5E3
    mode="extract"
    no_incompletes=False
    sections = ["mantle"]
    # =========================== Parameter set ================================

    # Create checkers with varying levels of perturbation
    for section in sections:
        print(f"section = {section}")
        fid = fid_template.format(section)

        # Save the outputs with a new tag
        tag = "_checker"
        fid_parts = os.path.splitext(fid)
        fid_out = f"{fid_parts[0]}{tag}{fid_parts[-1]}"

        # Set the Z-axis checkering based on section
        checker_z = dict_z[section]

        # Create the checkerboard data
        overlay, checkerboard_data = checkerboardiphy(
            xyz_fid=fid,  plot_fid=fid_out, spacing_x=spacing_x,
            spacing_y=spacing_y, checker_z=checker_z,
            perturbation=perturbation, taper_signal=taper_signal,
            no_incompletes=no_incompletes, mode=mode, std=std
        )
        checkerboard_header = parse_data_to_header(checkerboard_data)

        # Save the new data to a numpy array for easy i/o
        np.save(fid_out, checkerboard_data)

        # Write the new data to a new xyz file
        write_xyz(header=checkerboard_header, data=checkerboard_data,
                  fid_out=fid_out)