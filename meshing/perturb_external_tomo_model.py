"""
Perturb an existing .xyz tomography model (SPECFEM3D_Cartesian).
Originally called 'checkerboardiphy.py'. Utilities included to visualize the
perturbations.

Perturbation options:
    1. Checkerboard
    2. von Karman Spatial Correlation
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


def vonkarman(xyz_fid, **kwargs):
    """
    Read in the data from an XYZ tomography file and create a checkerboard
    overlay which has +/- {perturbation} checkers overlain on the tomography
    file. The checkers also include tapering so that there are no sharp
    transitions inside the checkerboard

    :type xyz_fid: str
    :param xyz_fid: path to file to reads
    """
    # Read in the data
    header, data = xyz_reader(xyz_fid=xyz_fid, save=True)

    # Initialize starting values
    checker_overlay = np.zeros(len(data))



def plot_top_down(data, header, fid, checkers, apply_to=None):
    """
    Plot a top down view of the returned data created by checkerboardiphy
    :return:
    """
    apply_dict = {"vp": 3, "vs": 4, "rho": 5, "qp": 6, "qs": 7}
    if not apply_to:
        apply_to = ["vp", "vs"]

    # Plot each depth layer for each changed parameter
    for val in apply_to:
        idx = apply_dict[val]

        # Use the checker data to figure out what x, y value slices through
        # the central maximum of the checkers
        checker_max = np.where(checkers == checkers.max())[0]
        z_max = data[checker_max, 2].max()
        z_idx = np.where(data[:, 2] == z_max)
        c = data[z_idx, idx]

        plt.scatter(x=data[z_idx, 0], y=data[z_idx, 1], c=c, s=0.1)
        plt.xlabel("UTM-60 EASTING (m)")
        plt.ylabel("UTM-60 NORTHING (m)")
        plt.title(f"{fid}, z={z_max:.0f}km, v={val}")

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
        plt.savefig(f"{fid}_{abs(z_max)}_{val}.png", dpi=80)
        plt.close()


def plot_side_on(data, fid, checkers, apply_to=None):
    """
    Plot a side on view of the returned data created by checkerboardiphy
    :return:
    """
    apply_dict = {"vp": 3, "vs": 4, "rho": 5, "qp": 6, "qs": 7}
    if not apply_to:
        apply_to = ["vp", "vs"]

    # Determine where the maximum of the first set of checkers occurs
    # Switching between x and y values
    for i in [0, 1]:
        name = ["x", "y"][i]
        for val in apply_to:
            idx = apply_dict[val]  # parameter to plot, e.g. 'vp'

            # Use the checker data to figure out what x, y value slices through
            # the central maximum of the checkers
            checker_max = np.where(checkers == checkers.max())[0]
            # For a fixed value on one axis we need a plane in the opposite axis
            # e.g. a fixed value of Y=0 requires a plot in XZ
            opp_idx = abs(i - 1)
            opp_axis = data[checker_max, opp_idx]
            oa_name = ["x", "y"][opp_idx]
            oa_val = opp_axis.min()
            oa_idx = np.where(data[:, opp_idx] == oa_val)[0]

            # Plot the side view of a cut where the first maximum row occurs
            plt.scatter(x=data[oa_idx, i], y=data[oa_idx, 2],
                        c=data[oa_idx, idx], s=2)

            plt.xlabel(f"UTM60S {name.upper()} (m)")
            plt.ylabel("DEPTH (m)")
            plt.title(f"{fid}, {oa_name.upper()}={oa_val}, {val}")
            plt.gca().ticklabel_format(style='sci', axis='both')
            plt.savefig(f"{fid}_{oa_name}{oa_val}_{val}_depth.png",
                        dpi=80)
            plt.close()


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
        "vp_min": data[:, 3].min(), "vp_max": data[:, 3].max(),
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
    perturb = vonkarman
    fids = ["tomography_model_shallow.xyz", "tomography_model_crust.xyz"]

    if perturb.__name__ == "checkerboard":
        kwargs = {
            "taper_signal": signal.windows.hann,
            "spacing_x": 40E3,
            "spacing_y": None,
            "dict_z": {
                "shallow": {"origin": 3E3, "spacing": 10E3, "std": None},
                "crust": {"origin": -7E3, "spacing": 10E3, "std": None},
                "mantle": {"origin": -40E3, "spacing": 24E3, "std": None}
            },
            "perturbation": 1,
            "std": None,
            "apply_to": ["vp", "vs"],
            "zero_values": [3000, 1500],
            "mode": "return",
            "stagger": True,
            "odd_checkering": False,
            "no_incompletes": {"x": False, "y": False, "z": True},
            "sections": ["shallow"],  # , "crust", "shallow"
            "invert_dict": {"mantle": True, "crust": True, "shallow": False},
            "plot": True
        }
    

    # =========================== Parameter set ================================

    # Create checkers with varying levels of perturbation
    for section in sections:
        print(f"section = {section}")
        fid = fid_template.format(section)

        # Save the outputs with a new tag so we dont overwrite original data
        tag = "_checker"
        fid_parts = os.path.splitext(fid)
        fid_out = f"{fid_parts[0]}{tag}{fid_parts[-1]}"

        # Set the Z-axis checkering based on section
        checker_z = dict_z[section]
        invert = invert_dict[section]

        overlay, checkerboard_data = perturb(
            xyz_fid=fid, spacing_x=spacing_x, spacing_y=spacing_y,
            checker_z=checker_z, zero_values=zero_values,
            apply_to=apply_to, perturbation=perturbation, invert=invert,
            taper_signal=taper_signal, no_incompletes=no_incompletes, mode=mode,
            stagger=stagger, std=std, odd_checkering=odd_checkering
        )
        checkerboard_header = parse_data_to_header(checkerboard_data)

        if plot:
            plot_top_down(checkerboard_data, checkerboard_header,
                          fid_out, overlay)
            plot_side_on(checkerboard_data, fid_out, overlay)

        # Save the new data to a numpy array for easy i/o
        np.save(fid_out, checkerboard_data)

        # Write the new data to a new xyz file
        write_xyz(header=checkerboard_header, data=checkerboard_data,
                  fid_out=fid_out)
