"""
Manipulate a SPECFEM3D external tomography file (.xyz) by creating a point
localized perturbation, for use in point spread function resolution analysis.
Allows for multiple, independent, perturbations to be included in a single run.
Treats each velocity model individually so points may overlap if the
perturbations are too close to the edges. Also treats each perturbation
individually so the User must ensure that they do not overlap on another
spatially.
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from checkerboardiphy import xyz_reader, parse_data_to_header, write_xyz


def perturb(data, origin, radii, sign, window, **kwargs):
    """
    Define a perturbation that can be re-used for all points in the model
    Origin must define actual discrete points within the model otherwise this
    will not work. Kwargs passed to the underlying scipy window signal

    :type data: np.ndarray
    :param data: multidimensional array defining the velocity model
        in the format [x, y, z, ...]
    :type origin: list of floats
    :param origin: [x0, y0, z0]
    :type radii: list of floats
    :param radii: [xr, yr, zr], radius of the given signal, the total length
        of the signal will then be 2r centered at the origin point
    :type signs: list of int
    :param signs: +1 or -1 to define if the peturbation is pos. or neg.
    :return:
    """
    print(f"\tperturbing with {window.__name__} window")
    prtrb = sign * np.ones(len(data))
    if "std" in kwargs:
        std = kwargs["std"]

    for i, (o, r) in enumerate(zip(origin, radii)):
        datum = data[:, i]
        unqdata = np.unique(datum)
        space = abs(unqdata[1] - unqdata[0])

        # Ensure that radius value is a multiple of the discretization
        if r < space:
            print(f"\tradius {r} < discretization {space}, rounding up")
            r = space
        if r % space != 0:
            print(f"\tradius {r} not a factor of discretization {space}", 
                  end=", ")
            r = int(space * np.ceil(float(r) / space))
            print(f"rounded up to {r}")

        # Temporary arrays to define the perturbation along a given axis
        range_ = np.arange(o - r, o + r, space)
        if "std" in kwargs:
            # gaussian std is a fraction of the radii for given dir
            kwargs["std"] = 0.4 * r / space
        prtrb_ = window(len(range_), **kwargs)

        # For each coordinate, we multiple all applicable values by the given
        # perturbation and set everything else to 0. This should allow each
        # coordinate to carve out its correct domain leading a 3D perturbation.
        for r_, p_ in zip(range_, prtrb_):
            prtrb[np.where(datum == r_)] *= p_

        # Set everything outside the range to 0 to isolate the perturbation
        prtrb[np.where((datum < range_.min()) | (datum > range_.max()))] *= 0

    return prtrb


def verify_location(data, origin, return_index=False):
    """
    Check if the given 3D coordinate is present in the given velocity model and
    determine the closest point based on the discretization of the model.

    :type data: np.ndarray
    :param data: multidimensional array defining the velocity model
        in the format [x, y, z, ...]
    :type origin: list of floats
    :param origin: [x0, y0, z0]
    :type return_index: bool
    :param return_index: return the index within the data array of the given
        coordinate point
    :rtype: tuple
    :return: (points, index), points defines 3D coordinate that best matches
        the velocity model discretization. index defines where that point lies
        within the data array, but is NoneType if return_index is False
    """
    points = []
    for i, val in enumerate(origin):
        arr = data[:, i]

        # Quick check to ensure the point falls completely within model bounds
        assert (arr.min() < val < arr.max()), "coord {i} OOB"

        # Find the nearest matching value which may differ due to discretization
        vals = np.unique(arr)
        if val not in vals:
            print(f"\torigin {val} is not a discrete point...", end=" ")
            val = vals[abs(vals - val).argmin()]
            print(f"{val} is nearest discrete point")

        points.append(val)

    # Determine the location of the given point within the matrix
    if return_index:
        index = np.where((data[:, 0] == points[0]) & (data[:, 1] == points[1]) &
                         (data[:, 2] == points[2]))
        return points, index
    else:
        return points, None


def plot_origin(data, origin, name=None, choice="vs"):
    """
    Make cross sections and depth slice through each axis at the origin point
    to get a quick glance at the size and shape of the input perturbation
    """
    # Convert string to index
    choice_idx = {"vp": 3, "vs": 4, "rho": 5, "qp": 6, "qs": 7}[choice]

    for i, o in enumerate(origin):
        data_unique = data[np.where(data[:, i] == o)]
        x_, y_, z_ = data_unique[:, :3].T
        c = data_unique[:, choice_idx]
        # Define the coordinate system based on origin point of interest
        if i == 0:
            x, y = y_, z_
        elif i == 1:
            x, y = x_, z_
        elif i == 2:
            x, y = x_, y_
    
        plt.scatter(x, y, c=c, s=1, cmap=plt.cm.seismic_r, vmin=-1, vmax=1)
        plt.title(f"{origin}[{i}] {choice}")
        if i != 2:
            plt.gca().set_aspect(2)
        else:
            plt.gca().set_aspect(1)
        plt.colorbar()
        plt.grid()
        plt.savefig(f"output/{name}_{i}_{choice}.png")
        plt.close()

    # if input("Acceptable?: "):
    #     return
    # else:
    #     sys.exit()

    
def main(anomalies, files, mode="return", apply_to=None, zero_values=None, 
         **kwargs):
    """
    Apply point local perturbations at given points to the given velocity model
    Kwargs passed to scipy signal function underlying perturb function
    """
    if not os.path.exists("output"):
        os.makedirs("output")
    apply_dict = {"vp": 3, "vs": 4, "rho": 5, "qp": 6, "qs": 7}

    # Default values determining which values to manipulate with perturbation
    if not apply_to:
        apply_to = ["vp", "vs"]
    if not zero_values:
        zero_values = [3000, 1500]

    # Quick assertion checks to make sure this won't crash at the end
    assert (len(apply_to) == len(zero_values)), \
        "len(apply_to) != len(zero_values)"
    assert (all([_ in apply_dict.keys() for _ in apply_to])), \
        f"unacceptable 'apply_to' values, see 'apply_dict'"

    apply_idx = [apply_dict[_] for _ in apply_to]

    for j, fid in enumerate(files):
        print(fid.split("/")[-1])
        _, data = xyz_reader(fid, save=True)
        data_out = data.copy()

        # Set output model to 0 so we can fill it up with perturbations
        data_out[:, apply_idx] *= 0.

        for i, (name, values) in enumerate(anomalies.items()):
            origin = values["origin"]
            radii = values["radii"]
            sign = values["sign"]

            print(f"PERTURBATION {name.upper()} ({i}/{len(anomalies) - 1})")
            try:
                origin, _ = verify_location(data, origin, return_index=False)
            except AssertionError:
                print(f"\t!!! origin {i} lies outside data bounds, skipping")
                continue

            # Return a perturbation array the same length as the data
            print(f"\tperturbing velocity model ({sign})")
            perturbation = perturb(data, origin, radii, sign, **kwargs)

            # Add each perturbation ontop of the perturbed, empty model
            for apply, zeroval in zip(apply_to, zero_values):
                idx = apply_dict[apply]
                if mode == "return":
                    prtrb = perturbation
                elif mode == "apply":
                    prtrb = perturbation * data[:, idx]
                elif mode == "extract":
                    prtrb = perturbation * data[:, idx]

                data_out[:, idx] += prtrb
        
            # Visualize perturbations
            plot_origin(data_out, origin, name, choice="vs")

        # Apply the offset values after all perturbations have been added
        for apply_, zeroval in zip(apply_to, zero_values):
            idx = apply_dict[apply_]
            if mode == "apply":
                offset = data[:, idx]
            else:
                offset = zeroval
            data_out[:, idx] += offset

        # Finalize by saving the data into new files for SPECFEM
        print("writing new .xyz file")
        header = parse_data_to_header(data_out)
        fid_ = os.path.basename(fid)
        fid_out = os.path.join("output", fid)

        # np.save(fid_out, data_out)  # I delete these anyways
        write_xyz(header, data_out, fid_out)

        # this is going to get written thrice but too lazy to fix
        if j == 0:
            print("writing config files")
            for i, (name, values) in enumerate(anomalies.items()):
                with open(f"output/{name}_cfg.txt", "w") as f:
                    for key, val in values.items():
                        f.write(f"{key}: {val}")


if __name__ == "__main__":
    # Each anomaly should have three keys: origin, radii, sign
    # 1) origin [x0, y0, z0]: point defining the center of the anomaly
    # 2) radii [x_r, y_r, z_r]: radius defining 3D extent of anomaly    
    # 3) sign (int): +1 or -1, the relative sign for fast or slow anomaly
    #
    # Note: Origins and radii must match the units and directions of the
    # underlying model. There is no checking involved to determine if correct

    anomalies = {
            # "mahia": 
            #     {"origin": [578000., 5668000., -12E3],
            #      "radii": [30E3, 30E3, 7E3],
            #      "sign": 1},
            # "porangahau": 
            #     {"origin": [466855., 5538861., -15E3],
            #      "radii": [15E3, 15E3, 5E3],
            #      "sign": 1},
            # "cook_strait": 
            #     {"origin": [307699., 5384284., -3E3],
            #      "radii": [20E3, 20E3, 7E3],
            #      "sign": -1},
            # "lg_cook": 
            #     {"origin": [317155., 5362128., -3E3],
            #      "radii": [40E3, 40E3, 10E3],
            #      "sign": -1},
            # "okataina": 
            #     {"origin": [440915., 5775292., -3E3],
            #      "radii": [15E3, 15E3, 7E3],
            #      "sign": -1},
            # "south_tvz": 
            #     {"origin": [418171.,5746769., 0],
            #      "radii": [20E3, 20E3, 15E3],
            #      "sign": -1},
            # "hg_palmy": 
            #     {"origin": [406512.,5563836., -10E3],
            #      "radii": [25E3, 25E3, 10E3],
            #      "sign": 1},
            # "mg_ruapehu": 
            #   {"origin": [378927.,5651838., -5E3],
            #    "radii": [60E3, 60E3, 10E3],
            #    "sign": 1},
            # "bottom": 
            #     {"origin": [460981.,5537260., -5E3],
            #      "radii": [25E3, 25E3, 10E3],
            #      "sign": 1},
            # "neg_pos": 
            #     {"origin": [427587.,5675437., -5E3],
            #      "radii": [20E3, 20E3, 12E3],
            #      "sign": 1},
            #"sm_taupo": 
            #    {"origin": [404105.,5702640., -4E3],
            #     "radii": [20E3, 20E3, 12E3],
            #     "sign": 1},
            # "taupo": 
            #     {"origin": [404105.,5702640., 0],
            #      "radii": [30E3, 30E3, 15E3],
            #      "sign": 1},
            "tarawera": 
                {"origin": [462757.,5791097.53, 0],
                 "radii": [25E3, 25E3, 15E3],
                 "sign": 1},
            # "tvz_deeper": 
            #     {"origin": [441230.,5776443., 0],
            #      "radii": [30E3, 30E3, 15E3],
            #      "sign": -1},
            # "tvz": 
            #     {"origin": [441230.,5776443., -3E3],
            #      "radii": [20E3, 20E3, 7E3],
            #      "sign": -1},
            # "whakamaru":
            #     {"origin": [427595., 5741709., 0.],
            #      "radii": [20E3, 20E3, 10E3],a
            #      "sign": -1},
            # "intraplate":
            #     {"origin": [509140., 5515069., -17.5E3],
            #      "radii": [30E3, 30E3, 10E3],
            #      "sign": -1},
            # "aboveintra":
            #     {"origin": [509140., 5515069., -5E3],
            #      "radii": [20E3, 20E3, 7E3],
            #      "sign": -1},
            # "taranaki":
            #     {"origin": [246758., 5646174., -3E3],
            #      "radii": [20E3, 20E3, 7E3],
            #      "sign": -1},
            
                }
    kwargs = {"window": signal.gaussian,
              "std": True
              }
    files = ["tomography_model_mantle.xyz",
             "tomography_model_crust.xyz",
             "tomography_model_shallow.xyz"
             ]
    main(anomalies, files, **kwargs)

