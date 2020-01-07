"""
Read in raw tomodd output files (provided by Donna Eberhart-Phillips),
and interpolate them onto a regular grid with the grid spacing, and boudnaries
of User choice. This script has been ported from the matlab scripts of 
Carl Tape (GEOTOOLS/read_tomo_eberhartnz.m, write_tomo_xyz_full.m)
"""
import os
import numpy as np


def read_eberhart19(path):
    """
    Read the tomo files provided by Donna. This is hardcoded to read the 
    .txt files that she provided via Zenodo. New tomo files will need adjusted
    read functions

    :type path: str
    :param path: path to the Zenodo .txt files
    :rtype: dict
    :return: dictionaries of the read in values, velocities in m/s, distances
        in m and ordered by Z (depth) from bottom to top
    """
    vel_fid = os.path.join(path, "nzw2.2_zenodo/vlnzw2p2dnxyzltln.tbl.txt")
    qp_fid = os.path.join(path, "Qpnzw2p2_zenodo/Qpnzw2p2xyzltln.tbl.txt")
    qs_fid = os.path.join(path, "Qsnzw2_zenodo/Qsnzw2p2xyzltln.tbl.txt")

    # Read in each file while skipping headers. Values are known from headers
    velocities = np.loadtxt(vel_fid, skiprows=2).transpose()
    vp, vpvs, vs, rho, sf_vp, sf_vpvs, x, y, z, lat, lon = velocities

    attenuation_p = np.loadtxt(qp_fid, skiprows=2).transpose()
    qp, sf_qp, x_qp, y_qp, z_qp, lat_qp, lon_qp = attenuation_p
    
    attenuation_s = np.loadtxt(qs_fid, skiprows=2).transpose()
    qs, sf_qs, x_qs, y_qs, z_qs, lat_qs, lon_qs = attenuation_s

    # Sanity check that vel and attenuation values are on the same grid
    assert(x == x_qp == x_qs), "X values do not match"
    assert(y == y_qp == y_qs), "Y values do not match"
    assert(z == z_qp == z_qs), "Z values do not match"
    assert(lat == lat_qp == lon_qs), "Lat values do not match"
    assert(lon == lat_qp == lon_qs), "Lon values do not match"
    
    # Carl wraps longitudes by 360, here we just ensure they are less than 180
    assert((abs(min(lon)) and max(lon) <= 180), "Longitudes exceed +/-180"

    # Convert units and flip arrays so that Z layers are ordered bottom to top
    vp, vs, x, y, z = [_[::-1] * 1E3 for _ in [vp, vs, x, y, z]]
    qp, qs, rho, lat, lon = [_[::-1] for _ in [qp, qs, rho, lat, lon]]

    # Convert Z from BSL to positive depth    
    z *= -1

    return {"vp": vp, "vs": vs, "qp": qp, "qs":qs, "x": x, "y": y, "z": z,
            "rho": rho, "lat": lat, "lon": lon}


def write_xyz(path_in='./', path_out='./', tag):
    """
    Reads in velocities, attenuations and grid points from an external
    file, and generates .xyz files and headers for input to Specfem3D Cartesian

    :type path_in: str
    :param path_in: path to the input files, passed to the read function
    :type path_out: str
    :param path_out: path to save the output .xyz files
    """
    vel_dict = read_eberhart19(path_in)

    # Make sure the read function produces the correct outputs
    for key in ["vp", "vs", "qp", "qs", "x", "y", "z", "rho"]:
        assert(key in val_dict), f"{key} not in dict"
        assert(np.isnan(val_dict[key]).all()), f"NaNs found in {key}"
    
    return
