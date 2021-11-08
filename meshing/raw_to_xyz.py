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
    qp_fid = os.path.join(path, "Qpnzw2_zenodo/Qpnzw2p2xyzltln.tbl.txt")
    qs_fid = os.path.join(path, "Qsnzw2_zenodo/Qsnzw2p2xyzltln.tbl.txt")
    
    for path_ in [vel_fid, qp_fid, qs_fid]:
        assert(os.path.exists(path_))

    # Read in each file while skipping headers. Values are known from headers
    velocities = np.loadtxt(vel_fid, skiprows=2).transpose()
    vp, vpvs, vs, rho, sf_vp, sf_vpvs, x, y, z, lat, lon = velocities

    attenuation_p = np.loadtxt(qp_fid, skiprows=2).transpose()
    qp, sf_qp, x_qp, y_qp, z_qp, lat_qp, lon_qp = attenuation_p
    
    attenuation_s = np.loadtxt(qs_fid, skiprows=2).transpose()
    qs, sf_qs, x_qs, y_qs, z_qs, lat_qs, lon_qs = attenuation_s

    # Sanity check that vel and attenuation values are on the same grid
    assert(np.array_equal(x, x_qp)), "X values do not match"
    assert(np.array_equal(x, x_qs)), "X values do not match"
    assert(np.array_equal(y, y_qp)), "Y values do not match"
    assert(np.array_equal(y, y_qs)), "Y values do not match"
    assert(np.array_equal(lat, lat_qp)), "Lat values do not match"
    assert(np.array_equal(lat, lat_qs)), "Lat values do not match"
    assert(np.array_equal(lon, lon_qp)), "X values do not match"
    assert(np.array_equal(lon, lon_qs)), "X values do not match"
    
    # Carl wraps longitudes by 360, here we just ensure they are less than 180
    assert((abs(min(lon)) and max(lon)) <= 180), "Longitudes exceed +/-180"

    # Convert units and flip arrays so that Z layers are ordered bottom to top
    vp, vs, x, y, z = [_[::-1] * 1E3 for _ in [vp, vs, x, y, z]]
    qp, qs, rho, lat, lon = [_[::-1] for _ in [qp, qs, rho, lat, lon]]

    # Convert Z from BSL to positive depth    
    z *= -1

    return {"vp": vp, "vs": vs, "qp": qp, "qs":qs, "x": x, "y": y, "z": z,
            "rho": rho, "lat": lat, "lon": lon}


def check_qkappa(path_in):
    """
    Reads in velocities, attenuations and grid points from an external
    file, and generates .xyz files and headers for input to Specfem3D Cartesian

    :type path_in: str
    :param path_in: path to the input files, passed to the read function
    :type path_out: str
    :param path_out: path to save the output .xyz files
    """
    def qkappa(qp, qs, vp, vs):
        """
        Define Qkappa
        """ 
        f = 4/3 * (vs/vp) ** 2
        numer = qp * qs * (1 - f)
        denom = qs - f * qp
        return numer / denom
    
    def qkappa_2(qp, qs, vp, vs):
        """
        Check the calucaltions are the same for sanity
        """
        f = 4/3 * (vs/vp) ** 2
        numer = (1 - f) 
        denom = ((1 / qp) - (f / qs))
        return numer / denom

    def qkappa_reassign(qp, qs, vp, vs):
        """
        Equation 4 from Anderson and Hart 1978
        """
        f = 4/3 * (vs/vp) ** 2
        new_qp = []
        for i, (qp_, qs_, f_) in enumerate(zip(qp, qs, f)):
            qp_hold = qp_
            while qs_ / qp_ < f_:
                qp_ -= 1
                if qp_ <= 0:
                   break 
            # if qp_hold != qp_:
                # percent_change.append(((qp_hold - qp_) / qp_hold) * 100)
                # print(f"{qp_hold} -> {qp_}, "
                #       f"{((qp_hold-qp_)/qp_hold * 100):.2f}%")
            new_qp.append(qp_)

        return new_qp

    # Read in the data
    vel_dict = read_eberhart19(path_in)

    # Make sure the read function produces the correct outputs
    for key in ["vp", "vs", "qp", "qs", "x", "y", "z", "rho"]:
        assert(key in vel_dict), f"{key} not in dict"
        assert(not np.isnan(vel_dict[key]).all()), f"NaNs found in {key}"

    qp = vel_dict["qp"]
    qs = vel_dict["qs"]
    vp = vel_dict["vp"]
    vs = vel_dict["vs"]
    
    # Assert that the two formulations are identical
    qk = qkappa(qp=qp, qs=qs, vp=vp, vs=vs)
    qk2 = qkappa_2(qp=qp, qs=qs, vp=vp, vs=vs)
   
    # Rounding errors will mean the two formualations aren't exactly the same
    assert(np.allclose(qk, qk2)), "qkappa formulations inconsistent"

    qp_new = qkappa_reassign(qp=qp, qs=qs, vp=vp, vs=vs)
    qk_new = qkappa(qp=qp_new, qs=qs, vp=vp, vs=vs)
    
    
    import ipdb;ipdb.set_trace()  
    
if __name__ == "__main__":
    check_qkappa("./") 
