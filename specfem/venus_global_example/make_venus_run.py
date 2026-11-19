#!/usr/bin/env python3
"""
Generate ready-to-run SPECFEM2D directories for the Venus global AXISYM
experiment (with and without atmosphere).

The script only writes files. It never runs the mesher or the solver and
never submits jobs. Copy the run folders to the cluster and submit
``run_venus.sh`` yourself.

Usage::

    python make_venus_run.py venus_config.py
    python make_venus_run.py --extract-table path/to/define_external_model.f90

Only the Python standard library is needed (Python >= 3.7).
"""
import argparse
import math
import os
import re
import runpy
import shutil
import stat
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

# Region flags written by the mesher (create_mesh_Venus.F90)
REG_UPPER_MANTLE = 1   # mesher: IREGION_MANTLE_BELOW_d670
REG_LOWER_MANTLE = 2   # mesher: IREGION_OUTER_CORE
REG_CORE = 3           # mesher: IREGION_INNER_CORE   (fluid)
REG_AIR = 4            # mesher: IREGION_MANTLE_CRUST_ABOVE_d670 (fluid)
FLUID_REGIONS = (REG_CORE, REG_AIR)

Q_NONE = 9999.0        # SPECFEM value meaning "no attenuation"
GLL_MIN_SPACING = 0.1727   # smallest GLL spacing / element width, NGLL=5
# Ratio of the solver's "Maximum suggested time step" to the simple
# estimate below, measured on a test run of this mesh (m=10, 18 Nov 2026)
DT_CALIBRATION = 1.44

DEFAULTS = dict(
    SPECFEM_BIN="/path/to/specfem2d/bin",
    FC="gfortran",
    FC_FLAGS="-O3",
    MESHER_TEMPLATE="templates/create_mesh_Venus.F90",
    PARFILE_TEMPLATE="templates/Par_file",
    SOURCE_TEMPLATE="templates/SOURCE",
    INTERIOR_TABLE="data/venus_interior.txt",
    AIR_PROFILE="data/venus_air_profile.txt",
    OUTPUT_DIR="runs",
    EXPERIMENT="venus",
    MIN_PERIOD_S=20.0,
    POINTS_PER_WAVELENGTH=5.0,
    MULTIPLICATION_FACTOR=None,
    INTERIOR_RADIAL_FACTOR=8,
    AIR_TOP_KM=100.0,
    AIR_TABLE_STEP_KM=1.0,
    SPONGE=True,
    SPONGE_BOTTOM_KM=70.0,
    SPONGE_Q_MIN=3.0,
    CASES={"air": 1.0, "noair": 1.0e-3},
    SOLID_ATTENUATION=False,
    SOURCE_TYPE="force",
    SOURCE_DEPTH_KM=10.0,
    SOURCE_FACTOR=1.0e10,
    SIM_DURATION_S=2000.0,
    DT_OVERRIDE=None,
    DT_SAFETY=0.8,
    SEISMOTYPE="2,4",
    GROUND_STATION_DEPTH_M=20.0,
    GROUND_STATION_DEG=[5, 10, 20, 30, 45, 60, 90],
    AIR_STATION_DEG=[10, 30, 60],
    AIR_STATION_ALT_KM=[1, 10, 30, 50, 60],
    NPROC=168,
    SLURM={},
    US_PER_ELEMENT_STEP=(1.0, 3.0),
)


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------
def die(msg):
    """
    Print an error message and stop.

    :type msg: str
    :param msg: message to print
    """
    sys.exit(f"ERROR: {msg}")


def path_from_config(cfg, key):
    """
    Resolve a path from the config relative to the config file location.

    :type cfg: dict
    :param cfg: configuration dictionary
    :type key: str
    :param key: config key holding a path
    :rtype: str
    :return: absolute path
    """
    value = cfg[key]
    if os.path.isabs(value):
        return value
    return os.path.join(cfg["_CONFIG_DIR"], value)


def load_config(config_path):
    """
    Load a Python config file and merge it with the defaults.

    :type config_path: str
    :param config_path: path to the config file (e.g., venus_config.py)
    :rtype: dict
    :return: configuration with defaults filled in
    """
    user = runpy.run_path(config_path)
    cfg = dict(DEFAULTS)
    for key, value in user.items():
        if key.isupper():
            cfg[key] = value
    unknown = [k for k in user if k.isupper() and k not in DEFAULTS]
    if unknown:
        die(f"unknown setting(s) in config: {', '.join(unknown)}")
    cfg["_CONFIG_DIR"] = os.path.dirname(os.path.abspath(config_path))
    return cfg


# ---------------------------------------------------------------------------
# model tables
# ---------------------------------------------------------------------------
def read_fortran_arrays(f90_path):
    """
    Read hard-coded model arrays from an old define_external_model.f90.

    :type f90_path: str
    :param f90_path: path to the old define_external_model.f90
    :rtype: dict
    :return: arrays keyed by name, each a list ordered by index
    """
    text = open(f90_path).read()
    names = ["radius", "density", "vp", "vs", "Qkappa", "Qmu"]
    arrays = {}
    for name in names:
        pattern = name + r"_ak135\(\s*(\d+)\)\s*=\s*([-+\dEe.]+)"
        vals = {int(i): float(v) for i, v in re.findall(pattern, text)}
        if not vals:
            die(f"no {name}_ak135 entries found in {f90_path}")
        arrays[name] = [vals[i] for i in sorted(vals)]
    nrow = len(arrays["radius"])
    if any(len(a) != nrow for a in arrays.values()):
        die("model arrays in the Fortran file have different lengths")
    return arrays


def extract_interior(f90_path, out_path):
    """
    Convert the solid-planet part of the Fortran table to a region-tagged
    SI text table (radius in m, rho in kg/m3, velocities in m/s).

    Rows are split at the two discontinuities that the mesher honors
    (core-mantle boundary and ground surface, both stored as duplicated
    radii in the table) and at the mid-mantle mesh boundary.

    :type f90_path: str
    :param f90_path: path to the old define_external_model.f90
    :type out_path: str
    :param out_path: where to write the interior table
    """
    a = read_fortran_arrays(f90_path)
    rad = a["radius"]
    dups = [i for i in range(1, len(rad)) if rad[i] == rad[i - 1]]
    if len(dups) < 2:
        die("expected duplicated radii at the CMB and at the surface")
    i_cmb, i_surf = dups[0], dups[-1]       # index of the upper duplicate
    r_cmb, r_surf = rad[i_cmb], rad[i_surf]
    # mid-mantle mesh boundary: the table radius closest to 5313.57 km
    i_mid = min(range(i_cmb, i_surf), key=lambda i: abs(rad[i] - 5313568.))
    r_mid = rad[i_mid]

    def row(i, region):
        """Build one SI row for table index i."""
        return (region, rad[i], a["density"][i] * 1e3, a["vp"][i] * 1e3,
                a["vs"][i] * 1e3, a["Qkappa"][i], a["Qmu"][i])

    rows = [row(i, REG_CORE) for i in range(0, i_cmb)]
    rows += [row(i, REG_LOWER_MANTLE) for i in range(i_cmb, i_mid + 1)]
    rows += [row(i, REG_UPPER_MANTLE) for i in range(i_mid, i_surf)]
    with open(out_path, "w") as f:
        f.write(f"# Venus interior extracted from {os.path.basename(f90_path)}"
                "\n# Region flags: 1 upper mantle, 2 lower mantle, 3 core\n"
                f"# CMB radius {r_cmb:.1f} m, mid-mantle boundary "
                f"{r_mid:.1f} m, surface {r_surf:.1f} m\n"
                "# region radius_m rho_kg_m3 vp_m_s vs_m_s qkappa qmu\n")
        for r in rows:
            f.write("{:d} {:.3f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}\n"
                    .format(*r))
    print(f"wrote {out_path}: {len(rows)} rows, CMB {r_cmb:.0f} m, "
          f"mid-mantle {r_mid:.0f} m, surface {r_surf:.0f} m")


def read_interior(path):
    """
    Read the region-tagged interior table.

    :type path: str
    :param path: interior table written by extract_interior
    :rtype: list of tuple
    :return: rows (region, r, rho, vp, vs, qk, qm)
    """
    rows = []
    for line in open(path):
        if line.strip() and not line.lstrip().startswith("#"):
            p = line.split()
            rows.append((int(p[0]),) + tuple(float(v) for v in p[1:7]))
    if not rows:
        die(f"interior table {path} is empty")
    return rows


def region_bounds(rows):
    """
    Get the radius range of each region in a table.

    :type rows: list of tuple
    :param rows: model rows
    :rtype: dict
    :return: {region: (r_min, r_max)}
    """
    out = {}
    for row in rows:
        lo, hi = out.get(row[0], (row[1], row[1]))
        out[row[0]] = (min(lo, row[1]), max(hi, row[1]))
    return out


def read_air_profile(path):
    """
    Read the atmosphere profile (altitude, density, sound speed).

    :type path: str
    :param path: air profile file (columns alt_km T P rho c)
    :rtype: list of tuple
    :return: rows (alt_m, rho, c) sorted by altitude
    """
    prof = []
    for line in open(path):
        if line.strip() and not line.lstrip().startswith("#"):
            p = [float(v) for v in line.split()]
            prof.append((p[0] * 1e3, p[3], p[4]))
    prof.sort()
    if len(prof) < 2:
        die("air profile needs at least two rows")
    return prof


def interp_air(prof, z):
    """
    Interpolate the air profile: log-linear density, linear sound speed.

    :type prof: list of tuple
    :param prof: rows (alt_m, rho, c)
    :type z: float
    :param z: altitude in m
    :rtype: tuple
    :return: (rho, c)
    """
    if z <= prof[0][0]:
        return prof[0][1], prof[0][2]
    for (z0, r0, c0), (z1, r1, c1) in zip(prof, prof[1:]):
        if z0 <= z <= z1:
            f = (z - z0) / (z1 - z0)
            rho = math.exp(math.log(r0) + f * (math.log(r1) - math.log(r0)))
            return rho, c0 + f * (c1 - c0)
    die(f"air profile does not reach {z / 1e3:.1f} km; extend the file")


def air_rows(cfg, r_surf, density_scale):
    """
    Build table rows for the atmosphere, including the Q sponge.

    :type cfg: dict
    :param cfg: configuration
    :type r_surf: float
    :param r_surf: ground radius in m
    :type density_scale: float
    :param density_scale: factor on air density (1e-3 = "no atmosphere")
    :rtype: list of tuple
    :return: rows (region, r, rho, vp, vs, qk, qm)
    """
    prof = read_air_profile(path_from_config(cfg, "AIR_PROFILE"))
    top = cfg["AIR_TOP_KM"] * 1e3
    step = cfg["AIR_TABLE_STEP_KM"] * 1e3
    sb = cfg["SPONGE_BOTTOM_KM"] * 1e3
    n = int(math.ceil(top / step))
    rows = []
    for k in range(n + 1):
        z = min(k * step, top)
        rho, c = interp_air(prof, z)
        qk = Q_NONE
        if cfg["SPONGE"] and z > sb:
            inv_q = (z - sb) / (top - sb) / cfg["SPONGE_Q_MIN"]
            qk = min(Q_NONE, 1.0 / max(inv_q, 1.0e-4))
        rows.append((REG_AIR, r_surf + z, rho * density_scale, c, 0.0,
                     qk, Q_NONE))
    return rows


def write_model_table(rows, path, header):
    """
    Write DATA/venus_model.txt in the format the Fortran routine reads.

    :type rows: list of tuple
    :param rows: model rows grouped by region, sorted by radius
    :type path: str
    :param path: output path
    :type header: str
    :param header: comment text describing the model
    """
    with open(path, "w") as f:
        for line in header.splitlines():
            f.write(f"# {line}\n")
        f.write("# region radius_m rho_kg_m3 vp_m_s vs_m_s qkappa qmu\n")
        f.write(f"{len(rows)}\n")
        for r in rows:
            f.write("{:d} {:.3f} {:.6e} {:.4f} {:.4f} {:.4f} {:.4f}\n"
                    .format(*r))


# ---------------------------------------------------------------------------
# mesh parameters, time step, cost
# ---------------------------------------------------------------------------
def even_up(n):
    """
    Round up to the next even integer.

    :type n: float
    :param n: value to round
    :rtype: int
    :return: even integer >= n
    """
    n = int(math.ceil(n))
    return n + (n % 2)


def mesh_parameters(cfg, r_cmb, r_mid, r_surf, c_min):
    """
    Choose the mesher parameters for the requested minimum period.

    :type cfg: dict
    :param cfg: configuration
    :type r_cmb: float
    :param r_cmb: core-mantle boundary radius [m] (mesher RICB)
    :type r_mid: float
    :param r_mid: mid-mantle boundary radius [m] (mesher RCMB)
    :type r_surf: float
    :param r_surf: ground radius [m] (mesher R670)
    :type c_min: float
    :param c_min: slowest sound speed in the air [m/s]
    :rtype: dict
    :return: mesher parameters and resolution numbers
    """
    r_top = r_surf + cfg["AIR_TOP_KM"] * 1e3
    h_max = 4.0 / cfg["POINTS_PER_WAVELENGTH"] * c_min * cfg["MIN_PERIOD_S"]
    m = cfg["MULTIPLICATION_FACTOR"]
    if m is None:
        m = even_up(2 * math.pi * r_top / (240.0 * h_max))
    h_lat = 2 * math.pi * r_top / (240.0 * m)
    n_air = max(4, int(math.ceil(cfg["AIR_TOP_KM"] * 1e3 / (0.5 * h_lat))))
    k = cfg["INTERIOR_RADIAL_FACTOR"]
    p = dict(
        multiplication_factor=m,
        nspec_surf_whole_circle=480 * m,
        nspec_rad_670_surf=2 * n_air,
        nspec_rad_CMB_670=k * m,
        nspec_rad_doubling_OC_to_CMB=k * m,
        nspec_rad_ICB_to_doubling_OC=k * m,
        nspec_rad_Cube_ICB=16 * m,
        R_EARTH=r_top, R670=r_surf, RCMB=r_mid, RICB=r_cmb,
        h_air_lateral=h_lat, h_air_radial=cfg["AIR_TOP_KM"] * 1e3 / n_air,
        h_max=h_max,
    )
    # same checks the Fortran mesher makes, so it never stops on the cluster
    checks = [
        (p["nspec_rad_670_surf"] % 2 == 0, "nspec_rad_670_surf must be even"),
        (p["nspec_rad_CMB_670"] % 4 == 0, "nspec_rad_CMB_670 % 4 != 0"),
        (p["nspec_rad_ICB_to_doubling_OC"] % 4 == 0,
         "nspec_rad_ICB_to_doubling_OC % 4 != 0"),
        (p["nspec_rad_doubling_OC_to_CMB"] % 4 == 0,
         "nspec_rad_doubling_OC_to_CMB % 4 != 0"),
        (p["nspec_rad_Cube_ICB"] % 4 == 0, "nspec_rad_Cube_ICB % 4 != 0"),
        (p["nspec_surf_whole_circle"] % 16 == 0,
         "nspec_surf_whole_circle % 16 != 0"),
        (p["nspec_rad_CMB_670"] >= 8, "nspec_rad_CMB_670 too small"),
        (p["nspec_rad_ICB_to_doubling_OC"] >= 8,
         "nspec_rad_ICB_to_doubling_OC too small"),
    ]
    for ok, msg in checks:
        if not ok:
            die(f"mesh parameter check failed: {msg} (change "
                "MULTIPLICATION_FACTOR or INTERIOR_RADIAL_FACTOR)")
    if h_lat > 1.001 * h_max:
        print(f"WARNING: air elements are {h_lat / 1e3:.2f} km wide but "
              f"{h_max / 1e3:.2f} km is needed for {cfg['MIN_PERIOD_S']} s")
    p["nspec"] = count_elements(p)
    return p


def count_elements(p):
    """
    Number of elements the Fortran mesher will create (half mesh).

    :type p: dict
    :param p: mesher parameters
    :rtype: int
    :return: element count
    """
    s = p["nspec_surf_whole_circle"]
    a, b = p["nspec_rad_670_surf"], p["nspec_rad_CMB_670"]
    c, d = p["nspec_rad_doubling_OC_to_CMB"], p["nspec_rad_ICB_to_doubling_OC"]
    e = p["nspec_rad_Cube_ICB"]
    total = ((s // 2) * (a // 2) + (s // 4) * (b // 4 - 1) + (s // 4) * 3
             + (s // 4) * (c // 4) + (s // 8) * (d // 4 - 1) + (s // 8) * 3
             + (s // 8) * (e // 4) + (s // 32) ** 2)
    return total // 2


def estimate_dt(p, rows):
    """
    Rough stable time step from the thinnest, fastest elements.

    The solver prints the real limit; use DT_OVERRIDE to adopt it.

    :type p: dict
    :param p: mesher parameters
    :type rows: list of tuple
    :param rows: full model table
    :rtype: float
    :return: time step estimate [s]
    """
    vp_max = {}
    for row in rows:
        vp_max[row[0]] = max(vp_max.get(row[0], 0.0), row[3])
    m = p["multiplication_factor"]
    r_dbl = p["RICB"] + 0.56 * (p["RCMB"] - p["RICB"])
    widths = [
        ((p["R670"] - p["RCMB"]) / (p["nspec_rad_CMB_670"] / 2),
         vp_max[REG_UPPER_MANTLE]),
        ((r_dbl - p["RICB"]) / (p["nspec_rad_ICB_to_doubling_OC"] / 2),
         vp_max[REG_LOWER_MANTLE]),
        ((p["RCMB"] - r_dbl) / (p["nspec_rad_doubling_OC_to_CMB"] / 2),
         vp_max[REG_LOWER_MANTLE]),
        (0.7 * 1300e3 / (15 * m), vp_max[REG_CORE]),      # central cube
        (p["h_air_radial"], vp_max[REG_AIR]),
    ]
    return DT_CALIBRATION * min(0.5 * GLL_MIN_SPACING * h / v
                                for h, v in widths)


# ---------------------------------------------------------------------------
# file writers
# ---------------------------------------------------------------------------
def replace_once(text, pattern, repl, what):
    """
    Replace exactly one regex match, or stop with a clear message.

    :type text: str
    :param text: text to edit
    :type pattern: str
    :param pattern: regex with one group for the part to keep
    :type repl: str
    :param repl: replacement value
    :type what: str
    :param what: name used in the error message
    :rtype: str
    :return: edited text
    """
    new, n = re.subn(pattern, lambda mt: mt.group(1) + repl, text,
                     flags=re.MULTILINE)
    if n != 1:
        die(f"expected one '{what}' line in the template, found {n}")
    return new


def write_mesher(cfg, p, out_path):
    """
    Write the mesher source with the chosen parameters patched in.

    :type cfg: dict
    :param cfg: configuration
    :type p: dict
    :param p: mesher parameters
    :type out_path: str
    :param out_path: output .F90 path
    """
    text = open(path_from_config(cfg, "MESHER_TEMPLATE")).read()
    ints = ["multiplication_factor", "nspec_surf_whole_circle",
            "nspec_rad_670_surf", "nspec_rad_CMB_670",
            "nspec_rad_doubling_OC_to_CMB", "nspec_rad_ICB_to_doubling_OC",
            "nspec_rad_Cube_ICB"]
    for name in ints:
        text = replace_once(
            text, r"^(\s*integer,\s*parameter\s*::\s*" + name + r"\s*=\s*)"
            r"[^!\n]*", f"{p[name]} ", name)
    for name in ["R_EARTH", "R670", "RCMB", "RICB"]:
        text = replace_once(
            text, r"^(\s*double precision,\s*parameter\s*::\s*" + name +
            r"\s*=\s*)[^!\n]*", f"{p[name]:.3f}d0 ", name)
    text = replace_once(
        text, r"^(\s*logical,\s*parameter\s*::\s*output_gnuplot_grid\s*=\s*)"
        r"[^!\n]*", ".false. ", "output_gnuplot_grid")
    with open(out_path, "w") as f:
        f.write(text)


def edit_keys(template, updates, out_path, label):
    """
    Copy a SPECFEM parameter file, replacing the values of given keys.

    :type template: str
    :param template: path to the template file
    :type updates: dict
    :param updates: {key: new value as string}
    :type out_path: str
    :param out_path: where to write the edited file
    :type label: str
    :param label: file name used in error messages
    """
    lines = open(template).read().splitlines(True)
    found = set()
    for i, line in enumerate(lines):
        mt = re.match(r"^(\s*)([A-Za-z_0-9]+)(\s*=\s*)([^#\n]*?)(\s*(#.*)?)$",
                      line.rstrip("\n"))
        if mt and mt.group(2) in updates:
            key = mt.group(2)
            if key in found:
                die(f"key {key} appears twice in the {label} template")
            found.add(key)
            lines[i] = (f"{mt.group(1)}{key}{mt.group(3)}{updates[key]}"
                        f"{mt.group(5) or ''}\n")
    missing = set(updates) - found
    if missing:
        die(f"keys missing from the {label} template: {sorted(missing)}")
    with open(out_path, "w") as f:
        f.writelines(lines)


def fbool(flag):
    """
    Format a Python bool as a Fortran logical.

    :type flag: bool
    :param flag: value
    :rtype: str
    :return: '.true.' or '.false.'
    """
    return ".true." if flag else ".false."


def write_stations(cfg, r_surf, out_path):
    """
    Write DATA/STATIONS: ground stations just below the surface and air
    stations at fixed altitudes. Angles are measured from the source
    (north pole, +z axis) toward the equator.

    :type cfg: dict
    :param cfg: configuration
    :type r_surf: float
    :param r_surf: ground radius [m]
    :type out_path: str
    :param out_path: output path
    :rtype: int
    :return: number of stations
    """
    sb = cfg["SPONGE_BOTTOM_KM"]
    n = 0
    with open(out_path, "w") as f:
        for deg in cfg["GROUND_STATION_DEG"]:
            r = r_surf - cfg["GROUND_STATION_DEPTH_M"]
            th = math.radians(deg)
            f.write(f"G{deg:03.0f}  VE  {r * math.sin(th):.2f}  "
                    f"{r * math.cos(th):.2f}  0.0  0.0\n")
            n += 1
        for deg in cfg["AIR_STATION_DEG"]:
            for alt in cfg["AIR_STATION_ALT_KM"]:
                if cfg["SPONGE"] and alt >= sb:
                    print(f"WARNING: air station at {alt} km is inside the "
                          "sponge; its amplitudes are not physical")
                r = r_surf + alt * 1e3
                th = math.radians(deg)
                f.write(f"A{deg:03.0f}H{alt:03.0f}  VE  "
                        f"{r * math.sin(th):.2f}  {r * math.cos(th):.2f}  "
                        "0.0  0.0\n")
                n += 1
    return n


def write_nummaterial(path):
    """
    Write the material file that tells SPECFEM which regions are fluid.
    Values other than the domain flag are overwritten by the external
    model, but vs must be zero for the fluid regions.

    :type path: str
    :param path: output path
    """
    with open(path, "w") as f:
        f.write("# domain_id material_id rho vp vs Qk Qmu aniso\n"
                "# domain 1 = acoustic (fluid), 2 = elastic (solid)\n"
                "2 1 3300.d0 8000.d0 4500.d0 9999.d0 9999.d0 0\n"
                "2 2 5000.d0 12000.d0 6500.d0 9999.d0 9999.d0 0\n"
                "1 3 11000.d0 9000.d0 0.d0 9999.d0 9999.d0 0\n"
                "1 4 65.d0 400.d0 0.d0 9999.d0 9999.d0 0\n")


def write_run_script(cfg, case_dir):
    """
    Write a job script that builds the mesh and runs the solver. It is not
    submitted; edit the #SBATCH lines if needed and submit it yourself.

    :type cfg: dict
    :param cfg: configuration
    :type case_dir: str
    :param case_dir: run directory
    """
    sl = dict(cfg["SLURM"])
    head = ["#!/bin/bash"]
    for key, value in sl.items():
        head.append(f"#SBATCH --{key.replace('_', '-')}={value}")
    body = f"""
set -e
cd "$(dirname "$0")"
mkdir -p OUTPUT_FILES MESH bin
rm -rf OUTPUT_FILES/*
ln -sf {cfg['SPECFEM_BIN']}/xmeshfem2D bin/xmeshfem2D
ln -sf {cfg['SPECFEM_BIN']}/xspecfem2D bin/xspecfem2D

echo "building the mesh: $(date)"
{cfg['FC']} {cfg['FC_FLAGS']} -o bin/xcreate_mesh create_mesh_venus.F90
./bin/xcreate_mesh > OUTPUT_FILES/output_create_mesh.txt

cp DATA/Par_file DATA/SOURCE DATA/STATIONS DATA/venus_model.txt \\
   OUTPUT_FILES/
echo "partitioning: $(date)"
mpirun -np {cfg['NPROC']} ./bin/xmeshfem2D
echo "solver: $(date)"
mpirun -np {cfg['NPROC']} ./bin/xspecfem2D
echo "done: $(date)"
"""
    path = os.path.join(case_dir, "run_venus.sh")
    with open(path, "w") as f:
        f.write("\n".join(head) + "\n" + body)
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR)


# ---------------------------------------------------------------------------
# main driver
# ---------------------------------------------------------------------------
def make_runs(cfg):
    """
    Generate one run directory per case in cfg['CASES'].

    :type cfg: dict
    :param cfg: configuration
    """
    interior = read_interior(path_from_config(cfg, "INTERIOR_TABLE"))
    if not cfg["SOLID_ATTENUATION"]:
        interior = [r[:5] + (Q_NONE, Q_NONE) for r in interior]
    bounds = region_bounds(interior)
    r_cmb = bounds[REG_CORE][1]
    r_mid = bounds[REG_LOWER_MANTLE][1]
    r_surf = bounds[REG_UPPER_MANTLE][1]
    if cfg["SPONGE"] and not 0 < cfg["SPONGE_BOTTOM_KM"] < cfg["AIR_TOP_KM"]:
        die("SPONGE_BOTTOM_KM must be between 0 and AIR_TOP_KM")

    base_air = air_rows(cfg, r_surf, 1.0)
    c_min = min(r[3] for r in base_air)
    tau_air = sum((base_air[i + 1][1] - base_air[i][1]) /
                  (0.5 * (base_air[i][3] + base_air[i + 1][3]))
                  for i in range(len(base_air) - 1))
    p = mesh_parameters(cfg, r_cmb, r_mid, r_surf, c_min)

    f_max = 1.0 / cfg["MIN_PERIOD_S"]
    f0 = f_max / 2.5
    dt_est = estimate_dt(p, interior + base_air)
    dt = cfg["DT_OVERRIDE"] or float(f"{cfg['DT_SAFETY'] * dt_est:.2g}")
    nstep = int(math.ceil(cfg["SIM_DURATION_S"] / dt))
    nstep += (-nstep) % 100
    lo, hi = (p["nspec"] * nstep * us * 1e-6 / 3600 * 1.15
              for us in cfg["US_PER_ELEMENT_STEP"])

    exp_dir = os.path.join(path_from_config(cfg, "OUTPUT_DIR"),
                           cfg["EXPERIMENT"])
    for case, scale in cfg["CASES"].items():
        d = os.path.join(exp_dir, case)
        for sub in ("DATA", "MESH"):
            os.makedirs(os.path.join(d, sub), exist_ok=True)

        rows = interior + air_rows(cfg, r_surf, scale)
        write_model_table(
            rows, os.path.join(d, "DATA", "venus_model.txt"),
            f"Venus model for case '{case}' (air density x {scale:g})\n"
            f"sponge: {cfg['SPONGE']}, from {cfg['SPONGE_BOTTOM_KM']} km "
            f"to {cfg['AIR_TOP_KM']} km, Q_min {cfg['SPONGE_Q_MIN']}")
        write_mesher(cfg, p, os.path.join(d, "create_mesh_venus.F90"))
        write_nummaterial(os.path.join(
            d, "MESH", "Nummaterial_velocity_file_NO_MUD"))
        nsta = write_stations(cfg, r_surf, os.path.join(d, "DATA",
                                                         "STATIONS"))
        edit_keys(path_from_config(cfg, "PARFILE_TEMPLATE"), {
            "title": f"Venus {cfg['EXPERIMENT']} {case}",
            "NPROC": str(cfg["NPROC"]),
            "NSTEP": str(nstep),
            "DT": f"{dt:g}",
            "time_stepping_scheme": "1",
            "P_SV": ".true.",
            "AXISYM": ".true.",
            "MODEL": "external",
            "ATTENUATION_VISCOELASTIC": fbool(cfg["SOLID_ATTENUATION"]),
            "ATTENUATION_VISCOACOUSTIC": fbool(cfg["SPONGE"]),
            "N_SLS": "3",
            "PML_BOUNDARY_CONDITIONS": ".false.",
            "STACEY_ABSORBING_CONDITIONS": ".false.",
            "read_external_mesh": ".true.",
            "use_existing_STATIONS": ".true.",
            "rec_normal_to_surface": ".false.",
            "seismotype": cfg["SEISMOTYPE"],
            "NTSTEP_BETWEEN_OUTPUT_SEISMOS": str(nstep // 10),
            "NTSTEP_BETWEEN_OUTPUT_INFO": str(max(100, nstep // 50)),
            "output_color_image": ".false.",
            "NTSTEP_BETWEEN_OUTPUT_IMAGES": str(nstep),
        }, os.path.join(d, "DATA", "Par_file"), "Par_file")
        src = {"force": ("1", "0."), "explosion": ("2", "0.")}
        if cfg["SOURCE_TYPE"] not in src:
            die("SOURCE_TYPE must be 'force' or 'explosion'")
        stype, angle = src[cfg["SOURCE_TYPE"]]
        edit_keys(path_from_config(cfg, "SOURCE_TEMPLATE"), {
            "source_surf": ".false.",
            "xs": "0.",
            "zs": f"{r_surf - cfg['SOURCE_DEPTH_KM'] * 1e3:.1f}",
            "source_type": stype,
            "time_function_type": "1",
            "f0": f"{f0:g}",
            "tshift": "0.0",
            "anglesource": angle,
            "Mxx": "1.", "Mzz": "1.", "Mxz": "0.",
            "factor": f"{cfg['SOURCE_FACTOR']:g}",
        }, os.path.join(d, "DATA", "SOURCE"), "SOURCE")
        write_run_script(cfg, d)

        summary = f"""Venus run '{cfg['EXPERIMENT']}/{case}'
air density scale       : {scale:g}
minimum period          : {cfg['MIN_PERIOD_S']} s  (Ricker f0 = {f0:g} Hz)
multiplication_factor   : {p['multiplication_factor']}
radial counts           : air {p['nspec_rad_670_surf']}, upper mantle \
{p['nspec_rad_CMB_670']}, lower mantle {p['nspec_rad_doubling_OC_to_CMB']}/\
{p['nspec_rad_ICB_to_doubling_OC']}, core {p['nspec_rad_Cube_ICB']}
air elements            : {p['h_air_lateral'] / 1e3:.2f} km wide x \
{p['h_air_radial'] / 1e3:.2f} km tall (need <= {p['h_max'] / 1e3:.2f} km)
slowest air sound speed : {c_min:.0f} m/s
elements (half mesh)    : {p['nspec']:,}
DT                      : {dt:g} s (max estimate {dt_est:.4f} s; compare \
with 'Maximum suggested time step' in the solver output)
NSTEP                   : {nstep:,}  ({nstep * dt:.0f} s simulated)
stations                : {nsta}
sponge                  : {cfg['SPONGE']} ({cfg['SPONGE_BOTTOM_KM']}-\
{cfg['AIR_TOP_KM']} km, Q_min {cfg['SPONGE_Q_MIN']})
one-way vertical air travel time to the top: {tau_air:.0f} s
  -> echoes off the top return ~{2 * tau_air:.0f} s after an arrival at a
     ground station; analyze before that if the sponge is off
estimated cost          : {lo:,.1f}-{hi:,.1f} core-hours; \
{lo / cfg['NPROC']:.2f}-{hi / cfg['NPROC']:.2f} h on {cfg['NPROC']} cores
"""
        with open(os.path.join(d, "run_summary.txt"), "w") as f:
            f.write(summary)
        print(f"--- {d}\n{summary}")
    shutil.copy(cfg["_CONFIG_PATH"], os.path.join(exp_dir,
                                                   "venus_config_used.py"))


def main(argv=None):
    """
    Command-line entry point.

    :type argv: list
    :param argv: command-line arguments (defaults to sys.argv[1:])
    """
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[1])
    ap.add_argument("config", nargs="?", help="config file, e.g. "
                    "venus_config.py")
    ap.add_argument("--extract-table", metavar="F90",
                    help="convert an old hard-coded define_external_model.f90 "
                         "into data/venus_interior.txt and exit")
    args = ap.parse_args(argv)
    if args.extract_table:
        extract_interior(args.extract_table,
                         os.path.join(HERE, "data", "venus_interior.txt"))
        return
    if not args.config:
        ap.error("give a config file")
    cfg = load_config(args.config)
    cfg["_CONFIG_PATH"] = os.path.abspath(args.config)
    make_runs(cfg)


if __name__ == "__main__":
    main()
