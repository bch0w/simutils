"""
Call PvPyPlot repeatedly to generate figures that will be used to interpret
the velocity model. Just so I don't have to remember which flags every time
"""
import os
from glob import glob
from subprocess import call

# COMMANDS
def tile(fids, out, ncol, nrow, null=False):
    for fid in fids:
        assert(os.path.exists(fid)), f"{fid} does not exist for tile"
    if null:
        fids.insert(null, "null:")
    call([f"montage -tile {ncol}x{nrow} {' '.join(fids)} -geometry +0+0 {out}"],
         shell=True)


def callpv(argstr):
    run = ("pvpython "
           "/Users/Chow/Documents/academic/vuw/packages/simutils/paraview/"
           "pvpyplot.py ")
    call([run + argstr], shell=True)


def wspace():
    for fid in glob("*"):
        call([f"convert {fid} -fuzz 1% -trim +repage {fid}"], shell=True)


# ================================ SWITCHES ====================================
# ACTIONS
MK = 1  # Make the figures using pvpyplot
TZ = 0  # Tile Z slice images using Imagemagick
TT = 0  # Tile Trench images using Imagemagick
TI = 1  # Tile Interface images using Imagemagick
# MODELS
VS = 0    # S-Wave Velocity
VP = 1    # P-Wave Velocity
PS = 0    # Vp/Vs Ratio
PO = 0    # Poisson's Ratio
MU = 0    # Shear Modulus
UP = 0    # Net Model Update
# SLICES
T = 0    # Trench
X = 0     # X Slices
Y = 0     # Y Slices
Z = 0     # Z Slices
M = 0     # Manual Z slices for Vp and Vs
I = 1     # Interface
F = 0     # Add active faults
S = 0     # Add source epicenter locations
R = 0     # Add receiver locations
C = 1     # Add contour lines to figures
K = 0     # Add interface contours lines
# MANUAL
B = " -b 3000,7000 "  # Set colorbar bounds manually
# ==============================================================================

s_flags = ""
if T:
    s_flags += "t"
if X:
    s_flags += "x"
if Y:
    s_flags += "y"
if I:
    s_flags += "i"
if C:
    s_flags += "c"
if S:
    s_flags += "s"
if R:
    s_flags += "r"
if F:
    s_flags += "f"
if s_flags:
    s_flags = f"-{s_flags}"
if B:
    s_flags += B

if Z:
    s_flags += " -Z 5 15 25 "


# Set paths yo
basepath = "/Users/Chow/Documents/academic/vuw/forest/figures/model_comparisons"
scratch = os.path.join(basepath, "output")

# Standard default flags
flags = f"-o {scratch} -l"
if K:
    flags += " --intcont "

os.chdir(basepath)
for dir_ in ["initial_model", "current_model"]:
    if not MK:
        continue
    os.chdir(dir_)
    print(f"{dir_.replace('_', ' ')}")

    # Vp Model
    if VP:
        fid = os.path.abspath(glob("model_????_vp.vtk")[0])
        assert(fid and os.path.exists(fid)), f"{dir_} model vs does not exist"
        
        print("making Vp figures")
        callpv(f"{fid} {flags} {s_flags}")
        if T: 
            callpv(f"{fid} {flags} -tc -p trench_vp")

        # Manually set colorbounds for depths
        if M:
            callpv(f"{fid} {flags} -Z s -b 2000,5200")
            callpv(f"{fid} {flags} -Z 5 -b 3200,6000")
            callpv(f"{fid} {flags} -Z 10 -b 4300,6500")
            callpv(f"{fid} {flags} -Z 15 -b 5000,6500")
            callpv(f"{fid} {flags} -Z 20 -b 5600,7300")
            callpv(f"{fid} {flags} -Z 25 -b 6000,8200")
            callpv(f"{fid} {flags} -Z 30 -b 6100,8700")
            callpv(f"{fid} {flags} -Z 35 -b 6500,9000")
            callpv(f"{fid} {flags} -Z 40 -b 6900,9200")
            callpv(f"{fid} {flags} -Z 45 -b 7250,9200")
            callpv(f"{fid} {flags} -Z 50 -b 7400,9200")

    # Vs Model
    if VS:
        fid = os.path.abspath(glob("model_????_vs.vtk")[0])
        assert(fid and os.path.exists(fid)), f"{dir_} model vs does not exist"

        print("making Vs figures")
        callpv(f"{fid} {flags} {s_flags}")
        if T: 
            callpv(f"{fid} {flags} -tc -p trench_vs")

        # Manually set colorbounds for depths
        if M:
            callpv(f"{fid} {flags} -Z s -b 1000,3100")
            callpv(f"{fid} {flags} -Z 5 -b 1600,3850")
            callpv(f"{fid} {flags} -Z 10 -b 2000,3900")
            callpv(f"{fid} {flags} -Z 15 -b 2500,4000")
            callpv(f"{fid} {flags} -Z 20 -b 3000,4500")
            callpv(f"{fid} {flags} -Z 25 -b 3200,4500")
            callpv(f"{fid} {flags} -Z 30 -b 3200,4850")
            callpv(f"{fid} {flags} -Z 35 -b 3500,5000")
            callpv(f"{fid} {flags} -Z 40 -b 3750,5000")
            callpv(f"{fid} {flags} -Z 45 50 -b 4000,5000")
            

    # Vp/Vs Ratio
    if PS:
        fid = os.path.abspath(glob("ratio_????_vpvs.vtk")[0])
        assert(fid and os.path.exists(fid)), f"{dir_} ratio vpvs does not exist"

        print("making Vp/Vs figures")
        if M:
            callpv(f"{fid} {flags} -Z s 5-50,5")
        if s_flags:
            callpv(f"{fid} {flags} {s_flags}")

    # Shear Modulus on interface
    if MU:
        fid = os.path.abspath(glob("modulus_????_shear.vtk")[0])
        assert(fid and os.path.exists(fid)), f"{dir_} shear modulus not exist"
        print("making shear modulus figures")
        callpv(f"{fid} -o {scratch} -l {s_flags}")

    # Poisson's Ratio on interface
    if PO:
        fid = os.path.abspath(glob("ratio_????_poissons.vtk")[0])
        assert (fid and os.path.exists(fid)), f"{dir_} poisson's does not exist"
        print("making poisson's ratio figures")
        callpv(f"{fid} -o {scratch} -l {s_flags}")

    os.chdir("..")

if MK and UP:
    os.chdir("update")
    for check, val in zip([VS, VP, PS, PO, MU], 
                          ["vs", "vp", "vpvs", "poisson", "shear"]):
        if check:
            for fid in glob(f"update_????_{val}.vtk"):
                print(f"making {fid.split('.')[0]} figures")
                callpv(f"{fid} -o {scratch} -l {s_flags}")
                if M:
                    callpv(f"{fid} -o {scratch} -Z s 5-50,5 -l")


# Set up for tile generation
tile_list = []
if VP:
    tile_list.append(["model_init_vp", "model_0028_vp", "update_0028_vp"])
if VS:
    tile_list.append(["model_init_vs", "model_0028_vs", "update_0028_vs"])
if PS:
    tile_list.append(["ratio_init_vpvs", "ratio_0028_vpvs", 
                      "update_0028_vpvs"])
if MU:
    tile_list.append(["modulus_init_shear", "modulus_0028_shear", 
                      "update_0028_shear"])
if PO:
    tile_list.append(["ratio_init_poissons", "ratio_0028_poissons", 
                      "update_0028_poissons"])

# I-Tile images
if TI:
    i_tile_dir = os.path.join(scratch, "i_tiles")
    if not os.path.exists(i_tile_dir):
        os.mkdir(i_tile_dir)
    os.chdir(scratch)

    print("making interface tiles")
    for tile_dirs in tile_list:
        tile_files = [os.path.join(_, "interface.png") for _ in tile_dirs]
        name = tile_dirs[0].replace("init_", "")
        tile(tile_files, out=os.path.join(i_tile_dir, f"tile_{name}.png"), 
             ncol=3, nrow=1)

# T-Tile images
if TT:
    t_tile_dir = os.path.join(scratch, "t_tiles")
    if not os.path.exists(t_tile_dir):
        os.mkdir(t_tile_dir)
    os.chdir(scratch)

    print("making trench tiles")
    tile_dirs = ["model_init_vs", "model_0028_vs", "update_0028_vs"]
    # tile_dirs = ["ratio_init_vpvs", "ratio_0028_vpvs", "update_0028_vpvs"]
    os.chdir(tile_dirs[0])
    t_files = sorted(glob("t_*.png"))
    os.chdir("..")
    for tf in t_files:
        tile_files = [os.path.join(_, tf) for _ in tile_dirs]
        tile(tile_files, out=os.path.join(t_tile_dir, f"tile_{tf}"),
             ncol=2, nrow=2, null=4)

# Z-Tile images
if TZ:
    z_tile_dir = os.path.join(scratch, "z_tiles")
    if not os.path.exists(z_tile_dir):
        os.mkdir(z_tile_dir)
    os.chdir(scratch)
    
    print("making z-slice tiles")
    for tile_dirs in tile_list:
        qty = tile_dirs[0].split("_")[-1]
        os.chdir(tile_dirs[0])
        z_files = sorted(glob("z_*.png"))
        os.chdir("..")
        for zf in z_files:
            tile_files = [os.path.join(_, zf) for _ in tile_dirs]
            tile(tile_files, out=os.path.join(z_tile_dir, f"tile_{qty}_{zf}"),
                 ncol=len(tile_dirs), nrow=1)







