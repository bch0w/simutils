"""
Call PvPyPlot repeatedly to generate figures that will be used to interpret
the velocity model. Just so I don't have to remember which flags every time
"""
import os
from glob import glob
from subprocess import call

# COMMANDS
def tile(fids, out, ncol, nrow):
    for fid in fids:
        assert(os.path.exists(fid)), f"{fid} does not exist for tile"
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
WS = 1  # Remove border whitespace using Imagemagick
TZ = 1  # Tile Z slice images using Imagemagick
TT = 0  # Tile Trench images using Imagemagick
TI = 0  # Tile Interface images using Imagemagick
# MODELS
VS = 0    # S-Wave Velocity
VP = 0    # P-Wave Velocity
PS = 1    # Vp/Vs Ratio
PO = 0    # Poisson's Ratio
MU = 0    # Shear Modulus
UP = 0    # Net Model Update
# SLICES
T = 0     # Trench
X = 0     # X Slices
Y = 0     # Y Slices
Z = 0     # Z Slices
M = 1     # Manual Z slices for Vp and Vs
I = 0     # Interface
# ==============================================================================

s_flags = ""
if T:
    s_flags += "t"
if X:
    s_flags += "x"
if Y:
    s_flags += "y"
if Z:
    s_flags += "z"
if I:
    s_flags += "i"
if s_flags:
    s_flags = f"-{s_flags}"

# DEPTH SLICES
basepath = "/Users/Chow/Documents/academic/vuw/forest/figures/model_comparisons"
scratch = os.path.join(basepath, "scratch")

# import ipdb; ipdb.set_trace()
# from IPython import embed; embed()

os.chdir(basepath)
for dir_ in ["initial_model", "current_model"]:
    if not MK:
        continue
    os.chdir(dir_)

    # Vp Model
    if VP:
        fid = os.path.abspath(glob("model_????_vp.vtk")[0])
        assert(fid and os.path.exists(fid)), f"{dir_} model vs does not exist"
        flags = f"-o {scratch} -lc"
        
        print("making Vp figures")
        callpv(f"{fid} {flags} {s_flags}")
        if T: 
            callpv(f"{fid} {flags} -t -p trench_vp")

        # Manually set colorbounds for depths
        if M:
            callpv(f"{fid} {flags} -Z s -b 2000,5500")
            callpv(f"{fid} {flags} -Z 1-5,1 -b 2500,6000")
            callpv(f"{fid} {flags} -Z 6-10,1 -b 3500,6500")
            callpv(f"{fid} {flags} -Z 11-15,1 -b 4500,6500")
            callpv(f"{fid} {flags} -Z 16-20,1 -b 5000,7500")
            callpv(f"{fid} {flags} -Z 21-25,1 -b 5500,8250")
            callpv(f"{fid} {flags} -Z 26-30,1 -b 6000,8750")
            callpv(f"{fid} {flags} -Z 31-40,1 -b 6000,9000")
            callpv(f"{fid} {flags} -Z 41-50,1 -b 6500,9250")

    # Vs Model
    if VS:
        fid = os.path.abspath(glob("model_????_vs.vtk")[0])
        assert(fid and os.path.exists(fid)), f"{dir_} model vs does not exist"
        flags = f"-o {scratch} -lc"

        print("making Vs figures")
        callpv(f"{fid} {flags} {s_flags}")
        if T: 
            callpv(f"{fid} {flags} -t -p trench_vs")

        # Manually set colorbounds for depths
        if M:
            callpv(f"{fid} {flags} -Z s -b 1000,3000")
            callpv(f"{fid} {flags} -Z 1-5,1 -b 1500,3500")
            callpv(f"{fid} {flags} -Z 6-10,1 -b 1750,4000")
            callpv(f"{fid} {flags} -Z 11-15,1 -b 2250,4000")
            callpv(f"{fid} {flags} -Z 16-20,1 -b 2750,4250")
            callpv(f"{fid} {flags} -Z 21-25,1 -b 3000,4500")
            callpv(f"{fid} {flags} -Z 26-30,1 -b 3250,4750")
            callpv(f"{fid} {flags} -Z 31-40,1 -b 3250,5000")
            callpv(f"{fid} {flags} -Z 41-50,1 -b 3750,5000")

        # Custom set
        # callpv(f"{fid} {flags} -Z s -b 1000,3000")
        # callpv(f"{fid} {flags} -Z 1-5,1 -b 1500,3500")

    # Vp/Vs Ratio
    if PS:
        fid = os.path.abspath(glob("ratio_????_vpvs.vtk")[0])
        assert(fid and os.path.exists(fid)), f"{dir_} ratio vpvs does not exist"

        print("making Vp/Vs figures")
        flags = f"-o {scratch} -lc"
        if M:
            callpv(f"{fid} {flags} -s -Z s 1-5,1")
            # callpv(f"{fid} {flags} -Z 6-50,1")
        if s_flags:
            callpv(f"{fid} {flags} {s_flags}")

    # Shear Modulus on interface
    if MU:
        fid = os.path.abspath(glob("modulus_????_shear.vtk")[0])
        assert(fid and os.path.exists(fid)), f"{dir_} shear modulus not exist"
        print("making shear modulus figures")
        callpv(f"{fid} -o {scratch} -c {s_flags}")

    # Poisson's Ratio on interface
    if PO:
        fid = os.path.abspath(glob("ratio_????_poissons.vtk")[0])
        assert (fid and os.path.exists(fid)), f"{dir_} poisson's does not exist"
        print("making poisson's ratio figures")
        callpv(f"{fid} -o {scratch} -c {s_flags}")

    os.chdir("..")

if MK and UP:
    os.chdir("update")
    for fid in glob("update_????_*.vtk"):
        print(f"making {fid.split('.')[0]} figures")
        callpv(f"{fid} -o {scratch} -l {s_flags}")
        if M:
            callpv(f"{fid} -o {scratch} -Z s 1-5,1 -l")
            callpv(f"{fid} -o {scratch} -Z 6-50,1 -l")

# Remove whitespace from all figures
if WS:
    os.chdir(scratch)
    print("removing whitespace from all figures in scratch dir")
    for dir_ in glob("*"):
        if os.path.isdir(dir_):
            os.chdir(dir_)
            wspace()
            os.chdir("..")


# I-Tile images
if TI:
    os.chdir(scratch)
    print("making interface tiles")
    tile_dirs = ["model_init_vs", "model_0011_vs", "update_0011_vs",
                 "ratio_init_vpvs", "ratio_0011_vpvs", "update_0011_vpvs",
                 "ratio_init_poissons", "ratio_0011_poissons", "update_0011_poissons",
                 "modulus_init_shear", "modulus_0011_shear", "update_0011_shear"]
    tile_files = [os.path.join(_, "interface.png") for _ in tile_dirs]
    tile(tile_files, out="tile_interface.png", ncol=3, nrow=4)

# T-Tile images
if TT:
    t_tile_dir = os.path.join(scratch, "t_tiles")
    if not os.path.exists(t_tile_dir):
        os.mkdir(t_tile_dir)
    os.chdir(scratch)

    print("making trench tiles")
    tile_dirs = ["model_init_vs", "ratio_init_vpvs",
                 "model_0011_vs", "ratio_0011_vpvs", 
                 "update_0011_vs", "update_0011_vpvs"]
    os.chdir(tile_dirs[0])
    t_files = sorted(glob("t_*.png"))
    os.chdir("..")
    for tf in t_files:
        tile_files = [os.path.join(_, tf) for _ in tile_dirs]
        tile(tile_files, out=os.path.join(t_tile_dir, f"tile_{tf}"),
             ncol=2, nrow=3)

# Z-Tile images
if TZ:
    z_tile_dir = os.path.join(scratch, "z_tiles")
    if not os.path.exists(z_tile_dir):
        os.mkdir(z_tile_dir)
    os.chdir(scratch)
    
    print("making z-slice tiles")
    tile_dirs = ["model_init_vs", "model_0011_vs", "update_0011_vs",
                 "ratio_init_vpvs", "ratio_0011_vpvs", "update_0011_vpvs",
                 "model_init_vp", "model_0011_vp", "update_0011_vp",]
    os.chdir(tile_dirs[0])
    z_files = sorted(glob("z_*.png"))
    os.chdir("..")
    for zf in z_files:
        tile_files = [os.path.join(_, zf) for _ in tile_dirs]
        tile(tile_files, out=os.path.join(z_tile_dir, f"tile_{zf}"),
             ncol=3, nrow=3)







