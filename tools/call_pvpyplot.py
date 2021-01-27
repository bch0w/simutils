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


# SWITCHES
MAKE = 1
WHITESPACE = 1
TILE = 1

VS = 1
VPVS = 1
POI = 0
MU = 0
UPD = 0

# DEPTH SLICES
basepath = "/Users/Chow/Documents/academic/vuw/forest/relative"
scratch = os.path.join(basepath, "scratch")

# import ipdb; ipdb.set_trace()
# from IPython import embed; embed()

os.chdir(basepath)
for dir_ in ["initial_model", "current_model"]:
    if not MAKE:
        continue
    os.chdir(dir_)

    # Vs Model
    if VS:
        fid = os.path.abspath(glob("model_????_vs.vtk")[0])
        assert(fid and os.path.exists(fid)), f"{dir_} model vs does not exist"
        flags = f"-o {scratch} -c"

        callpv(f"{fid} {flags} -f -Z s -b 1000,3000")
        callpv(f"{fid} {flags} -f -Z 1-5,1 -b 1500,3500")
        # callpv(f"{fid} {flags} -Z 6-10,1 -b 1750,4000")
        # callpv(f"{fid} {flags} -Z 11-15,1 -b 2250,4000")
        # callpv(f"{fid} {flags} -Z 16-20,1 -b 2750,4250")
        # callpv(f"{fid} {flags} -Z 21-25,1 -b 3000,4500")
        # callpv(f"{fid} {flags} -Z 26-30,1 -b 3250,4750")
        # callpv(f"{fid} {flags} -Z 31-40,1 -b 3250,5000")
        # callpv(f"{fid} {flags} -Z 41-50,1 -b 3750,5000")
        # callpv(f"{fid} {flags} -xyti")

    # Vp/Vs Ratio
    if VPVS:
        fid = os.path.abspath(glob("ratio_????_vpvs.vtk")[0])
        assert(fid and os.path.exists(fid)), f"{dir_} ratio vpvs does not exist"
        flags = f"-o {scratch} -c"
        callpv(f"{fid} {flags} -xytif -Z s 1-5,1")
        callpv(f"{fid} {flags} -xyti -Z 6-50,1")

    # Shear Modulus
    if MU:
        fid = os.path.abspath(glob("modulus_????_shear.vtk")[0])
        assert(fid and os.path.exists(fid)), f"{dir_} shear modulus not exist"
        callpv(f"{fid} -o {scratch} -ci")

    # Poisson's Ratio
    if POI:
        fid = os.path.abspath(glob("ratio_????_poissons.vtk")[0])
        assert (fid and os.path.exists(fid)), f"{dir_} poisson's does not exist"
        callpv(f"{fid} -o {scratch} -ci")

    os.chdir("..")

if MAKE and UPD:
    os.chdir("update")
    for fid in glob("update_????_*.vtk"):
        callpv(f"{fid} -o {scratch} -xytic")
        callpv(f"{fid} -o {scratch} -Z s 1-5,1 -fc")
        callpv(f"{fid} -o {scratch} -Z 6-50,1 -c")

# Remove whitespace from all figures
if WHITESPACE:
    os.chdir(scratch)
    for dir_ in glob("*"):
        os.chdir(dir_)
        wspace()
        os.chdir("..")

# Z-Tile images
if TILE:
    z_tile_dir = os.path.join(scratch, "z_tiles")
    if not os.path.exists(z_tile_dir):
        os.mkdir(z_tile_dir)

    os.chdir(scratch)
    tile_dirs = ["model_init_vs", "model_0011_vs", "update_0011_vs",
                 "ratio_init_vpvs", "ratio_0011_vpvs", "update_0011_vpvs"]
    os.chdir(tile_dirs[0])
    z_files = sorted(glob("z_*.png"))
    os.chdir("..")
    for zf in z_files:
        tile_files = [os.path.join(_, zf) for _ in tile_dirs]
        tile(tile_files, out=os.path.join(z_tile_dir, f"tile_{zf}"),
             ncol=3, nrow=2)






