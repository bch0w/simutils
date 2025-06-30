"""
From a list of .adj files, cut down a STATIONS file to make STATIONS_ADJOINT
"""
import os
from glob import glob


stations = "./DATA/STATIONS"
adjsrcs = glob("./SEM/Z/*.adj")

adjsrcs = [os.path.basename(_) for _ in adjsrcs]
adjstas = [_.split(".")[0] + "." + _.split(".")[1] for _ in adjsrcs]

with open("./DATA/STATIONS_ADJOINT", "w") as f:
    for line in open(stations, "r").readlines():
        sta, net, *_ = line.split()
        netsta = f"{net}.{sta}"
        if netsta in adjstas:
            f.write(line)


