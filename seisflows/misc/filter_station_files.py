"""
Ensure that only source receiver pairs that have measurements are used in this
synthetic-synthetic example
"""
import os
from glob import glob
from pyasdf import ASDFDataSet 


dir_out = "/home/chowbr/current/resolution/sfDATA"
stations = "/home/chowbr/current/resolution/sfDATA/STATIONS"
# datasets = "/scale_wlg_nobackup/filesets/nobackup/nesi00263/bchow/tomo/seisflows/observations/forest/birch/scratch/preprocess/datasets/i11_copies"
datasets = "/home/chowbr/current/aspen/scratch/preprocess/datasets"
iteration = "i17"
step_count = "s01"

for dsfid in glob(os.path.join(datasets, "*.h5")):
    event_id = os.path.basename(dsfid).split(".")[0]
    print(event_id)
    with ASDFDataSet(dsfid) as ds:
        windows = ds.auxiliary_data.MisfitWindows[iteration][step_count].list()        

    sta_list = list(set(["_".join(_.split("_")[:2]) for _ in windows]))
    with open(stations, "r") as fr:
        fid_out = os.path.join(dir_out, f"STATIONS_{event_id}")
        with open(fid_out, "w") as fw:
            for line in fr.readlines():
                sta, net, _, _, _, _ = line.split()
                sta_check = f"{net}_{sta}"
                if sta_check in sta_list:
                    print(f"\t{sta_check}")
                    fw.write(line)
                

