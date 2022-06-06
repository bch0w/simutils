"""
Quick utility to convert units of mseed data stored in SEED format
"""
import os
from glob import glob
from obspy import read


path = "/scale_wlg_persistent/filesets/project/gns03247/bchow/data/waveforms"

# MAke sure we don't overwrite
files_changed = open("files_changed.txt", "r").readlines()
files_changed = [_.strip() for _ in files_changed]

with open("files_changed.txt", "w") as f:
    for dep in glob(os.path.join(path, "*")):
        if os.path.basename(dep) not in ["GA", "HD"]:
            continue
        for year in glob(os.path.join(dep, "*")):
            for net in glob(os.path.join(year, "Z?")):
                if os.path.basename(net) not in ["ZX", "Z8"]:
                    continue    
                for sta in glob(os.path.join(net, "*")):
                    for comp in glob(os.path.join(sta, "*")):
                        for fid in glob(os.path.join(comp, "*")):
                            if fid in files_changed:
                                print(f"SKIP: {fid}")
                                continue
                            st = read(fid)
                            for tr in st:
                                tr.data *= 1E-9
                            f.write(f"{fid}\n")
                            st.write(fid, format="MSEED")
                            print(f"CHANGE: {fid}")
                    
