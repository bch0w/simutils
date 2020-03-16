import numpy as np
from glob import glob


def split_values(array):
    """split loaded array into key and value"""
    key, value = array
    key.replace(" ", "_")
    value = float(value.strip())
    return key, value

cmtsolutions = glob("CMTSOLUTION*")
dict_out = {key: [] for key in 
              ["latitude", "longitude", "depth", "Mrr", "Mtt", 
               "Mpp", "Mrt", "Mrp", "Mtp"]}

for cmtfid in cmtsolutions:
    cmt = np.loadtxt(cmtfid, skiprows=1, dtype="str", delimiter=":")
    for idx in range(3, len(cmt)): 
        key, value = split_values(cmt[idx])
        dict_out[key].append(value)

np.savez("all_events", **dict_out)
    
