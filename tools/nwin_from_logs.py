"""
Collect the summary log information that is outputted by Pyaflowa during a 
SeisFlows Inversion
"""
from glob import glob

iterstep = "i06s00"

nwin = 0
for fid in glob(f"{iterstep}_*.log"):
    with open(fid, "r") as f:
        line = f.readlines()[-3]
        nwin += int(line.split(" ")[-1]) 
        # lines = f.readlines()[-5:]
        # source_name = lines[0].strip().split(" ")[-1]
        # return_dict[source_name] = {}
        # for line in lines[1:]:
        #     key, val = line.strip().split(":")
        #     return_dict[source_name][key.strip()] = val.strip()

print(f"{nwin} windows") 

