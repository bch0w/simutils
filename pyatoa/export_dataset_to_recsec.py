"""
Take an ASDFDataSet and export the observed and synthetic data to SAC files that
can be read by RecSec for record section plotting
"""
import os
from glob import glob
from pyasdf import ASDFDataSet 


workdir = "./"
output_dir = "./output"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)
fids = glob(os.path.join(workdir, "*.h5"))
print(f"converting {len(fids)} files")
for fid in fids:
    print(fid)
    with ASDFDataSet(fid) as ds:

