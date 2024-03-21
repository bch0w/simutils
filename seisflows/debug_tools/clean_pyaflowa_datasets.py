import os
from glob import glob
from pyasdf import ASDFDataSet
from pyatoa.utils.asdf.clean import clean_dataset

for fid in glob("*.h5"):
    print(fid)
    with ASDFDataSet(fid) as ds:
        for step_count in ["s01", "s02", "s03"]:
            clean_dataset(ds, iteration="i01", step_count=step_count)
