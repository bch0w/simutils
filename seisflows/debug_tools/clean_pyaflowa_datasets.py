import os
from glob import glob
from pyasdf import ASDFDataSet
from pyatoa.utils.asdf.clean import clean_dataset

for fid in glob("*.h5"):
    print(fid)
    with ASDFDataSet(fid) as ds:
        for iter_ in [4]:
            iteration = f"i{iter_:0>2}"
            for step in [1, 2, 3, 4]:
                step_count = f"s{step:0>2}"
                print(f"{iteration}{step_count}")
                clean_dataset(ds, iteration=iteration, step_count=step_count)
