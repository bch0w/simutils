from glob import glob
from pyasdf import ASDFDataSet
from pyatoa.utils.asdf.clean import clean_dataset

for fid in glob("*.h5"):
    print(fid)
    with ASDFDataSet(fid) as ds:
        for iter_ in ["i08", "i09", "i10", "i11"]:
            clean_dataset(ds=ds, iteration=iter_)
