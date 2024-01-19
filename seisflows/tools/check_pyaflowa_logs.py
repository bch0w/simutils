import os
import sys
from glob import glob

evaluation = sys.argv[1]
globstr = f"*_{evaluation}.log"
nwin = 0
misfit = 0
print(globstr)
assert(glob(globstr))
for fid in glob(globstr):
    with open(fid) as f:
        lines = f.readlines()
        line = lines[5]
        assert("WINDOWS" in line)
        nval = int(line.strip().split(":")[1])
        nwin += nval

        line = lines[6]
        assert("MISFIT" in line)
        misfit_ = float(line.strip().split(":")[1])
        misfit += misfit_ / 2


print(f"nwin={nwin}; misfit={misfit}; scaled={misfit/nwin}")
