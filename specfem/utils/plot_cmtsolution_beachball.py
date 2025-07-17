import sys
from obspy.imaging.beachball import beachball

with open(sys.argv[1], "r") as f:
    lines = f.readlines()

for line in lines:
    print(line)

fm = []
for line in lines[7:]:
    key, val = line.strip().split(":")
    fm.append(float(val))

beachball(fm)
