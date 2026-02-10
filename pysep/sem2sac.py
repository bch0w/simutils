import sys
from glob import glob
from pathlib import Path
from obspy import UTCDateTime
from pysep import read_sem, read_events_plus, read_stations
from pysep.utils.cap_sac import append_sac_headers

starttime = UTCDateTime("2017-09-03T03:30:01.769538Z")
endtime = UTCDateTime("2017-09-03T03:37:28.940951Z")
sampling_rate = 20

files = Path(".")
for fid in files.glob("*.sem?"):
    st = read_sem(fid, source="./CMTSOLUTION", stations="./STATIONS")
    st.trim(starttime, endtime)
    st.resample(sampling_rate)

    fid_out = str(fid).replace("semv", "sac")
    st.write(f"../{fid_out}", format="SAC")
