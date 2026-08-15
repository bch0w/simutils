"""
Read a CMTSOLUTION and STATIONS file and for a given station or for all the 
stations in the file, list the source-receiver distance
"""
import sys
from obspy.geodetics import gps2dist_azimuth, kilometers2degrees


src_fid = "CMTSOLUTION"
sta_fid = "STATIONS"
try:
    sta_choice = sys.argv[1]
except IndexError:
    sta_choice = None

with open(src_fid, "r") as f:
    lines = f.readlines()
    src_lat = float(lines[4].split(":")[1])
    src_lon = float(lines[5].split(":")[1])

with open(sta_fid, "r") as f:
    lines = f.readlines()

for line in lines:
    sta, net, rcv_lat, rcv_lon, *_ = line.split()
    if sta_choice and sta != sta_choice:
        continue

    dist_m, *_ = gps2dist_azimuth(src_lat, src_lon, 
                                  float(rcv_lat), float(rcv_lon))

    dist_km = dist_m * 1E-3
    dist_deg = kilometers2degrees(dist_km)
    print(f"{sta} = {dist_km:.2f}km = {dist_deg:.2f}deg")
    




