"""
Rough draft attempt to implement Ruan et al. (2019) source receiver weighting
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from obspy.geodetics import gps2dist_azimuth
from pysep.utils.io import read_stations


def plot_weights(x, y, c, title, save):
    plt.scatter(x, y, c=c, cmap="jet_r", marker="^", # vmin=0, vmax=10,
                norm=mpl.colors.LogNorm())
    plt.colorbar()
    plt.xlabel("lon")
    plt.ylabel("lat")
    plt.title(title)
    plt.savefig(save)
    plt.close("all")

# First get station geographic information
inv = read_stations("/home/bchow/REPOSITORIES/spectral/nzatom/STATIONS")
codes, lons, lats = [], [], []
for net in inv:
    for sta in net:
        codes.append(f"{net.code}_{sta.code}")
        lons.append(sta.longitude)
        lats.append(sta.latitude)

# Calculate receiver-receiver distances
dists = []
_ref_dists = []
for code_i, lon_i, lat_i in zip(codes, lons, lats):
    dists_i = []
    for code_j, lon_j, lat_j in zip(codes, lons, lats):
        # Don't consider w_ii
        if code_i == code_j:
            continue

        dist_m, *_ = gps2dist_azimuth(lat1=lat_i, lon1=lon_i, 
                                      lat2=lat_j, lon2=lon_j)
        dists_i.append(dist_m * 1E-3)  # units: km
        _ref_dists.append(dist_m * 1E-3)
    dists.append(dists_i)

# Determine minimum/maximum reference distancmin_ref_dist = np.min(dists)
min_ref_dist = np.floor(min(_ref_dists)) or 1
max_ref_dist = max(_ref_dists) 

min_ref_dist = 10
max_ref_dist = 100

# Loop through reference distances and calculate weightsS
ratios = []
ref_dists = np.arange(min_ref_dist, max_ref_dist, 1)

for ref_dist in ref_dists:
    weights = []
    for dist_i in dists:
        weight_i = 0
        for dist_ in dist_i:
            weight_i += np.e ** (-1 * (dist_ / ref_dist)**2)
        weights.append(1 / weight_i)
    # Plot weights
    plot_weights(x=lons, y=lats, c=weights, title=f"ref dist={ref_dist}",
                 save=f"./figures/weight_{ref_dist}.png")

    # Calculate ratio of min and max weights
    ratios.append(max(weights) / min(weights))

plt.scatter(ref_dists, ratios)
plt.ylim([0, 100])
plt.xlabel("Reference distance")
plt.ylabel("Condition number")
plt.show()

    

                    



