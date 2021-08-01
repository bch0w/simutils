"""
Short utility script to generate an adjoint source from a cut SPECFEM 
ASCII seismogram for use in generating 'banana-doughnut' style sensitivity
kernels using SPECFEM's adjoint simulation

Calculate half duration using a scaling for moment magnitude
"""
import sys
from glob import glob
import numpy as np
from scipy.signal.windows import tukey
import matplotlib.pyplot as plt
from pyatoa import read_sem
from pyatoa.utils.adjoint import traveltime_adjoint_source
from pyatoa.utils.process import stf_convolve



def preprocess(st_original, tmin, tmax, half_duration=None):
    """
    Quick preprocessing to ensure that only the periods of interest included
    """
    st = st_original.copy()

    st.taper(max_percentage=0.1)
    st.detrend("linear")
    st.filter("bandpass", freqmin=1/tmax, freqmax=1/tmin, zerophase=True)
    st.taper(max_percentage=0.1)
    st.detrend("linear")

    if half_duration is not None:
        st = stf_convolve(st, half_duration)

    return st


def plot(tr, fp):
    """
    Simple plot to visualize the waveform cut and the resulting adjoint source
    """
    # Plot up the resulting waveform before time reversing
    plt.plot(tr.times() + float(tr.stats.starttime), tr.data/tr.data.max(), 
             c="k", label=tr.get_id())
    plt.plot(fp[::-1, 0], fp[:, 1]/fp[:, 1].max(), c="r", 
             label=f"Adjoint Source")
    plt.xlabel("time [s]")
    plt.ylabel("amplitude")
    plt.savefig(f"{tr.get_id()}_phase.png", figsize=(20,8), dpi=200)
    plt.close()
    


if __name__ == "__main__":
    # ====================================
    # DEFINE PARAMETERS HERE
    tmin = 6
    tmax = 30
    windows = {"Z": [235, 245],}
    half_duration = .7
    # =====================================

    if len(sys.argv) < 1:
        sys.exit("Seismogram file id required as argument")
    fids = sys.argv[1:]

    for fid in fids:
        # Define naming schema to save the output files
        net, sta, cha, suffix = fid.split(".")
        comp = cha[-1]
        fid_out = fid.replace(suffix, "adj")
       
        st = read_sem(fid)

        # Write out zeroes if component doesn't match list
        if comp not in windows.keys():
            print(f"{cha[-1]} does not match components, writing zeros")
            zeros = True 
            window = None
        else:
            zeros = False
            window = windows[comp]
            st = preprocess(st, tmin, tmax, half_duration)

        tr = st[0]
        fp = traveltime_adjoint_source(tr=tr, time_window=window, zeros=zeros,
                                       reverse=True, save=fid_out)
        if not zeros:
            plot(tr, fp) 
        
