"""
Short utility script to generate an adjoint source from a cut SPECFEM 
ASCII seismogram for use in generating 'banana-doughnut' style sensitivity
kernels using SPECFEM's adjoint simulation
"""
import sys
from glob import glob
import numpy as np
from scipy.signal.windows import tukey
import matplotlib.pyplot as plt
from pyatoa import read_sem
from pyatoa.utils.process import stf_convolve



def select_phase(st_original, t_start, t_end):
    """
    Cut out a portion of the waveform using a Tukey window, replace the rest 
    with zeros to isolate a specific phase
    """
    st = st_original.copy()
     
    for tr in st: 
        overlay = np.zeros(tr.stats.npts) 
        samp_start = int((t_start - float(tr.stats.starttime)) / tr.stats.delta)
        samp_end = int((t_end - float(tr.stats.starttime)) / tr.stats.delta)

        window = tukey(samp_end - samp_start, alpha=0) 
        overlay[samp_start:samp_end] = window

        tr.data = np.multiply(tr.data, overlay)

    return st


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


def plot(tr1, tr2):
    """
    Simple plot to visualize the waveform cut and the resulting adjoint source
    """
    # Plot up the resulting waveform before time reversing
    plt.plot(tr1.times() + float(tr.stats.starttime), tr1.data, c="k", 
             label=tr1.get_id())
    plt.plot(tr2.times() + float(tr.stats.starttime), tr2.data, c="r", 
             label=f"{tr2.get_id()}_phase")
    plt.xlabel("time [s]")
    plt.ylabel("displacement [m]")
    plt.savefig(f"{tr1.get_id()}_phase.png")
    plt.close()
    


if __name__ == "__main__":
    # ====================================
    # DEFINE PARAMETERS HERE
    tmin = 5
    tmax = 30
    windows = {"Z": [225, 240],}
    half_duration = 1.74
    # =====================================

    if len(sys.argv) < 1:
        sys.exit("Seismogram file id required as argument")
    fids = sys.argv[1:]

    for fid in fids:
        # Define naming schema to save the output files
        net, sta, cha, suffix = fid.split(".")
        fid_out = fid.replace(suffix, "adj")
       
        st = read_sem(fid)

        # Write out zeroes if component doesn't match list
        if cha[-1] not in windows.keys():
            print(f"{cha[-1]} does not match components, writing zeros")
            tr = st[0]
            fp = np.zeros(len(st[0].data))

        # Else, process and cut the data
        else:
            t_start, t_end = windows[cha[-1]]
            st_proc = preprocess(st, tmin, tmax, half_duration)
            st_phase = select_phase(st_proc, t_start, t_end)
            plot(st_proc[0], st_phase[0]) 
        
            # Adjoint sources need to be time reversed
            tr = st_phase[0]
            fp = tr.data[::-1]

        # Ensure that the times has the correct offset present in the synthetics
        times = tr.times() + float(tr.stats.starttime)
       
        # Write out the data in two-column ASCII for SPECFEM 
        data = np.vstack((times, fp)).T
        np.savetxt(fid_out, data, "%14.7f %20.8E")

